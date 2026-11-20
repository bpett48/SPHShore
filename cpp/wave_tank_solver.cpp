//This is a weakly compressible 3D SPH wave tank by Ben Pett
//This version has a moving piston, a sloping shoreline and two fixed cylindrical pillars
//Fixed: floor-first pillar collision projection, including particles that overshoot below the seabed
//Optional file output, proportional geometry, and unused grid propulsion/velocity-map helpers
//One pillar is submerged and the other rises above the initial water level
//The fluid uses a 3D cubic spline kernel, density diffusion and artificial viscosity
//Only fluid particles are saved to LogBook.txt; the geometry is saved separately

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <chrono>
#include <random>
#include <string>
#include <limits>
#include <cstdint>
#include <memory>
#include <cstring>
using namespace std;

struct Vector
{
    double X = 0.0;
    double Y = 0.0;
    double Z = 0.0;
};

struct Pillar
{
    //A vertical cylinder fixed to the flat floor
    double X = 0.0;
    double Y = 0.0;
    double Radius = 0.04;
    double Height = 0.10;
};

struct parameters
{
    //These settings control the detail and duration of the simulation
    bool WriteFiles = true; //False: do not create or open any output files
    int N = 12000;
    int Iterations = 50000;
    double TimeStep = 0.00010;
    double SaveEveryTime = 0.025;
    double DiagnosticsEveryTime = 0.5;
    double MaxRunHours = 2.8;
    unsigned int RandomSeed = 42069;

    //Change MaxDimmensions to scale the tank and its attached geometry.
    //The reference geometry corresponds to {1.20, 0.45, 0.30} metres.
    //Uniform scaling of all three components preserves exact shape similarity.
    //Non-uniform scaling preserves relative positions along each tank axis.
    //Call ConfigureGeometry after changing MinDimmensions or MaxDimmensions.
    //A longer flume provides space for waves to develop around the pillars
    Vector MinDimmensions = { 0.0, 0.0, 0.0 };
    Vector MaxDimmensions = { 1.20, 0.45, 0.30 };
    double WaterLevel = 0.13;

    //Particle properties are calculated automatically from N and the available water volume
    double ParticleSpacing = 0.0;
    double SmoothingLength = 0.0;
    double ParticleMass = 0.0;
    double SmoothingLengthFactor = 1.40;

    //Fluid properties
    double ReferenceDensity = 1000.0;
    double SoundSpeed = 25.0;
    double Gamma = 7.0;
    double GravAccelaration = -9.81;
    double Pi = 3.141592653589793;

    //Numerical stability
    double DensityDiffusion = 0.10;
    double ArtificialViscosity = 0.055;
    double XSPHFactor = 0.05;
    double WallRestitution = 0.0;
    double InitialJitter = 0.025;

    //The sloping shoreline begins after both pillars
    double RampStartX = 0.64;
    double RampAngle = 16.0;

    //The piston creates continuous waves at the left end
    double PistonBaseX = 0.03;
    double PistonAmplitude = 0.007;
    double WaveFrequency = 1.6;
    double WaveStartTime = 0.25;
    double WaveRampTime = 0.9;
    bool UseCustomPiston = false;
    double CustomPistonPosition = 0.0;
    double CustomPistonVelocity = 0.0;

    //Two pillars are staggered across the width so their wakes are visible
    Pillar UnderwaterPillar = { 0.30, 0.145, 0.045, 0.075 };
    Pillar AboveWaterPillar = { 0.48, 0.305, 0.045, 0.190 };
};


//Update every dimensional geometry value from the tank extents.
//This function is called once before validation and particle initialization.
//X and Y coordinates scale with tank length and width respectively; Z
//coordinates scale with tank height. Cylinder radii scale with the smaller
//horizontal factor so they stay circular and inside the flume.
//For exact geometrical similarity, multiply all three MaxDimmensions by
//the SAME factor relative to the reference {1.20,0.45,0.30}.
void ConfigureGeometry(parameters& Par)
{
    const double SX = (Par.MaxDimmensions.X - Par.MinDimmensions.X) / 1.20;
    const double SY = (Par.MaxDimmensions.Y - Par.MinDimmensions.Y) / 0.45;
    const double SZ = (Par.MaxDimmensions.Z - Par.MinDimmensions.Z) / 0.30;
    if (!isfinite(SX) || !isfinite(SY) || !isfinite(SZ) ||
        SX <= 0.0 || SY <= 0.0 || SZ <= 0.0)
    {
        throw runtime_error("Tank dimensions must be positive and finite");
    }
    const double SR = min(SX, SY);
    Par.WaterLevel = Par.MinDimmensions.Z + 0.13 * SZ;
    Par.RampStartX = Par.MinDimmensions.X + 0.64 * SX;
    //Preserve the physical slope shape when X and Z are scaled differently.
    Par.RampAngle = atan(tan(16.0 * Par.Pi / 180.0) * SZ / SX) * 180.0 / Par.Pi;
    Par.PistonBaseX = Par.MinDimmensions.X + 0.03 * SX;
    Par.PistonAmplitude = 0.007 * SX;
    Par.UnderwaterPillar = {
        Par.MinDimmensions.X + 0.30 * SX,
        Par.MinDimmensions.Y + 0.145 * SY,
        0.045 * SR, 0.075 * SZ
    };
    Par.AboveWaterPillar = {
        Par.MinDimmensions.X + 0.48 * SX,
        Par.MinDimmensions.Y + 0.305 * SY,
        0.045 * SR, 0.190 * SZ
    };
}

Vector DisplacmentOfVectors(Vector A, Vector B)
{
    return { A.X - B.X, A.Y - B.Y, A.Z - B.Z };
}

double GetRampHeight(double X, const parameters& Par)
{
    if (X <= Par.RampStartX) { return Par.MinDimmensions.Z; }
    return Par.MinDimmensions.Z + (X - Par.RampStartX) * tan(Par.RampAngle * Par.Pi / 180.0);
}

Vector GetRampNormal(const parameters& Par)
{
    double Angle = Par.RampAngle * Par.Pi / 180.0;
    return { -sin(Angle), 0.0, cos(Angle) };
}

double GetRampDistance(Vector Position, const parameters& Par)
{
    Vector Normal = GetRampNormal(Par);
    return (Position.X - Par.RampStartX) * Normal.X +
        (Position.Z - Par.MinDimmensions.Z) * Normal.Z;
}

double GetWaveEnvelope(double Time, const parameters& Par)
{
    if (Time <= Par.WaveStartTime) { return 0.0; }
    if (Par.WaveRampTime <= 0.0) { return 1.0; }
    double Q = min(1.0, (Time - Par.WaveStartTime) / Par.WaveRampTime);
    return Q * Q * (3.0 - 2.0 * Q);
}

double GetWaveEnvelopeRate(double Time, const parameters& Par)
{
    if (Time <= Par.WaveStartTime || Par.WaveRampTime <= 0.0 ||
        Time >= Par.WaveStartTime + Par.WaveRampTime) { return 0.0; }
    double Q = (Time - Par.WaveStartTime) / Par.WaveRampTime;
    return 6.0 * Q * (1.0 - Q) / Par.WaveRampTime;
}

double GetPistonPosition(double Time, const parameters& Par)
{
    if (Par.UseCustomPiston) return Par.CustomPistonPosition;
    double Phase = 2.0 * Par.Pi * Par.WaveFrequency * (Time - Par.WaveStartTime);
    return Par.PistonBaseX + Par.PistonAmplitude * GetWaveEnvelope(Time, Par) * sin(Phase);
}

double GetPistonVelocity(double Time, const parameters& Par)
{
    if (Par.UseCustomPiston) return Par.CustomPistonVelocity;
    double Phase = 2.0 * Par.Pi * Par.WaveFrequency * (Time - Par.WaveStartTime);
    return Par.PistonAmplitude * (GetWaveEnvelopeRate(Time, Par) * sin(Phase) +
        GetWaveEnvelope(Time, Par) * 2.0 * Par.Pi * Par.WaveFrequency * cos(Phase));
}

bool IsInsidePillar(Vector Position, const Pillar& Obstacle, const parameters& Par)
{
    double DX = Position.X - Obstacle.X;
    double DY = Position.Y - Obstacle.Y;
    return DX * DX + DY * DY < Obstacle.Radius * Obstacle.Radius &&
        Position.Z >= Par.MinDimmensions.Z &&
        Position.Z < Par.MinDimmensions.Z + Obstacle.Height;
}

bool IsValidWaterPosition(Vector Position, const parameters& Par)
{
    return Position.X > Par.PistonBaseX && Position.X < Par.MaxDimmensions.X &&
        Position.Y > Par.MinDimmensions.Y && Position.Y < Par.MaxDimmensions.Y &&
        Position.Z > GetRampHeight(Position.X, Par) && Position.Z < Par.WaterLevel &&
        !IsInsidePillar(Position, Par.UnderwaterPillar, Par) &&
        !IsInsidePillar(Position, Par.AboveWaterPillar, Par);
}

double GetPressure(double Density, const parameters& Par)
{
    //Tait equation of state; negative pressure is disabled
    double B = Par.ReferenceDensity * Par.SoundSpeed * Par.SoundSpeed / Par.Gamma;
    return max(0.0, B * (pow(Density / Par.ReferenceDensity, Par.Gamma) - 1.0));
}

double GetDensityFromPressure(double Pressure, const parameters& Par)
{
    double B = Par.ReferenceDensity * Par.SoundSpeed * Par.SoundSpeed / Par.Gamma;
    return Par.ReferenceDensity * pow(max(0.01, 1.0 + Pressure / B), 1.0 / Par.Gamma);
}

double CSKernel(double Distance, const parameters& Par)
{
    double Q = Distance / Par.SmoothingLength;
    if (Q >= 2.0) { return 0.0; }
    double H = Par.SmoothingLength;
    double FC = 1.0 / (Par.Pi * H * H * H);
    if (Q < 1.0) { return FC * (1.0 - 1.5 * Q * Q + 0.75 * Q * Q * Q); }
    double T = 2.0 - Q;
    return FC * T * T * T / 4.0;
}

Vector CSKernelGradient(Vector D, double Distance, const parameters& Par)
{
    if (Distance < 1e-12 || Distance >= 2.0 * Par.SmoothingLength) { return {}; }
    double H = Par.SmoothingLength;
    double Q = Distance / H;
    double FC = 1.0 / (Par.Pi * H * H * H * H);
    double dWdr = Q < 1.0 ? FC * (-3.0 * Q + 2.25 * Q * Q) :
        FC * (-0.75 * (2.0 - Q) * (2.0 - Q));
    double Scale = dWdr / Distance;
    return { D.X * Scale, D.Y * Scale, D.Z * Scale };
}

struct Particle
{
    int ID = -1;
    Vector Position;
    Vector Velocity;
    Vector Acceleration;
    Vector XSPH;
    double Density = 1000.0;
    double Pressure = 0.0;
    double DensityRate = 0.0;
    bool IsBoundary = false;

    Particle() = default;
    Particle(int StartID, Vector P, Vector V, double Rho, bool Boundary = false)
    {
        ID = StartID;
        Position = P;
        Velocity = V;
        Density = Rho;
        IsBoundary = Boundary;
    }
};

struct Cell
{
    int X;
    int Y;
    int Z;
};

struct SpatialGrid
{
    double CellSize;
    Vector Origin;
    int NX;
    int NY;
    int NZ;
    vector<vector<int>> Cells;

    SpatialGrid(const parameters& Par)
    {
        CellSize = 2.0 * Par.SmoothingLength;
        Origin = { Par.MinDimmensions.X - 3.0 * CellSize,
                   Par.MinDimmensions.Y - 3.0 * CellSize,
                   Par.MinDimmensions.Z - 3.0 * CellSize };
        NX = (int)ceil((Par.MaxDimmensions.X - Origin.X) / CellSize) + 4;
        NY = (int)ceil((Par.MaxDimmensions.Y - Origin.Y) / CellSize) + 4;
        NZ = (int)ceil((Par.MaxDimmensions.Z - Origin.Z) / CellSize) + 4;
        Cells.resize((size_t)NX * NY * NZ);
    }

    Cell GetCell(Vector P) const
    {
        return { (int)floor((P.X - Origin.X) / CellSize),
                 (int)floor((P.Y - Origin.Y) / CellSize),
                 (int)floor((P.Z - Origin.Z) / CellSize) };
    }

    bool IsInside(Cell C) const
    {
        return C.X >= 0 && C.X < NX && C.Y >= 0 && C.Y < NY && C.Z >= 0 && C.Z < NZ;
    }

    size_t GetIndex(Cell C) const
    {
        return ((size_t)C.Z * NY + C.Y) * NX + C.X;
    }

    void Build(const vector<Particle>& Population)
    {
        for (auto& CurrentCell : Cells) { CurrentCell.clear(); }
        for (int i = 0; i < (int)Population.size(); i++)
        {
            Cell C = GetCell(Population[i].Position);
            if (!IsInside(C)) { throw runtime_error("Particle outside spatial grid"); }
            Cells[GetIndex(C)].push_back(i);
        }
    }
};


//A velocity vector sampled at the centre of an occupied spatial-hash cell.
struct MapVector
{
    Vector Position;
    Vector Velocity;
    int ParticleCount = 0;
};

//UNUSED HOOK: apply a total force (newtons) to the real fluid particles
//currently located in the hashing cell containing Position.
//The same total force is shared equally across all real particles in that
//cell. It adds acceleration = Force / (ParticleMass * Count) to each one.
//Call only when Grid.Build(AllParticles) has already been performed for
//the CURRENT positions. Boundary/mirror particles are deliberately ignored.
//Returns the number of real particles affected; empty cells return zero.
int AddPropulsion(Vector Position, Vector Force, vector<Particle>& Population,
    const SpatialGrid& Grid, const parameters& Par)
{
    const Cell C = Grid.GetCell(Position);
    if (!Grid.IsInside(C)) { return 0; }
    if (!isfinite(Force.X) || !isfinite(Force.Y) || !isfinite(Force.Z) ||
        Par.ParticleMass <= 0.0) { throw runtime_error("Invalid propulsion force or mass"); }
    const auto& Indices = Grid.Cells[Grid.GetIndex(C)];
    int Count = 0;
    for (int I : Indices)
    {
        if (I >= 0 && I < (int)Population.size() && !Population[I].IsBoundary)
            ++Count;
    }
    if (Count == 0) { return 0; }
    const double Factor = 1.0 / (Par.ParticleMass * Count);
    for (int I : Indices)
    {
        if (I < 0 || I >= (int)Population.size() || Population[I].IsBoundary) { continue; }
        Population[I].Acceleration.X += Force.X * Factor;
        Population[I].Acceleration.Y += Force.Y * Factor;
        Population[I].Acceleration.Z += Force.Z * Factor;
    }
    return Count;
}

//UNUSED HOOK: outputMap constructs a velocity-vector map in MEMORY, without
//writing any file. Each occupied grid cell contributes its centre coordinate
//and the average velocity of its real fluid particles.
//Boundary particles are excluded even when the grid contains mirrors.
//To write a map to disk in a future version, serialize its returned vector
//only inside a block guarded by Par.WriteFiles.
vector<MapVector> outputMap(const vector<Particle>& Population,
    const SpatialGrid& Grid)
{
    vector<MapVector> Result;
    for (int Z = 0; Z < Grid.NZ; ++Z)
    for (int Y = 0; Y < Grid.NY; ++Y)
    for (int X = 0; X < Grid.NX; ++X)
    {
        Cell C{ X, Y, Z };
        const auto& Indices = Grid.Cells[Grid.GetIndex(C)];
        if (Indices.empty()) { continue; }
        MapVector Sample;
        for (int I : Indices)
        {
            if (I < 0 || I >= (int)Population.size() || Population[I].IsBoundary)
                continue;
            const Vector V = Population[I].Velocity;
            Sample.Velocity.X += V.X;
            Sample.Velocity.Y += V.Y;
            Sample.Velocity.Z += V.Z;
            ++Sample.ParticleCount;
        }
        if (Sample.ParticleCount == 0) { continue; }
        const double Inv = 1.0 / Sample.ParticleCount;
        Sample.Velocity.X *= Inv;
        Sample.Velocity.Y *= Inv;
        Sample.Velocity.Z *= Inv;
        Sample.Position = {
            Grid.Origin.X + (X + 0.5) * Grid.CellSize,
            Grid.Origin.Y + (Y + 0.5) * Grid.CellSize,
            Grid.Origin.Z + (Z + 0.5) * Grid.CellSize
        };
        Result.push_back(Sample);
    }
    return Result;
}

void AddMirrorParticle(const Particle& P, vector<Particle>& BoundaryParticles,
    Vector MirrorPosition, Vector MirrorVelocity, const parameters& Par)
{
    //Estimate virtual density and pressure by extrapolating gravity through the wall
    double MirrorPressure = P.Pressure + Par.ReferenceDensity * Par.GravAccelaration *
        (MirrorPosition.Z - P.Position.Z);
    MirrorPressure = max(0.0, MirrorPressure);
    BoundaryParticles.emplace_back(-1, MirrorPosition, MirrorVelocity,
        GetDensityFromPressure(MirrorPressure, Par), true);
    BoundaryParticles.back().Pressure = MirrorPressure;
}

void AddPillarBoundaryParticles(const Particle& P, const Pillar& Obstacle,
    vector<Particle>& BoundaryParticles, const parameters& Par)
{
    //Mirror nearby water particles across the curved cylinder side or its top cap
    double Support = 2.0 * Par.SmoothingLength;
    double DX = P.Position.X - Obstacle.X;
    double DY = P.Position.Y - Obstacle.Y;
    double R = sqrt(DX * DX + DY * DY);
    double Top = Par.MinDimmensions.Z + Obstacle.Height;

    //The curved side is a solid boundary up to the top of the pillar
    if (P.Position.Z <= Top && R >= Obstacle.Radius &&
        R - Obstacle.Radius < Support && R > 1e-12)
    {
        double Distance = R - Obstacle.Radius;
        Vector Normal = { DX / R, DY / R, 0.0 };
        Vector Mirror = { P.Position.X - 2.0 * Distance * Normal.X,
                          P.Position.Y - 2.0 * Distance * Normal.Y,
                          P.Position.Z };
        double VN = P.Velocity.X * Normal.X + P.Velocity.Y * Normal.Y;
        Vector MirrorVelocity = { P.Velocity.X - 2.0 * VN * Normal.X,
                                  P.Velocity.Y - 2.0 * VN * Normal.Y,
                                  P.Velocity.Z };
        AddMirrorParticle(P, BoundaryParticles, Mirror, MirrorVelocity, Par);
    }

    //The underwater pillar needs a horizontal top surface as well
    if (R <= Obstacle.Radius && P.Position.Z >= Top &&
        P.Position.Z - Top < Support)
    {
        Vector Mirror = { P.Position.X, P.Position.Y, 2.0 * Top - P.Position.Z };
        Vector MirrorVelocity = { P.Velocity.X, P.Velocity.Y, -P.Velocity.Z };
        AddMirrorParticle(P, BoundaryParticles, Mirror, MirrorVelocity, Par);
    }
}

void CreateBoundaryParticles(const vector<Particle>& Population,
    vector<Particle>& BoundaryParticles, double Time, const parameters& Par)
{
    //Create virtual SPH particles near the piston, walls, ramp and both pillars
    BoundaryParticles.clear();
    double Support = 2.0 * Par.SmoothingLength;
    double Left = GetPistonPosition(Time, Par);
    double LeftVelocity = GetPistonVelocity(Time, Par);
    Vector Normal = GetRampNormal(Par);

    for (const Particle& P : Population)
    {
        Vector X = P.Position;
        Vector V = P.Velocity;

        if (X.X - Left < Support)
        {
            AddMirrorParticle(P, BoundaryParticles, { 2.0 * Left - X.X, X.Y, X.Z },
                { 2.0 * LeftVelocity - V.X, V.Y, V.Z }, Par);
        }
        if (Par.MaxDimmensions.X - X.X < Support)
        {
            AddMirrorParticle(P, BoundaryParticles,
                { 2.0 * Par.MaxDimmensions.X - X.X, X.Y, X.Z }, { -V.X, V.Y, V.Z }, Par);
        }
        if (X.Y - Par.MinDimmensions.Y < Support)
        {
            AddMirrorParticle(P, BoundaryParticles,
                { X.X, 2.0 * Par.MinDimmensions.Y - X.Y, X.Z }, { V.X, -V.Y, V.Z }, Par);
        }
        if (Par.MaxDimmensions.Y - X.Y < Support)
        {
            AddMirrorParticle(P, BoundaryParticles,
                { X.X, 2.0 * Par.MaxDimmensions.Y - X.Y, X.Z }, { V.X, -V.Y, V.Z }, Par);
        }
        if (X.Z - Par.MinDimmensions.Z < Support && X.X <= Par.RampStartX)
        {
            AddMirrorParticle(P, BoundaryParticles,
                { X.X, X.Y, 2.0 * Par.MinDimmensions.Z - X.Z }, { V.X, V.Y, -V.Z }, Par);
        }
        if (Par.MaxDimmensions.Z - X.Z < Support)
        {
            AddMirrorParticle(P, BoundaryParticles,
                { X.X, X.Y, 2.0 * Par.MaxDimmensions.Z - X.Z }, { V.X, V.Y, -V.Z }, Par);
        }

        double Distance = GetRampDistance(X, Par);
        double ContactX = X.X - Distance * Normal.X;
        if (Distance >= 0.0 && Distance < Support &&
            ContactX >= Par.RampStartX && ContactX <= Par.MaxDimmensions.X)
        {
            Vector Mirror = { X.X - 2.0 * Distance * Normal.X, X.Y,
                              X.Z - 2.0 * Distance * Normal.Z };
            double VN = V.X * Normal.X + V.Z * Normal.Z;
            Vector MirrorVelocity = { V.X - 2.0 * VN * Normal.X, V.Y,
                                      V.Z - 2.0 * VN * Normal.Z };
            AddMirrorParticle(P, BoundaryParticles, Mirror, MirrorVelocity, Par);
        }

        AddPillarBoundaryParticles(P, Par.UnderwaterPillar, BoundaryParticles, Par);
        AddPillarBoundaryParticles(P, Par.AboveWaterPillar, BoundaryParticles, Par);
    }
}

void CalculateTotalForce(vector<Particle>& Population,
    const vector<Particle>& BoundaryParticles, const SpatialGrid& Grid,
    const parameters& Par)
{
    //Calculate pressure, artificial viscosity, continuity and XSPH smoothing
    int ParticleCount = (int)Population.size();
    double H2 = Par.SmoothingLength * Par.SmoothingLength;

#pragma omp parallel for schedule(static)
    for (int i = 0; i < ParticleCount; i++)
    {
        Particle& P = Population[i];
        P.Acceleration = { 0.0, 0.0, Par.GravAccelaration };
        P.DensityRate = 0.0;
        P.XSPH = {};
        Cell C = Grid.GetCell(P.Position);

        for (int X = -1; X <= 1; X++)
        for (int Y = -1; Y <= 1; Y++)
        for (int Z = -1; Z <= 1; Z++)
        {
            Cell SearchCell = { C.X + X, C.Y + Y, C.Z + Z };
            if (!Grid.IsInside(SearchCell)) { continue; }
            for (int Index : Grid.Cells[Grid.GetIndex(SearchCell)])
            {
                if (Index == i) { continue; }
                const Particle& Neighbor = Index < ParticleCount ?
                    Population[Index] : BoundaryParticles[Index - ParticleCount];
                Vector R = DisplacmentOfVectors(P.Position, Neighbor.Position);
                double R2 = R.X * R.X + R.Y * R.Y + R.Z * R.Z;
                if (R2 < 1e-24 || R2 >= 4.0 * H2) { continue; }

                double Distance = sqrt(R2);
                Vector Grad = CSKernelGradient(R, Distance, Par);
                Vector DV = DisplacmentOfVectors(P.Velocity, Neighbor.Velocity);
                double DotVG = DV.X * Grad.X + DV.Y * Grad.Y + DV.Z * Grad.Z;
                double DotVR = DV.X * R.X + DV.Y * R.Y + DV.Z * R.Z;

                //The continuity equation and density diffusion
                P.DensityRate += Par.ParticleMass * DotVG;
                if (!Neighbor.IsBoundary)
                {
                    double DotRG = R.X * Grad.X + R.Y * Grad.Y + R.Z * Grad.Z;
                    P.DensityRate += -Par.DensityDiffusion * Par.SmoothingLength *
                        Par.SoundSpeed * (Par.ParticleMass / Neighbor.Density) * 2.0 *
                        (Neighbor.Density - P.Density) * DotRG / (R2 + 0.01 * H2);
                }

                double FullConstant = Par.ParticleMass *
                    (P.Pressure / (P.Density * P.Density) +
                     Neighbor.Pressure / (Neighbor.Density * Neighbor.Density));
                if (DotVR < 0.0)
                {
                    double Mu = Par.SmoothingLength * DotVR / (R2 + 0.01 * H2);
                    double MeanDensity = 0.5 * (P.Density + Neighbor.Density);
                    FullConstant += Par.ParticleMass *
                        (-Par.ArtificialViscosity * Par.SoundSpeed * Mu / MeanDensity);
                }

                P.Acceleration.X -= FullConstant * Grad.X;
                P.Acceleration.Y -= FullConstant * Grad.Y;
                P.Acceleration.Z -= FullConstant * Grad.Z;

                if (!Neighbor.IsBoundary)
                {
                    double Weight = Par.XSPHFactor *
                        (Par.ParticleMass / Neighbor.Density) * CSKernel(Distance, Par);
                    P.XSPH.X -= DV.X * Weight;
                    P.XSPH.Y -= DV.Y * Weight;
                    P.XSPH.Z -= DV.Z * Weight;
                }
            }
        }
    }
}

//Project water particles onto the accessible side of every solid surface.
//The ramp is treated as a height field so the ramp/floor junction cannot
//repeatedly push a particle back and forth across RampStartX.
//The tiny clearance is numerical (not an additional physical wall).
void RemoveIncomingNormalVelocity(Particle& P, Vector Normal,
    double WallNormalVelocity, const parameters& Par)
{
    double Relative = P.Velocity.X * Normal.X +
        P.Velocity.Y * Normal.Y + P.Velocity.Z * Normal.Z - WallNormalVelocity;
    if (Relative < 0.0)
    {
        double Correction = (1.0 + Par.WallRestitution) * Relative;
        P.Velocity.X -= Correction * Normal.X;
        P.Velocity.Y -= Correction * Normal.Y;
        P.Velocity.Z -= Correction * Normal.Z;
    }
}

void ReflectPillar(Particle& P, const Pillar& Obstacle,
    const parameters& Par, double Clearance)
{
    //Cylinder points along Z and starts on the flat seabed.
    const double Top = Par.MinDimmensions.Z + Obstacle.Height;
    double DX = P.Position.X - Obstacle.X;
    double DY = P.Position.Y - Obstacle.Y;
    double Radius = sqrt(DX * DX + DY * DY);

    //Only points inside the solid need correcting. Clearance is applied
    //to the chosen escape surface, not to detection of an intersection.
    if (Radius >= Obstacle.Radius || P.Position.Z >= Top ||
        P.Position.Z < Par.MinDimmensions.Z - Clearance) { return; }

    double ThroughSide = Obstacle.Radius - Radius;
    double ThroughTop = Top - P.Position.Z;

    if (ThroughTop < ThroughSide)
    {
        //Closest escape is upwards through the top cap.
        P.Position.Z = Top + Clearance;
        RemoveIncomingNormalVelocity(P, { 0.0, 0.0, 1.0 }, 0.0, Par);
    }
    else
    {
        //Closest escape is radially out through the curved surface.
        //The deterministic normal at the exact axis avoids division by zero.
        double NX = Radius > 1e-12 ? DX / Radius : 1.0;
        double NY = Radius > 1e-12 ? DY / Radius : 0.0;
        P.Position.X = Obstacle.X + (Obstacle.Radius + Clearance) * NX;
        P.Position.Y = Obstacle.Y + (Obstacle.Radius + Clearance) * NY;
        RemoveIncomingNormalVelocity(P, { NX, NY, 0.0 }, 0.0, Par);
    }
}

void ReflectParticle(Particle& P, double Time, const parameters& Par)
{
    //The floor MUST be resolved BEFORE the cylinders. A fast-moving water
    //particle may cross both a pillar side and the floor during one time step.
    //In the earlier version the cylinder routine ignored particles below Z=0,
    //then floor correction lifted them INSIDE the solid cylinder.
    const double Epsilon = max(1e-8, Par.ParticleSpacing * 1e-4);
    double Left = GetPistonPosition(Time, Par);
    double LeftVelocity = GetPistonVelocity(Time, Par);

    if (!isfinite(P.Position.X) || !isfinite(P.Position.Y) ||
        !isfinite(P.Position.Z) || !isfinite(P.Velocity.X) ||
        !isfinite(P.Velocity.Y) || !isfinite(P.Velocity.Z))
        throw runtime_error("Non-finite particle before boundary correction, ID=" +
            to_string(P.ID));

    //Clip the horizontal tank footprint, removing only incoming velocity.
    if (P.Position.X < Left + Epsilon)
    {
        P.Position.X = Left + Epsilon;
        RemoveIncomingNormalVelocity(P, { 1.0, 0.0, 0.0 }, LeftVelocity, Par);
    }
    if (P.Position.X > Par.MaxDimmensions.X - Epsilon)
    {
        P.Position.X = Par.MaxDimmensions.X - Epsilon;
        RemoveIncomingNormalVelocity(P, { -1.0, 0.0, 0.0 }, 0.0, Par);
    }
    if (P.Position.Y < Par.MinDimmensions.Y + Epsilon)
    {
        P.Position.Y = Par.MinDimmensions.Y + Epsilon;
        RemoveIncomingNormalVelocity(P, { 0.0, 1.0, 0.0 }, 0.0, Par);
    }
    if (P.Position.Y > Par.MaxDimmensions.Y - Epsilon)
    {
        P.Position.Y = Par.MaxDimmensions.Y - Epsilon;
        RemoveIncomingNormalVelocity(P, { 0.0, -1.0, 0.0 }, 0.0, Par);
    }

    //FIRST: Lift points above the flat floor or the sloping ramp.
    //Vertical correction avoids crossing the corner at RampStartX.
    double Floor = GetRampHeight(P.Position.X, Par);
    if (P.Position.Z < Floor + Epsilon)
    {
        P.Position.Z = Floor + Epsilon;
        Vector Normal = (P.Position.X > Par.RampStartX) ?
            GetRampNormal(Par) : Vector{ 0.0, 0.0, 1.0 };
        RemoveIncomingNormalVelocity(P, Normal, 0.0, Par);
    }
    if (P.Position.Z > Par.MaxDimmensions.Z - Epsilon)
    {
        P.Position.Z = Par.MaxDimmensions.Z - Epsilon;
        RemoveIncomingNormalVelocity(P, { 0.0, 0.0, -1.0 }, 0.0, Par);
    }

    //SECOND: Expel water from the vertical cylinders. Since the particles are
    //already above the floor, a cylinder side or cap projection is valid.
    ReflectPillar(P, Par.UnderwaterPillar, Par, Epsilon);
    ReflectPillar(P, Par.AboveWaterPillar, Par, Epsilon);

    //THIRD: Reapply seabed constraints for numerical roundoff. Both pillars
    //sit on the flat seabed, so this cannot conflict with their projections.
    Floor = GetRampHeight(P.Position.X, Par);
    if (P.Position.Z < Floor + Epsilon)
    {
        P.Position.Z = Floor + Epsilon;
        RemoveIncomingNormalVelocity(P,
            (P.Position.X > Par.RampStartX) ? GetRampNormal(Par) :
            Vector{ 0.0, 0.0, 1.0 }, 0.0, Par);
    }

    if (P.Position.Z < Floor - 1e-9 ||
        IsInsidePillar(P.Position, Par.UnderwaterPillar, Par) ||
        IsInsidePillar(P.Position, Par.AboveWaterPillar, Par))
    {
        throw runtime_error(
            "Boundary error at particle " + to_string(P.ID) +
            " X=" + to_string(P.Position.X) +
            " Y=" + to_string(P.Position.Y) +
            " Z=" + to_string(P.Position.Z) +
            " floor=" + to_string(Floor) +
            ". Inspect time step / velocities / geometry.");
    }
}

void SaveParticles(ofstream& SFile, const vector<Particle>& Population,
    int Iteration, double Time)
{
    for (const Particle& P : Population)
    {
        SFile << Iteration << "," << Time << "," << P.ID << "," <<
            P.Position.X << "," << P.Position.Y << "," << P.Position.Z << "," <<
            P.Velocity.X << "," << P.Velocity.Y << "," << P.Velocity.Z << "," <<
            P.Density << "," << P.Pressure << "\n";
    }
}

int CountInitialParticles(double Spacing, const parameters& Par, int StopAfter)
{
    //Count valid initial water positions while leaving the solid cylinders empty
    int Count = 0;
    for (double Z = Par.MinDimmensions.Z + 0.5 * Spacing;
         Z < Par.WaterLevel - 0.25 * Spacing; Z += Spacing)
    {
        for (double Y = Par.MinDimmensions.Y + 0.5 * Spacing;
             Y < Par.MaxDimmensions.Y - 0.5 * Spacing; Y += Spacing)
        {
            for (double X = Par.PistonBaseX + 0.5 * Spacing;
                 X < Par.MaxDimmensions.X - 0.5 * Spacing; X += Spacing)
            {
                Vector Position = { X, Y, Z };
                if (Z < GetRampHeight(X, Par) + 0.5 * Spacing) { continue; }
                if (IsInsidePillar(Position, Par.UnderwaterPillar, Par) ||
                    IsInsidePillar(Position, Par.AboveWaterPillar, Par)) { continue; }
                Count++;
                if (Count >= StopAfter) { return Count; }
            }
        }
    }
    return Count;
}

double GetParticleSpacing(const parameters& Par)
{
    //Estimate initial water volume, subtracting the ramp and submerged cylinder volumes
    double Width = Par.MaxDimmensions.Y - Par.MinDimmensions.Y;
    double Height = Par.WaterLevel - Par.MinDimmensions.Z;
    double Gradient = tan(Par.RampAngle * Par.Pi / 180.0);
    double WetEnd = min(Par.MaxDimmensions.X, Par.RampStartX + Height / Gradient);
    double InclinedLength = max(0.0, WetEnd - Par.RampStartX);
    double Volume = Width * (Height * (WetEnd - Par.PistonBaseX) -
        0.5 * Gradient * InclinedLength * InclinedLength);

    auto SubmergedCylinderVolume = [&](const Pillar& Obstacle)
    {
        double SubmergedHeight = min(Height, Obstacle.Height);
        return Par.Pi * Obstacle.Radius * Obstacle.Radius * SubmergedHeight;
    };
    Volume -= SubmergedCylinderVolume(Par.UnderwaterPillar);
    Volume -= SubmergedCylinderVolume(Par.AboveWaterPillar);
    if (Volume <= 0.0) { throw runtime_error("No space for initial water"); }

    double Estimate = cbrt(Volume / Par.N);
    double Lower = 0.65 * Estimate;
    double Upper = 1.35 * Estimate;
    for (int i = 0; i < 25 && CountInitialParticles(Lower, Par, Par.N) < Par.N; i++)
    {
        Lower *= 0.8;
    }
    for (int i = 0; i < 25 && CountInitialParticles(Upper, Par, Par.N) >= Par.N; i++)
    {
        Upper *= 1.2;
    }
    if (CountInitialParticles(Lower, Par, Par.N) < Par.N)
    {
        throw runtime_error("Could not determine particle spacing");
    }
    for (int i = 0; i < 32; i++)
    {
        double Middle = 0.5 * (Lower + Upper);
        if (CountInitialParticles(Middle, Par, Par.N) >= Par.N) { Lower = Middle; }
        else { Upper = Middle; }
    }
    return Lower * 0.999999;
}

void ValidatePillar(const Pillar& Obstacle, const parameters& Par,
    const string& Name)
{
    if (!isfinite(Obstacle.X) || !isfinite(Obstacle.Y) ||
        !isfinite(Obstacle.Radius) || !isfinite(Obstacle.Height) ||
        Obstacle.Radius <= 0.0 || Obstacle.Height <= 0.0 ||
        Obstacle.Height >= Par.MaxDimmensions.Z - Par.MinDimmensions.Z ||
        Obstacle.X - Obstacle.Radius <= Par.PistonBaseX + Par.PistonAmplitude ||
        Obstacle.X + Obstacle.Radius >= Par.RampStartX ||
        Obstacle.Y - Obstacle.Radius <= Par.MinDimmensions.Y ||
        Obstacle.Y + Obstacle.Radius >= Par.MaxDimmensions.Y)
    {
        throw runtime_error("Invalid position or size for " + Name);
    }
}

void ValidateParameters(const parameters& Par)
{
    if (Par.N <= 0 || Par.TimeStep <= 0.0 ||
        Par.SoundSpeed <= 0.0 || Par.ReferenceDensity <= 0.0 ||
        Par.Gamma <= 0.0 || Par.SmoothingLengthFactor <= 0.0 ||
        Par.SaveEveryTime < 0.0 || Par.DiagnosticsEveryTime < 0.0 ||
        Par.MaxRunHours < 0.0 || Par.PistonAmplitude < 0.0 ||
        Par.WaveFrequency < 0.0 || Par.WaveStartTime < 0.0 || Par.WaveRampTime < 0.0 ||
        Par.MaxDimmensions.X <= Par.MinDimmensions.X ||
        Par.MaxDimmensions.Y <= Par.MinDimmensions.Y ||
        Par.MaxDimmensions.Z <= Par.MinDimmensions.Z ||
        Par.WaterLevel <= Par.MinDimmensions.Z ||
        Par.WaterLevel >= Par.MaxDimmensions.Z ||
        Par.RampAngle <= 0.0 || Par.RampAngle >= 90.0 ||
        Par.RampStartX <= Par.PistonBaseX + Par.PistonAmplitude ||
        Par.RampStartX >= Par.MaxDimmensions.X ||
        Par.PistonBaseX - Par.PistonAmplitude <= Par.MinDimmensions.X ||
        GetRampHeight(Par.MaxDimmensions.X, Par) >= Par.MaxDimmensions.Z ||
        Par.InitialJitter < 0.0 || Par.InitialJitter >= 0.25 ||
        Par.WallRestitution < 0.0 || Par.WallRestitution > 1.0)
    {
        throw runtime_error("Invalid simulation parameters or tank geometry");
    }

    ValidatePillar(Par.UnderwaterPillar, Par, "UnderwaterPillar");
    ValidatePillar(Par.AboveWaterPillar, Par, "AboveWaterPillar");
    if (Par.UnderwaterPillar.Height >= Par.WaterLevel - Par.MinDimmensions.Z)
    {
        throw runtime_error("UnderwaterPillar must finish below the initial water level");
    }
    if (Par.AboveWaterPillar.Height <= Par.WaterLevel - Par.MinDimmensions.Z)
    {
        throw runtime_error("AboveWaterPillar must extend above the initial water level");
    }
    double DX = Par.UnderwaterPillar.X - Par.AboveWaterPillar.X;
    double DY = Par.UnderwaterPillar.Y - Par.AboveWaterPillar.Y;
    double MinimumGap = Par.UnderwaterPillar.Radius + Par.AboveWaterPillar.Radius;
    if (DX * DX + DY * DY <= MinimumGap * MinimumGap)
    {
        throw runtime_error("The two pillars overlap");
    }
}



//Python-facing library. The numerical methods above are retained from the
//original wave tank. There is no executable main(), step count, or wall-time
//limit in this library. Python supplies each batch of steps.
#if defined(_WIN32)
#define WT_API extern "C" __declspec(dllexport)
#else
#define WT_API extern "C" __attribute__((visibility("default")))
#endif

//Keep this C layout synchronized with _NativeConfig in wave_tank.py.
struct WT_Config
{
    int N;
    double TimeStep;
    double MinX, MinY, MinZ;
    double MaxX, MaxY, MaxZ;
    double PistonBaseX, PistonAmplitude, WaveFrequency;
    double WaveStartTime, WaveRampTime;
    double SoundSpeed, SmoothingLengthFactor;
    double DensityDiffusion, ArtificialViscosity, XSPHFactor;
    unsigned int RandomSeed;
    int WriteFiles;
    double SaveEveryTime;
};
struct WT_MapRow
{
    double X,Y,Z,VX,VY,VZ;
    int Count;
};
struct WT_ParticleRow
{
    int ID;
    double X,Y,Z,VX,VY,VZ,Density,Pressure;
};
struct WT_Info
{
    uint64_t Iteration;
    double Time, Spacing, SmoothingLength, ParticleMass;
};
static thread_local string WT_ErrorMessage;

struct WaveTank
{
    parameters Par;
    vector<Particle> Population;
    vector<Particle> BoundaryParticles;
    vector<Particle> AllParticles;
    unique_ptr<SpatialGrid> Grid;
    ofstream SFile;
    uint64_t Iteration = 0;
    double Time = 0.0;
    double NextSave = 0.0;

    explicit WaveTank(const WT_Config& C)
    {
        Par.N = C.N;
        Par.TimeStep = C.TimeStep;
        Par.MinDimmensions = {C.MinX,C.MinY,C.MinZ};
        Par.MaxDimmensions = {C.MaxX,C.MaxY,C.MaxZ};
        Par.SoundSpeed = C.SoundSpeed;
        Par.SmoothingLengthFactor = C.SmoothingLengthFactor;
        Par.DensityDiffusion = C.DensityDiffusion;
        Par.ArtificialViscosity = C.ArtificialViscosity;
        Par.XSPHFactor = C.XSPHFactor;
        Par.RandomSeed = C.RandomSeed;
        Par.WriteFiles = C.WriteFiles != 0;
        Par.SaveEveryTime = C.SaveEveryTime;
        ConfigureGeometry(Par);
        //The piston coordinates are supplied in METRES; they are not reset by
        //ConfigureGeometry. Python defaults to proportional scaled values.
        Par.PistonBaseX = C.PistonBaseX;
        Par.PistonAmplitude = C.PistonAmplitude;
        Par.WaveFrequency = C.WaveFrequency;
        Par.WaveStartTime = C.WaveStartTime;
        Par.WaveRampTime = C.WaveRampTime;
        ValidateParameters(Par);
        Par.ParticleSpacing = GetParticleSpacing(Par);
        Par.SmoothingLength = Par.SmoothingLengthFactor * Par.ParticleSpacing;
        Par.ParticleMass = Par.ReferenceDensity * pow(Par.ParticleSpacing, 3.0);
        double MaximumTimeStep = 0.22 * Par.SmoothingLength / Par.SoundSpeed;
        if (Par.TimeStep > MaximumTimeStep)
            throw runtime_error("TimeStep exceeds CFL safety bound " + to_string(MaximumTimeStep));

        vector<Vector> InitialPositions;
        InitialPositions.reserve((size_t)Par.N + Par.N/10);
        const double dx = Par.ParticleSpacing;
        for (double Z = Par.MinDimmensions.Z + 0.5*dx;
                    Z < Par.WaterLevel - 0.25*dx; Z += dx)
        for (double Y = Par.MinDimmensions.Y + 0.5*dx;
                    Y < Par.MaxDimmensions.Y - 0.5*dx; Y += dx)
        for (double X = Par.PistonBaseX + 0.5*dx;
                    X < Par.MaxDimmensions.X - 0.5*dx; X += dx)
        {
            Vector Position{X,Y,Z};
            if (Z < GetRampHeight(X,Par) + 0.5*dx ||
                IsInsidePillar(Position,Par.UnderwaterPillar,Par) ||
                IsInsidePillar(Position,Par.AboveWaterPillar,Par)) continue;
            InitialPositions.push_back(Position);
        }
        if (InitialPositions.size() < (size_t)Par.N)
            throw runtime_error("Not enough valid initial fluid positions");
        mt19937 Generator(Par.RandomSeed);
        shuffle(InitialPositions.begin(),InitialPositions.end(),Generator);
        uniform_real_distribution<double> Jitter(-Par.InitialJitter*dx,Par.InitialJitter*dx);
        Population.reserve(Par.N);
        for (int i=0;i<Par.N;i++)
        {
            Vector P=InitialPositions[(size_t)i];
            Vector Candidate{P.X+Jitter(Generator),P.Y+Jitter(Generator),P.Z+Jitter(Generator)};
            if (Candidate.X > Par.PistonBaseX + 0.25*dx &&
                Candidate.X < Par.MaxDimmensions.X - 0.25*dx &&
                Candidate.Y > Par.MinDimmensions.Y + 0.25*dx &&
                Candidate.Y < Par.MaxDimmensions.Y - 0.25*dx &&
                Candidate.Z > GetRampHeight(Candidate.X,Par) + 0.25*dx &&
                Candidate.Z < Par.WaterLevel - 0.05*dx &&
                !IsInsidePillar(Candidate,Par.UnderwaterPillar,Par) &&
                !IsInsidePillar(Candidate,Par.AboveWaterPillar,Par)) P=Candidate;
            double Pressure=Par.ReferenceDensity*(-Par.GravAccelaration)*(Par.WaterLevel-P.Z);
            Population.emplace_back(i,P,Vector{},GetDensityFromPressure(Pressure,Par));
            Population.back().Pressure=Pressure;
        }
        Grid=make_unique<SpatialGrid>(Par);
        BoundaryParticles.reserve(Population.size()*6);
        AllParticles.reserve(Population.size()*7);
        NextSave=Par.SaveEveryTime;
        if (Par.WriteFiles)
        {
            SFile.open("LogBook.txt");
            if (!SFile) throw runtime_error("Cannot open LogBook.txt");
            SFile << "Iteration,Time,ID,XPos,YPos,ZPos,XVel,YVel,ZVel,Density,Pressure\n" << setprecision(15);
            ofstream F("SimulationConfig.txt");
            if (!F) throw runtime_error("Cannot open SimulationConfig.txt");
            F << setprecision(15) << "Parameter,Value\n";
            auto Out=[&](const string& K,double V){ F << K << "," << V << "\n"; };
            Out("TankLength",Par.MaxDimmensions.X-Par.MinDimmensions.X);
            Out("TankWidth",Par.MaxDimmensions.Y-Par.MinDimmensions.Y);
            Out("TankHeight",Par.MaxDimmensions.Z-Par.MinDimmensions.Z);
            Out("MinX",Par.MinDimmensions.X);
            Out("MinY",Par.MinDimmensions.Y);
            Out("MinZ",Par.MinDimmensions.Z);
            Out("WaterLevel",Par.WaterLevel);
            Out("RampStartX",Par.RampStartX);
            Out("RampAngle",Par.RampAngle);
            Out("PistonBaseX",Par.PistonBaseX);
            Out("PistonAmplitude",Par.PistonAmplitude);
            Out("WaveFrequency",Par.WaveFrequency);
            Out("WaveStartTime",Par.WaveStartTime);
            Out("WaveRampTime",Par.WaveRampTime);
            Out("Pillar1X",Par.UnderwaterPillar.X);
            Out("Pillar1Y",Par.UnderwaterPillar.Y);
            Out("Pillar1Radius",Par.UnderwaterPillar.Radius);
            Out("Pillar1Height",Par.UnderwaterPillar.Height);
            Out("Pillar2X",Par.AboveWaterPillar.X);
            Out("Pillar2Y",Par.AboveWaterPillar.Y);
            Out("Pillar2Radius",Par.AboveWaterPillar.Radius);
            Out("Pillar2Height",Par.AboveWaterPillar.Height);
            Out("ParticleSpacing",Par.ParticleSpacing);
            Out("ParticleCount",Par.N);
            Out("TimeStep",Par.TimeStep);
            F.close();
            SaveParticles(SFile,Population,0,0.0);
        }
    }

    //No fixed end iteration. uint64_t counters are practically unbounded for
    //real simulation durations; Python may call Step repeatedly.
    void Step(uint64_t Steps)
    {
        for (uint64_t k=0;k<Steps;k++)
        {
            CreateBoundaryParticles(Population,BoundaryParticles,Time,Par);
            AllParticles.clear();
            AllParticles.insert(AllParticles.end(),Population.begin(),Population.end());
            AllParticles.insert(AllParticles.end(),BoundaryParticles.begin(),BoundaryParticles.end());
            Grid->Build(AllParticles);
            CalculateTotalForce(Population,BoundaryParticles,*Grid,Par);
            const double NextTime=Time+Par.TimeStep;
            for (auto& P:Population)
            {
                P.Velocity.X+=Par.TimeStep*P.Acceleration.X;
                P.Velocity.Y+=Par.TimeStep*P.Acceleration.Y;
                P.Velocity.Z+=Par.TimeStep*P.Acceleration.Z;
                P.Density+=Par.TimeStep*P.DensityRate;
                P.Density=max(0.95*Par.ReferenceDensity,min(1.05*Par.ReferenceDensity,P.Density));
                P.Position.X+=Par.TimeStep*(P.Velocity.X+P.XSPH.X);
                P.Position.Y+=Par.TimeStep*(P.Velocity.Y+P.XSPH.Y);
                P.Position.Z+=Par.TimeStep*(P.Velocity.Z+P.XSPH.Z);
                ReflectParticle(P,NextTime,Par);
                P.Pressure=GetPressure(P.Density,Par);
                if (!isfinite(P.Density) || !isfinite(P.Pressure) ||
                    !isfinite(P.Velocity.X) || !isfinite(P.Velocity.Y) ||
                    !isfinite(P.Velocity.Z))
                    throw runtime_error("Invalid particle state at step " + to_string(Iteration));
            }
            Time=NextTime;
            Iteration++;
            if (Par.WriteFiles && Par.SaveEveryTime>0.0 &&
                Time >= NextSave-0.5*Par.TimeStep)
            {
                //Iteration stored as integer string; SaveParticles originally
                //takes int, so use an explicit wide-ID serialization here.
                for (const Particle& P:Population)
                    SFile << Iteration << "," << Time << "," << P.ID << "," <<
                        P.Position.X << "," << P.Position.Y << "," << P.Position.Z << "," <<
                        P.Velocity.X << "," << P.Velocity.Y << "," << P.Velocity.Z << "," <<
                        P.Density << "," << P.Pressure << "\n";
                NextSave+=Par.SaveEveryTime;
            }
        }
    }

    //Refresh the spatial hash at the CURRENT post-step positions, without
    //virtual boundary particles: the returned map contains only real water.
    vector<MapVector> VelocityMap()
    {
        Grid->Build(Population);
        return outputMap(Population,*Grid);
    }
};

WT_API const char* WT_LastError() { return WT_ErrorMessage.c_str(); }
WT_API void* WT_Create(const WT_Config* C)
{
    try
    {
        if (!C) throw runtime_error("Null config");
        auto Result=make_unique<WaveTank>(*C);
        WT_ErrorMessage.clear();
        return Result.release();
    }
    catch(const exception& E) { WT_ErrorMessage=E.what(); return nullptr; }
}
WT_API void WT_Destroy(void* H)
{
    delete static_cast<WaveTank*>(H);
}
WT_API int WT_Step(void* H,uint64_t Steps)
{
    try
    {
        if (!H) throw runtime_error("Null simulation handle");
        static_cast<WaveTank*>(H)->Step(Steps);
        WT_ErrorMessage.clear();
        return 0;
    }
    catch(const exception& E) { WT_ErrorMessage=E.what(); return -1; }
}
WT_API int64_t WT_GetMap(void* H,WT_MapRow* Rows,uint64_t Capacity)
{
    try
    {
        if (!H) throw runtime_error("Null simulation handle");
        auto Map=static_cast<WaveTank*>(H)->VelocityMap();
        if (Rows && Capacity < Map.size()) throw runtime_error("Map buffer too small");
        if (Rows)
            for (size_t i=0;i<Map.size();i++)
            {
                const auto& M=Map[i];
                Rows[i]={M.Position.X,M.Position.Y,M.Position.Z,
                         M.Velocity.X,M.Velocity.Y,M.Velocity.Z,M.ParticleCount};
            }
        WT_ErrorMessage.clear();
        return (int64_t)Map.size();
    }
    catch(const exception& E) { WT_ErrorMessage=E.what(); return -1; }
}
WT_API int64_t WT_GetParticles(void* H,WT_ParticleRow* Rows,uint64_t Capacity)
{
    try
    {
        if (!H) throw runtime_error("Null simulation handle");
        const auto& Pop=static_cast<WaveTank*>(H)->Population;
        if (Rows && Capacity<Pop.size()) throw runtime_error("Particle buffer too small");
        if (Rows)
            for(size_t i=0;i<Pop.size();i++)
            {
                const auto& P=Pop[i];
                Rows[i]={P.ID,P.Position.X,P.Position.Y,P.Position.Z,
                         P.Velocity.X,P.Velocity.Y,P.Velocity.Z,P.Density,P.Pressure};
            }
        WT_ErrorMessage.clear();
        return (int64_t)Pop.size();
    }
    catch(const exception& E) { WT_ErrorMessage=E.what(); return -1; }
}
WT_API int WT_GetInfo(void* H,WT_Info* Info)
{
    if (!H || !Info) { WT_ErrorMessage="Invalid handle/info"; return -1; }
    const WaveTank& S=*static_cast<WaveTank*>(H);
    *Info={S.Iteration,S.Time,S.Par.ParticleSpacing,S.Par.SmoothingLength,S.Par.ParticleMass};
    WT_ErrorMessage.clear();
    return 0;
}

//For arbitrary Python-driven piston trajectories, set a position and velocity
//before each sim.step(1). The state is held constant for the batch.
WT_API int WT_SetPistonState(void* H,double Position,double Velocity)
{
    try
    {
        if (!H) throw runtime_error("Null simulation handle");
        WaveTank& S=*static_cast<WaveTank*>(H);
        const parameters& P=S.Par;
        if (!isfinite(Position) || !isfinite(Velocity) ||
            Position <= P.MinDimmensions.X ||
            Position >= P.RampStartX ||
            Position >= P.UnderwaterPillar.X-P.UnderwaterPillar.Radius)
            throw runtime_error("Custom piston position outside safe flume range");
        S.Par.UseCustomPiston=true;
        S.Par.CustomPistonPosition=Position;
        S.Par.CustomPistonVelocity=Velocity;
        WT_ErrorMessage.clear();
        return 0;
    }
    catch(const exception& E) { WT_ErrorMessage=E.what(); return -1; }
}
WT_API int WT_UseSinusoidalPiston(void* H)
{
    if (!H) { WT_ErrorMessage="Null simulation handle"; return -1; }
    static_cast<WaveTank*>(H)->Par.UseCustomPiston=false;
    WT_ErrorMessage.clear();
    return 0;
}
