"""Run a short SPHShore example and inspect the resulting NumPy arrays."""
from sphshore import Parameters, WaveTank

settings = Parameters(
    MaxDimmensions=(1.20, 0.45, 0.30),
    N=1200,
    TimeStep=0.00005,
    MapEveryIterations=20,
    WaveFrequency=1.6,
    PistonAmplitude=0.007,
    WriteFiles=False,
)

with WaveTank(settings) as simulation:
    for iteration, time, velocity_map in simulation.run_maps(100):
        print(f"step={iteration:4d} time={time:.5f}s cells={len(velocity_map)}")
    print("Positions:", simulation.positions().shape)
    print("Velocities:", simulation.velocities().shape)
    print("Map columns: x y z vx vy vz particle_count")
