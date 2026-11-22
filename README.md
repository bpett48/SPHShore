# SPHShore

SPHShore is a Python library for simulating water movement in a three-dimensional tank. It uses C++ to perform the calculations and Python to control the simulation.

The simulation includes a moving piston that creates waves, a sloping shoreline and two cylindrical pillars. All of which interact with the water.

This is for use of a third year project at university, things might still by buggy.

## Installation

Install SPHShore using pip:

```bash
pip install SPHShore
```

The library requires Python and NumPy. The simulation calculations are performed by a compiled C++ solver.

## What can SPHShore do?

This library can:

- **Create waves:** Use a moving piston to generate waves in the tank.
- **Change the tank size:** Set the tank dimensions in Python. The shoreline, pillars and piston are scaled to match.
- **Choose the number of particles:** Control how many particles are used to represent the water.
- **Control simulation speed:** Set the timestep and number of iterations.
- **Apply forces:** Add a force to a selected region of water using the spatial hashing grid.
- **View water movement:** Retrieve the average velocity of the particles in each occupied grid cell.
- **Access individual particles:** Retrieve the position and velocity of every water particle.

## Creating a simulation

Starting bits:

```python
import sphshore

Par = sphshore.Parameters(
    N=1200,
    TimeStep=0.00005,
    MaxDimmensions=(1.20, 0.45, 0.30)
)

with sphshore.WaveTank(Par) as Tank:
    Tank.step(100)

    Positions = Tank.positions()
    Velocities = Tank.velocities()

    print(Positions)
    print(Velocities)
```

`N` sets the number of particles, `TimeStep` controls the time between calculations, and `MaxDimmensions` sets the size of the tank in metres. There is no max number of iterations. 

## Adding propulsion

You can add force to a specific point. This force is distributed between the particles at this area.

```python
import sphshore

Par = sphshore.Parameters(
    N=1200,
    TimeStep=0.00005
)

with sphshore.WaveTank(Par) as Tank:

    Map = Tank.outputMap()

    #Choose an occupied cell
    Position = Map[0, :3]

    #Apply a total force of 0.01 N in the X direction
    Count = Tank.add_propulsion(
        Position,
        (0.01, 0.0, 0.0)
    )

    print("Particles affected:", Count)

    Tank.step(1)
```
The force is applied to the next timestep. To apply continuous propulsion, call "add_propulsion()" before each timestep.

## Getting simulation information

You can get information about the simulation in two ways:

### 1. Velocity map

```python
Map = Tank.outputMap()
```
This gives you the location and average particle velocity for each cell

Each row contains:

| X | Y | Z | VX | VY | VZ | Particle count |
|---|---|---|---|---|---|---|

This is the quicker way.

### 2. Individual particles

```python
Positions = Tank.positions()
Velocities = Tank.velocities()
```

`positions()` returns the X, Y and Z coordinates of every particle.

`velocities()` returns the X, Y and Z velocity components of every particle.

You can also use:

```python
Particles = Tank.particles()
```

This gives a list of all the particles and their properties, might be quicker idk.

## Controlling the piston

The piston is what creates the waves.
You can change its frequency, amplitude and starting time through the simulation parameters. 
You can set attributes directly or have it run on a sin wave.

```python
Tank.set_piston_state(0.032, 0.02)
```

To return to the automatic wavemaker:

```python
Tank.use_sinusoidal_piston()
```

## How the simulation works

SPHShore uses a smooth particle hydrodynamic model (SPH), specifically using a 3d spline kernel.
The fluid is represented using a large number of particles, each of which respond to the environment and each other.
The simulation uses a spatial hashing function to make things faster and transfer information.
The simulation uses spatial hashing to find nearby particles more efficiently.
All of the heavy lifting is done in the c++.

Email me in event of anything bpett48@gmail.com


