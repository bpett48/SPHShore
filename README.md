# SPHShore

**SPHShore** is an early-stage, C++17-accelerated Python package for simulating water waves in a three-dimensional smoothed particle hydrodynamics (SPH) tank. It retains the original wave-tank solver, with a moving piston wavemaker, an inclined shoreline and two fixed vertical pillars (one submerged, one extending above the still-water level).
## Install locally

Python 3.10+ and a C++17 compiler are required when building from source. On Windows, install Visual Studio Build Tools with the C++ workload; on Linux/macOS use a suitable C++ compiler.

```bash
python -m pip install .
```

For development:

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
```

Once a release is actually uploaded to PyPI and wheels are built for your platform, the installation command will be:

```bash
python -m pip install SPHShore
```

**Do not assume the PyPI command works until the owner has published a release.**

## Quick start

```python
import sphshore

config = sphshore.Parameters(
    MaxDimmensions=(1.20, 0.45, 0.30),
    N=1200,
    TimeStep=0.00005,
    MapEveryIterations=50,
    WriteFiles=False,
)

with sphshore.WaveTank(config) as tank:
    for iteration, time, velocity_map in tank.run_maps(200):
        print(iteration, time, velocity_map.shape)
    positions = tank.positions()   # shape (N, 3), metres
    velocities = tank.velocities() # shape (N, 3), metres/s
    particles = tank.particles()   # shape (N, 9)
```

### Velocity map

`outputMap()` or `velocity_map()` uses the C++ spatial hashing grid and returns an `(M, 7)` `float64` NumPy array with columns:

| 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| x | y | z | vx | vy | vz | particle count |

Cell positions are in metres and average velocities in metres per second. The map includes only **occupied** cells. It is produced on demand and held in memory; no files are needed.

### Fluid particles

`particles()` returns an `(N, 9)` `float64` array:

`id, x, y, z, vx, vy, vz, density, pressure`.

`positions()` and `velocities()` return copies with three columns each.

### Adjustable parameters

| Python parameter | Meaning |
|---|---|
| `MaxDimmensions=(L,W,H)` | Dimensions in metres; rescales the shoreline, water depth and both pillars |
| `MinDimmensions=(x,y,z)` | Tank origin |
| `N` | Number of water particles |
| `TimeStep` | Solver integration timestep in seconds, checked against a CFL bound |
| `MapEveryIterations` | Python-side interval between successive returned velocity maps |
| `PistonBaseX` | Piston neutral X-position (metres); `None` means a scale-dependent default |
| `PistonAmplitude` | Piston oscillation amplitude (metres); `None` means a scale-dependent default |
| `WaveFrequency` | Sinusoidal piston frequency, in Hz |
| `WaveStartTime` | Piston start, in simulated seconds |
| `WaveRampTime` | Smooth piston startup duration |
| `WriteFiles` | Opt-in CSV log and configuration output; default `False` |
| `SaveEveryTime` | Simulated seconds between log snapshots when writing files |

Changing all three dimensions by the same factor preserves geometric similarity. For unequal X, Y and Z scaling, horizontal pillar cross sections stay circular; the ramp slope adjusts accordingly. Changing the tank scale does *not* automatically tune the timestep, wave period or numerical stability parameters.

### Arbitrary run length

The C++ library has no fixed total-iteration cutoff. Python controls stepping:

```python
with sphshore.WaveTank(config) as tank:
    for _ in range(100):
        tank.step(100)
```

Calls to `step()` accept a nonnegative 64-bit integer count. The cumulative counter is also finite (`uint64_t`); practical limits include time, memory and numerical stability.

### Custom piston motion

```python
with sphshore.WaveTank(config) as tank:
    for _ in range(100):
        t = tank.time
        tank.set_piston_state(position=0.03, velocity=0.0)
        tank.step(1)
    tank.use_sinusoidal_piston()
```

Change the custom state before every individual timestep for non-constant motion. The position must stay within the safe flume region.

### No disk output

`WriteFiles=False` (default) means no `LogBook.txt` or `SimulationConfig.txt` is created by the solver. All state and maps can be requested as NumPy arrays. When `WriteFiles=True`, the existing CSV schema remains available for the older animation scripts.

### AddPropulsion

The original C++ `AddPropulsion` grid helper is retained in the solver source, but **not invoked or exposed** in this release, respecting the earlier request to keep it unused.

## Building packages for distribution

```bash
python -m pip install build
python -m build
```

The output is placed in `dist/`. A platform-specific wheel requires a matching OS/Python/architecture; GitHub Actions uses `cibuildwheel` to produce supported wheels and run tests.

### Release checklist

1. **Review the MIT licence and owner/author attribution**, and verify rights to all included code.
2. Put this directory into a GitHub repository you control.
3. Enable CI and review/test cross-platform wheels.
4. Create a PyPI account, verify ownership of the desired **SPHShore** project name, and configure a **Trusted Publisher** for the repository and `.github/workflows/release.yml` (environment `pypi`).
5. Update the version in `pyproject.toml` and `src/sphshore/__init__.py` before each release.
6. Create a GitHub Release; the release workflow builds wheels and source distribution, then requests publishing through PyPI Trusted Publishing.
7. Verify `python -m pip install SPHShore` in a fresh environment.

No GitHub repository or PyPI project is created automatically by this source code. Publishing requires the owner's account authorization and review.

## Scientific and implementation notes

- Weakly compressible SPH with a 3D cubic-spline kernel, pressure/viscous terms, density diffusion and spatial hashing.
- The water is constrained by a sloping rigid bed, moving piston, side walls and two fixed cylindrical obstacles.
- SPH boundary conditions are approximate; the geometry/solver has not been fully validated against experiments.
- Python calls a compiled C++ extension built from the existing solver source; it is not a Python rewrite of the physics.
- Copyright notices on the original source are retained. This package's licence should be approved by all code rights-holders before publication.

## Contributing

Bug reports, benchmark comparisons and reproducible test cases are welcome after a public repository has been created. Include Python version, platform, compiler, timestep and tank geometry for any numerical issue.
