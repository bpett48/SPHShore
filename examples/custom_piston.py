"""Control a piston trajectory from Python rather than the built-in sinusoid."""
import math
from sphshore import Parameters, WaveTank

parameters = Parameters(N=500, TimeStep=0.00005, WriteFiles=False)
with WaveTank(parameters) as simulation:
    for _ in range(100):
        t = simulation.time
        # Small, smooth, 1 Hz movement about the original piston base.
        pos = 0.03 + 0.002 * math.sin(2 * math.pi * t)
        vel = 0.002 * 2 * math.pi * math.cos(2 * math.pi * t)
        simulation.set_piston_state(pos, vel)
        simulation.step(1)
    simulation.use_sinusoidal_piston()
    print(simulation.info())
