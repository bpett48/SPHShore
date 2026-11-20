"""SPHShore — a C++-accelerated 3D SPH wave-tank simulator.

Example
-------
>>> from sphshore import Parameters, WaveTank
>>> with WaveTank(Parameters(N=500, TimeStep=0.00005)) as tank:
...     tank.step(10)
...     print(tank.positions().shape)
(500, 3)
"""
from .simulation import Parameters, WaveTank

__version__ = "0.1.0"
__all__ = ["Parameters", "WaveTank", "__version__"]
