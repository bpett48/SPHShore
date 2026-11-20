"""Python interface to the original C++17 SPHShore solver."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Optional, Tuple
import math
import numbers

import numpy as np
from . import _core


@dataclass
class Parameters:
    """Set up one fluid simulation; lengths in metres and time in seconds.

    The default geometry includes a moving piston, a beach ramp, and two
    vertical cylinders. Resizing MaxDimmensions rescales the whole geometry.
    All dimensional defaults are taken from the existing wave-tank solver.
    """

    MaxDimmensions: Tuple[float, float, float] = (1.20, 0.45, 0.30)
    MinDimmensions: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    N: int = 12000
    TimeStep: float = 0.00010

    PistonBaseX: Optional[float] = None
    PistonAmplitude: Optional[float] = None
    WaveFrequency: float = 1.6
    WaveStartTime: float = 0.25
    WaveRampTime: float = 0.9

    SoundSpeed: float = 25.0
    SmoothingLengthFactor: float = 1.40
    DensityDiffusion: float = 0.10
    ArtificialViscosity: float = 0.055
    XSPHFactor: float = 0.05
    RandomSeed: int = 42069
    WriteFiles: bool = False
    SaveEveryTime: float = 0.025

    MapEveryIterations: int = 100

    def _native_config(self) -> dict:
        """Validate arguments and create the configuration expected by C++."""
        def coordinates(name, values):
            if len(values) != 3:
                raise ValueError(f"{name} must have exactly three coordinates")
            out = tuple(float(x) for x in values)
            if not all(math.isfinite(x) for x in out):
                raise ValueError(f"{name} must contain finite values")
            return out

        minimum = coordinates("MinDimmensions", self.MinDimmensions)
        maximum = coordinates("MaxDimmensions", self.MaxDimmensions)
        if any(b <= a for a, b in zip(minimum, maximum)):
            raise ValueError("Every tank maximum must exceed its minimum")
        if isinstance(self.N, bool) or not isinstance(self.N, numbers.Integral) or not 1 <= self.N <= 2**31 - 1:
            raise ValueError("N must be an integer between 1 and 2^31-1")
        if isinstance(self.MapEveryIterations, bool) or not isinstance(self.MapEveryIterations, numbers.Integral) or self.MapEveryIterations <= 0:
            raise ValueError("MapEveryIterations must be positive")
        if isinstance(self.RandomSeed, bool) or not isinstance(self.RandomSeed, numbers.Integral) or not 0 <= self.RandomSeed <= 2**32 - 1:
            raise ValueError("RandomSeed must be an unsigned 32-bit integer")
        scale_x = (maximum[0] - minimum[0]) / 1.20
        piston_base = minimum[0] + 0.03 * scale_x if self.PistonBaseX is None else float(self.PistonBaseX)
        amplitude = 0.007 * scale_x if self.PistonAmplitude is None else float(self.PistonAmplitude)
        doubles = dict(
            TimeStep=self.TimeStep, PistonBaseX=piston_base, PistonAmplitude=amplitude,
            WaveFrequency=self.WaveFrequency, WaveStartTime=self.WaveStartTime,
            WaveRampTime=self.WaveRampTime, SoundSpeed=self.SoundSpeed,
            SmoothingLengthFactor=self.SmoothingLengthFactor,
            DensityDiffusion=self.DensityDiffusion,
            ArtificialViscosity=self.ArtificialViscosity, XSPHFactor=self.XSPHFactor,
            SaveEveryTime=self.SaveEveryTime,
        )
        doubles = {k: float(v) for k, v in doubles.items()}
        if not all(math.isfinite(v) for v in doubles.values()):
            raise ValueError("Simulation parameters must be finite")
        if doubles['TimeStep'] <= 0:
            raise ValueError("TimeStep must be positive")
        return dict(
            N=int(self.N), RandomSeed=int(self.RandomSeed), WriteFiles=int(bool(self.WriteFiles)),
            MinX=minimum[0], MinY=minimum[1], MinZ=minimum[2],
            MaxX=maximum[0], MaxY=maximum[1], MaxZ=maximum[2], **doubles,
        )


class WaveTank:
    """One persistent simulation powered by C++17.

    The simulation is not given a final iteration number. Call step() or
    run_maps() as often as needed. Use a context manager to release resources.
    """

    def __init__(self, parameters: Optional[Parameters] = None, **kwargs):
        if parameters is not None and kwargs:
            raise TypeError("Pass a Parameters instance OR keyword settings, not both")
        self.parameters = parameters if parameters is not None else Parameters(**kwargs)
        if not isinstance(self.parameters, Parameters):
            raise TypeError("parameters must be a Parameters instance")
        self._handle = _core.create(self.parameters._native_config())

    def _check_open(self):
        if self._handle is None:
            raise RuntimeError("This WaveTank has been closed")

    def close(self) -> None:
        if self._handle is not None:
            _core.close(self._handle)
            self._handle = None

    def __enter__(self):
        self._check_open()
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()

    def __del__(self):
        if getattr(self, "_handle", None) is not None:
            self.close()

    def step(self, iterations: int = 1) -> None:
        """Advance any nonnegative count of steps (bounded by uint64 per call)."""
        self._check_open()
        if isinstance(iterations, bool) or not isinstance(iterations, numbers.Integral) or not 0 <= iterations <= 2**64-1:
            raise ValueError("iterations must be a nonnegative uint64 integer")
        _core.step(self._handle, int(iterations))

    def info(self) -> dict:
        self._check_open()
        return _core.info(self._handle)

    @property
    def time(self) -> float:
        return self.info()["time"]

    @property
    def iteration(self) -> int:
        return self.info()["iteration"]

    def set_piston_state(self, position: float, velocity: float) -> None:
        """Override the piston with a custom X position and X velocity.

        For time-varying trajectories, set the state before each step(1).
        """
        self._check_open()
        _core.set_piston(self._handle, float(position), float(velocity))

    def use_sinusoidal_piston(self) -> None:
        """Return to the configured sinusoidal wave-maker motion."""
        self._check_open()
        _core.use_sinusoidal_piston(self._handle)

    def outputMap(self) -> np.ndarray:
        """Return an (M,7) NumPy array of occupied hash cells.

        Columns: x, y, z, vx, vy, vz, count; spatial units metres and m/s.
        Only real fluid particles contribute to the map.
        """
        self._check_open()
        return np.asarray(_core.output_map(self._handle), dtype=np.float64).reshape(-1, 7)

    def velocity_map(self) -> np.ndarray:
        """Pythonic alias for outputMap()."""
        return self.outputMap()

    def particles(self) -> np.ndarray:
        """Return an (N,9) NumPy array of current particle states.

        Columns: id, x, y, z, vx, vy, vz, density, pressure.
        """
        self._check_open()
        return np.asarray(_core.particles(self._handle), dtype=np.float64).reshape(-1, 9)

    def positions(self) -> np.ndarray:
        """Return an independent (N,3) position array."""
        return self.particles()[:, 1:4].copy()

    def velocities(self) -> np.ndarray:
        """Return an independent (N,3) velocity array."""
        return self.particles()[:, 4:7].copy()

    def run_maps(self, iterations: int, every: Optional[int] = None) -> Iterator[tuple[int, float, np.ndarray]]:
        """Advance and yield (iteration, simulated_time, velocity_map) batches.

        Includes the last partial batch. The solver has no overall run limit.
        """
        interval = self.parameters.MapEveryIterations if every is None else every
        if isinstance(interval, bool) or not isinstance(interval, numbers.Integral) or interval <= 0:
            raise ValueError("every must be a positive integer")
        if isinstance(iterations, bool) or not isinstance(iterations, numbers.Integral) or iterations < 0:
            raise ValueError("iterations must be a nonnegative integer")
        left = int(iterations)
        while left > 0:
            batch = min(left, int(interval))
            self.step(batch)
            left -= batch
            data = self.info()
            yield (data["iteration"], data["time"], self.outputMap())
