import math
from pathlib import Path

import numpy as np
import pytest

from sphshore import Parameters, WaveTank


def small(**overrides):
    defaults = dict(N=180, TimeStep=0.00005, MapEveryIterations=3, WriteFiles=False)
    defaults.update(overrides)
    return Parameters(**defaults)


def test_particles_and_map():
    with WaveTank(small()) as sim:
        assert sim.iteration == 0
        assert sim.time == 0
        first = sim.positions()
        assert first.shape == (180, 3)
        sim.step(3)
        assert sim.iteration == 3
        assert sim.time == pytest.approx(0.00015)
        assert sim.particles().shape == (180, 9)
        assert sim.velocities().shape == (180, 3)
        data = sim.outputMap()
        assert data.ndim == 2 and data.shape[1] == 7
        assert int(np.sum(data[:, 6])) == 180
        assert np.all(np.isfinite(data))
        assert np.isfinite(sim.particles()).all()
        assert np.linalg.norm(sim.positions() - first) > 0


def test_map_interval_and_no_file_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with WaveTank(small()) as sim:
        batches = list(sim.run_maps(7))
        assert [b[0] for b in batches] == [3, 6, 7]
        assert [b[2].shape[1] for b in batches] == [7, 7, 7]
        assert sim.iteration == 7
    assert list(tmp_path.iterdir()) == []


def test_file_output_opt_in(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with WaveTank(small(WriteFiles=True, SaveEveryTime=0.0001)) as sim:
        sim.step(3)
    assert (tmp_path / "LogBook.txt").is_file()
    assert (tmp_path / "SimulationConfig.txt").is_file()
    assert "Pillar2Height" in (tmp_path / "SimulationConfig.txt").read_text()


def test_scaled_geometry_with_config(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with WaveTank(small(MaxDimmensions=(0.6, 0.225, 0.15), WriteFiles=True)):
        pass
    lines = (tmp_path / "SimulationConfig.txt").read_text().splitlines()[1:]
    config = {line.split(",")[0]: float(line.split(",")[1]) for line in lines}
    assert config["WaterLevel"] == pytest.approx(0.065)
    assert config["Pillar1Height"] == pytest.approx(0.0375)
    assert config["Pillar2Height"] == pytest.approx(0.095)
    assert config["RampStartX"] == pytest.approx(0.32)
    assert config["PistonAmplitude"] == pytest.approx(0.0035)


def test_piston_override_and_restore():
    with WaveTank(small()) as sim:
        sim.set_piston_state(0.032, 0.02)
        sim.step(1)
        sim.use_sinusoidal_piston()
        sim.step(1)
        assert sim.iteration == 2


def test_bad_config_and_steps():
    with pytest.raises(ValueError):
        WaveTank(small(MaxDimmensions=(0.0, 0.45, 0.30)))
    with pytest.raises(ValueError):
        WaveTank(small(MapEveryIterations=0))
    with WaveTank(small()) as sim:
        with pytest.raises(ValueError):
            sim.step(-1)
        with pytest.raises(ValueError):
            sim.step(1.5)
        with pytest.raises(RuntimeError):
            sim.set_piston_state(-0.1, 0)
        with pytest.raises(ValueError):
            list(sim.run_maps(1, every=0))
        sim.close()
        with pytest.raises(RuntimeError):
            sim.step(1)
