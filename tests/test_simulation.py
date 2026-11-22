
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


def test_no_file_output_with_writefiles_enabled(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with WaveTank(small(WriteFiles=True, SaveEveryTime=0.0001)) as sim:
        sim.step(3)
        assert sim.iteration == 3
        assert sim.positions().shape == (180, 3)
        assert np.isfinite(sim.velocities()).all()

#No logbook files should be created
    assert not (tmp_path / "LogBook.txt").exists()
    assert not (tmp_path / "SimulationConfig.txt").exists()
    assert list(tmp_path.iterdir()) == []


def test_scaled_geometry():
    with WaveTank(small()) as original:
        OriginalPositions = original.positions()
        OriginalInfo = original.info()

    with WaveTank(small(
        MaxDimmensions=(0.6, 0.225, 0.15),
        TimeStep=0.000025
    )) as scaled:
        ScaledPositions = scaled.positions()
        ScaledInfo = scaled.info()

        assert ScaledPositions.shape == (180, 3)

        #The entire geometry and initial particle positions should scale
        assert np.allclose(
            ScaledPositions,
            OriginalPositions * 0.5,
            rtol=1e-6,
            atol=1e-9
        )

        assert ScaledInfo["particle_spacing"] == pytest.approx(
            OriginalInfo["particle_spacing"] * 0.5
        )

        assert ScaledInfo["smoothing_length"] == pytest.approx(
            OriginalInfo["smoothing_length"] * 0.5
        )

        assert ScaledInfo["particle_mass"] == pytest.approx(
            OriginalInfo["particle_mass"] * 0.125
        )

        scaled.step(3)
        assert scaled.iteration == 3
        assert np.isfinite(scaled.particles()).all()


def test_piston_override_and_restore():
    with WaveTank(small()) as sim:
        sim.set_piston_state(0.032, 0.02)
        sim.step(1)
        sim.use_sinusoidal_piston()
        sim.step(1)
        assert sim.iteration == 2


def test_add_propulsion():
    #Two identical simulations are used to measure only the added force
    with WaveTank(small()) as Normal, WaveTank(small()) as Propelled:

        Map = Normal.outputMap()

        #Choose the cell with the most particles
        Cell = Map[np.argmax(Map[:, 6])]
        Position = Cell[:3]
        CellCount = int(Cell[6])

        #This is the TOTAL force on the cell, not each particle
        Force = 0.01

        Count = Propelled.add_propulsion(
            Position,
            (Force, 0.0, 0.0)
        )

        assert Count == CellCount
        assert Count > 0

        Normal.step(1)
        Propelled.step(1)

        #Compare velocities to isolate the effect of propulsion
        Difference = Propelled.velocities() - Normal.velocities()

        ParticleMass = Propelled.info()["particle_mass"]
        Expected = Force / Count / ParticleMass * small().TimeStep

        Changed = np.abs(Difference[:, 0]) > 1e-12

        assert int(np.sum(Changed)) == Count

        #Every affected particle should receive an equal share
        assert np.allclose(
            Difference[Changed, 0],
            Expected,
            rtol=1e-5,
            atol=1e-10
        )

        #Particles outside the selected cell should be unaffected
        assert np.allclose(
            Difference[~Changed],
            0.0,
            atol=1e-10
        )

        assert np.all(np.isfinite(Propelled.velocities()))


def test_propulsion_empty_cell():
    with WaveTank(small()) as sim:
        #Position outside the spatial grid should affect no particles
        Count = sim.add_propulsion(
            (100.0, 100.0, 100.0),
            (0.01, 0.0, 0.0)
        )

        assert Count == 0

        sim.step(1)
        assert np.isfinite(sim.velocities()).all()


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

        with pytest.raises(ValueError):
            sim.add_propulsion(
                (0.2, 0.2),
                (0.01, 0.0, 0.0)
            )

        with pytest.raises(ValueError):
            sim.add_propulsion(
                (0.2, 0.2, 0.05),
                (float("nan"), 0.0, 0.0)
            )

        sim.close()

        with pytest.raises(RuntimeError):
            sim.step(1)

        with pytest.raises(RuntimeError):
            sim.add_propulsion(
                (0.2, 0.2, 0.05),
                (0.01, 0.0, 0.0)
            )
