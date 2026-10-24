import pathlib

import pytest

from flowcy import synth, load_sample


@pytest.fixture(scope="session")
def study(tmp_path_factory):
    folder = tmp_path_factory.mktemp("studie")
    truth = synth.make_study(folder, n_per_group=4, n=30_000)
    return folder, truth


@pytest.fixture(scope="session")
def sample_and_truth(study):
    folder, truth = study
    s = load_sample(folder / "K01_PBMC.fcs")
    return s, truth["K01_PBMC.fcs"]


@pytest.fixture()
def fresh_sample(study):
    folder, truth = study
    return load_sample(pathlib.Path(folder) / "K01_PBMC.fcs"), truth["K01_PBMC.fcs"]
