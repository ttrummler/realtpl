from contextlib import contextmanager
import glob
import os
import pandas as pd
import shutil
import subprocess
import tempfile
from unittest import TestCase
import yaml


def test_cyclopentane():
    config_file = os.path.join(
        os.path.dirname(__file__), "config_test_Cyclopentane.yaml"
    )
    actual_dfs = _run(config_file)
    expected_dfs = _load_expected(config_file)

    _compare(actual_dfs, expected_dfs)
    return


def test_n_dodecane():
    config_file = os.path.join(
        os.path.dirname(__file__), "config_test_nDodecane.yaml"
    )
    actual_dfs = _run(config_file)
    expected_dfs = _load_expected(config_file)

    _compare(actual_dfs, expected_dfs)
    return

def test_hydrogen():
    config_file = os.path.join(
        os.path.dirname(__file__), "config_test_hydrogen.yaml"
    )
    actual_dfs = _run(config_file)
    expected_dfs = _load_expected(config_file)

    _compare(actual_dfs, expected_dfs)
    return

def test_air():
    config_file = os.path.join(
        os.path.dirname(__file__), "config_test_air.yaml"
    )
    actual_dfs = _run(config_file)
    expected_dfs = _load_expected(config_file)

    _compare(actual_dfs, expected_dfs)
    return

def _run(config_file):
    with open(config_file, 'r') as f:
        cfg = yaml.safe_load(f)

    dfs = {}
    with _temp_dir() as tmpdir:
        shutil.copy(config_file, tmpdir)
        subprocess.run(
            ["realtpl", "--config-file", os.path.basename(config_file)],
            cwd=tmpdir,
            check=True
        )
        for x in glob.glob(os.path.join(tmpdir, cfg['output_dir'],
                                        cfg['fluid_name'], 'data', '*.csv')):
            dfs[os.path.basename(x)] = pd.read_csv(x, sep='\t')
    return dfs


@contextmanager
def _temp_dir():
    tmpdir = tempfile.mkdtemp()
    try:
        yield tmpdir
    finally:
        shutil.rmtree(tmpdir)


def _load_expected(config_file):
    with open(config_file, 'r') as f:
        cfg = yaml.safe_load(f)
    dfs = {}
    for x in glob.glob(os.path.join(os.path.dirname(__file__), 'reference',
                                    cfg['fluid_name'], 'data', '*.csv')):
        dfs[os.path.basename(x)] = pd.read_csv(x, sep='\t')
    return dfs


def _compare(actual_dfs, expected_dfs):
    TestCase().assertCountEqual(actual_dfs.keys(), expected_dfs.keys())

    for key in expected_dfs:
        print(key)
        pd.testing.assert_frame_equal(
            actual_dfs[key],
            expected_dfs[key],
        )
    return
