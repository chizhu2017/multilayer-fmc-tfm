import importlib.util
from pathlib import Path
import numpy as np
import pytest
from scipy.signal import hilbert
from tfm import Config, reconstruct
from tfm.data import Dataset
from tfm.reconstruction import interpolate, time_zero
from tfm.travel import travel_tables

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("generate_demo", ROOT/"tools/generate_demo.py")
demo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(demo)


def test_interpolation_one_based_and_final_sample():
    np.testing.assert_array_equal(interpolate(np.array([10, 20, 40], np.float32),
        np.array([.99, 1, 1.5, 2.5, 3, 4, np.nan], np.float32)), [0, 10, 15, 30, 0, 0, 0])


def test_lut_against_unrestricted_absolute_interface_search():
    c = Config(nx=11, h1=4, h2=3, h3=3, x_loc=2, gap=2)
    lut2, lut3 = travel_tables(c, 10)
    q = np.arange(1, 12, dtype=np.float64)
    for probe, target in [(1,11), (9,2), (5,5)]:
        d = abs(probe-target)
        first = np.hypot(q-probe, c.h1)/c.v1
        pre = np.min(first[:,None] + np.hypot(q[:,None]-q[None,:],c.h2)/c.v2, axis=0)
        for dz in range(c.h2):
            expected = np.min(first + np.hypot(target-q,dz)/c.v2)
            np.testing.assert_allclose(lut2[dz,d],expected,rtol=2e-7)
        for dz in range(c.h3+1):
            expected = np.min(pre + np.hypot(target-q,dz)/c.v3)
            np.testing.assert_allclose(lut3[dz,d],expected,rtol=2e-7)


def test_signed_first_peak_and_validation(tmp_path):
    t = np.arange(5)*1e-6
    np.savetxt(tmp_path/"1.txt",np.column_stack([t,[-100,2,2,1,0]]))
    ds = Dataset(tmp_path)
    assert ds.peak_indices()[0,0] == 2
    (tmp_path/"extra.txt").write_text("0 1\n1 2")
    with pytest.raises(ValueError,match="exactly"):
        Dataset(tmp_path)


def test_nonuniform_time_rejected(tmp_path):
    np.savetxt(tmp_path/"1.txt",np.column_stack([[0.,1e-6,3e-6],[1,2,3]]))
    with pytest.raises(ValueError,match="uniform"):
        Dataset(tmp_path)


def test_invalid_time_zero():
    with pytest.raises(ValueError,match="time-zero"):
        time_zero(np.zeros((1,1),np.float32),np.array([1],np.float32),np.float32(1e-6),Config())


def test_synthetic_end_to_end_and_target_localization(tmp_path):
    demo.generate(tmp_path,ROOT/"config/demo.json")
    c = Config.load(ROOT/"config/demo.json")
    result = reconstruct(tmp_path,c)
    assert result.image.shape == (32,32)
    assert result.image.dtype == np.float32
    np.testing.assert_allclose(result.envelope,np.abs(hilbert(result.image,axis=0)),rtol=1e-6)
    # Strongest target in layer 3 should focus near (x,z)=(16,27) mm.
    layer3 = result.envelope[22:]
    z,x = np.unravel_index(np.argmax(layer3),layer3.shape)
    assert abs((x+1)-16)<=1 and abs((z+23)-27)<=1


@pytest.mark.parametrize("kwargs",[{"h1":1},{"gap":0},{"v1":float("nan")},{"nx":3.5}])
def test_bad_config(kwargs):
    with pytest.raises(ValueError):
        Config(**kwargs)
