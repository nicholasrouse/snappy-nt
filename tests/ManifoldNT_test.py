import pytest

import snappynt.ManifoldNT as ManifoldNT


@pytest.fixture
def fig_eight_computed():
    mfld = ManifoldNT.ManifoldNT("4_1")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    return mfld


@pytest.fixture
def m010minusOneTwo_computed():
    mfld = ManifoldNT.ManifoldNT("m010(-1,2)")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    return mfld


def test_fix_names_error():
    with pytest.raises(ValueError):
        ManifoldNT.fix_names("bad name")


def test_fig_eight_homology(fig_eight_computed):
    mfld = fig_eight_computed
    assert mfld.homology_two_rank() == 0


def test_m010minusOneTwo_homology(m010minusOneTwo_computed):
    # Its trace field and invariant trace field are different.
    mfld = m010minusOneTwo_computed
    assert mfld.homology_two_rank() == 1


def test_str_magic_method(fig_eight_computed):
    assert fig_eight_computed.__str__() == "4_1(0,0)"


def test_repr_magic_method(fig_eight_computed):
    assert fig_eight_computed.__repr__() == "4_1(0,0)"
