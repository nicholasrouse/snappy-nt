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
    assert not mfld._has_two_torsion_in_homology()


def test_m010minusOneTwo_homology(m010minusOneTwo_computed):
    # Its trace field and invariant trace field are different.
    mfld = m010minusOneTwo_computed
    assert mfld._has_two_torsion_in_homology()


def test_str_magic_method(fig_eight_computed):
    assert fig_eight_computed.__str__() == "4_1(0,0)"


def test_repr_magic_method(fig_eight_computed):
    assert fig_eight_computed.__repr__() == "4_1(0,0)"


def test_eight_seventeen_trace_field():
    mfld = ManifoldNT.ManifoldNT("8_17")
    while mfld._trace_field is None:
        mfld.trace_field()
    assert mfld._trace_field.degree() == 18


def test_eight_seventeen_computation():
    mfld = ManifoldNT.ManifoldNT("8_17")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    mfld.next_prec_and_degree("trace field")
    assert mfld._arithmetic_invariants_known()


def test_nine_fourteen_computation():
    # This trace field is large enough to be outside the defaults.
    mfld = ManifoldNT.ManifoldNT("9_14")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    assert mfld._arithmetic_invariants_known()


def test_m010minusOneTwo_computation():
    mfld = ManifoldNT.ManifoldNT("m010(-1,2)")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    mfld.next_prec_and_degree("trace field")
    mfld.next_prec_and_degree("invariant trace field")
    mfld.next_prec_and_degree("quaternion algebra")
    mfld.next_prec_and_degree("invariant quaternion algebra")
    assert mfld._arithmetic_invariants_known()
