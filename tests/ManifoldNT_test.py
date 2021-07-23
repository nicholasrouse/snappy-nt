import pytest
from sage.all import NumberField, var

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


@pytest.fixture
def nine_fourteen_trace_field():
    x = var("x")
    poly = (
        x ** 18
        - 3 * x ** 17
        - 5 * x ** 16
        - 4 * x ** 15
        + 192 * x ** 14
        - 2375 * x ** 13
        + 14845 * x ** 12
        - 36312 * x ** 11
        + 36785 * x ** 10
        - 213737 * x ** 9
        + 1454725 * x ** 8
        - 4869252 * x ** 7
        + 10456460 * x ** 6
        - 16339634 * x ** 5
        + 19631702 * x ** 4
        - 21382208 * x ** 3
        + 24861005 * x ** 2
        - 22425163 * x
        + 9093699
    )
    field = NumberField(poly, "z")
    return field


def test_fix_names_error():
    with pytest.raises(ValueError):
        ManifoldNT.fix_names("bad name")


def test_fig_eight_homology(fig_eight_computed):
    mfld = fig_eight_computed
    assert not mfld._has_two_torsion_in_homology()


def test_fig_eight_arithmetic():
    mfld = ManifoldNT.ManifoldNT("4_1")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    assert mfld.is_arithmetic()


def test_fig_eight_arithmetic2():
    mfld = ManifoldNT.ManifoldNT("4_1")
    # is_arithmetic() doesn't compute anything on its own.
    # We should perhaps change this behavior though.
    assert mfld.is_arithmetic() is None


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


def test_eight_seventeen_not_arithmetic():
    mfld = ManifoldNT.ManifoldNT("8_17")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    assert not mfld.is_arithmetic()


def test_eight_seventeen_computation():
    mfld = ManifoldNT.ManifoldNT("8_17")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    mfld.next_prec_and_degree("trace field")
    assert mfld._arithmetic_invariants_known()


def test_nine_fourteen_computation():
    # This trace field is large enough to be outside the defaults.
    mfld = ManifoldNT.ManifoldNT("9_14")
    mfld.compute_arithmetic_invariants(prec=100, degree=5)
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    mfld.compute_arithmetic_invariants()
    assert mfld._arithmetic_invariants_known()


def test_nine_fourteen_itf_first(nine_fourteen_trace_field):
    # We tell it the invariant trace field to check the logic for knowing the this, but
    # not the trace field itself.
    mfld = ManifoldNT.ManifoldNT("9_14")
    mfld._invariant_trace_field = nine_fourteen_trace_field
    mfld.trace_field(prec=100, degree=5)
    while mfld._trace_field is None:
        mfld.trace_field()
    assert mfld._trace_field.is_isomorphic(mfld._invariant_trace_field)


def test_nine_fourteen_itf_first2():
    mfld = ManifoldNT.ManifoldNT("9_14")
    while mfld._invariant_trace_field is None:
        mfld.invariant_trace_field()
    mfld.trace_field(prec=100, degree=5)
    while mfld._trace_field is None:
        mfld.trace_field()
    assert mfld._trace_field.is_isomorphic(mfld._invariant_trace_field)


def test_nine_fourteen_tf_first(nine_fourteen_trace_field):
    mfld = ManifoldNT.ManifoldNT("9_14")
    mfld._trace_field = nine_fourteen_trace_field
    mfld.invariant_trace_field(prec=100, degree=5)
    while mfld._invariant_trace_field is None:
        mfld.invariant_trace_field()
    mfld.invariant_trace_field()
    assert mfld._trace_field.is_isomorphic(mfld._invariant_trace_field)


def test_nine_fourteen_denoms_without_tf():
    mfld = ManifoldNT.ManifoldNT("9_14")
    mfld.denominators()
    while mfld._denominators is None:
        mfld.denominators()
    assert mfld._denominators is not None


def test_m141twoThree_computation():
    mfld = ManifoldNT.ManifoldNT("m141(2,3)")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    assert mfld._arithmetic_invariants_known()


def test_m141twoThree_qa_small_prec():
    mfld = ManifoldNT.ManifoldNT("m141(2,3)")
    while not mfld._trace_field:
        mfld.trace_field()
    mfld.quaternion_algebra(prec=100)
    while not mfld._quaternion_algebra:
        mfld.quaternion_algebra()
    mfld.quaternion_algebra()
    assert mfld._quaternion_algebra is not None


def test_m010minusOneTwo_computation():
    mfld = ManifoldNT.ManifoldNT("m010(-1,2)")
    mfld.invariant_trace_field()
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    mfld.next_prec_and_degree("trace field")
    mfld.next_prec_and_degree("invariant trace field")
    mfld.next_prec_and_degree("quaternion algebra")
    mfld.next_prec_and_degree("invariant quaternion algebra")
    mfld.invariant_quaternion_algebra()
    assert mfld._arithmetic_invariants_known()


def test_m010minusOneTwo_iqa_computation():
    mfld = ManifoldNT.ManifoldNT("m010(-1,2)")
    mfld.invariant_quaternion_algebra(prec=25)  # low enough not to compute.
    while not mfld._invariant_quaternion_algebra:
        mfld.invariant_quaternion_algebra()
    mfld.invariant_quaternion_algebra()
    assert mfld._invariant_quaternion_algebra is not None


def test_m137_computation():
    mfld = ManifoldNT.ManifoldNT("m137")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    mfld.compute_arithmetic_invariants()
    assert mfld._denominators is not None and mfld._denominators != set()


def test_m137_denominator_residue_chars():
    mfld = ManifoldNT.ManifoldNT("m137")
    while mfld._denominator_residue_characteristics is None:
        mfld.denominator_residue_characteristics()
    mfld.denominator_residue_characteristics()
    assert 2 in mfld._denominator_residue_characteristics


def test_m137_denominator_residue_chars2():
    # This one simulates computing the residue characteristics from the ideals. Usually
    # the residue characteristics are computed in denominators().
    mfld = ManifoldNT.ManifoldNT("m137")
    while mfld._denominators is None:
        mfld.denominators()
    mfld._denominator_residue_characteristics = None
    mfld.denominator_residue_characteristics()
    assert 2 in mfld._denominator_residue_characteristics


def test_dehn_fill():
    mfld = ManifoldNT.ManifoldNT("7_4")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    mfld.dehn_fill(())
    assert not mfld._arithmetic_invariants_known()


# Testing _isomorphic_quaternion_algebras


@pytest.mark.skip(reason="Broken right now")
def test_fields_missing_qa():
    mfld1 = ManifoldNT.ManifoldNT("9_14")
    mfld2 = ManifoldNT.ManifoldNT("7_4")
    with pytest.raises(RuntimeError):
        mfld1._isomorphic_quaternion_algebras(mfld2)


@pytest.mark.skip(reason="Broken right now")
def test_fields_missing_iqa():
    mfld1 = ManifoldNT.ManifoldNT("9_14")
    mfld2 = ManifoldNT.ManifoldNT("7_4")
    with pytest.raises(RuntimeError):
        mfld1._isomorphic_quaternion_algebras(mfld2, _invariant_qa=True)


@pytest.mark.skip(reason="Broken right now")
def test_isomorphic_fields_matrix_algebras():
    # The fields are abstractly isomorphic and the algebras are both matrix.
    mfld1 = ManifoldNT.ManifoldNT("6_1")
    mfld2 = ManifoldNT.ManifoldNT("7_7")
    while not mfld1._arithmetic_invariants_known():
        mfld1.compute_arithmetic_invariants()
    while not mfld2._arithmetic_invariants_known():
        mfld2.compute_arithmetic_invariants()
    assert mfld1._isomorphic_quaternion_algebras(mfld2)


@pytest.mark.skip(reason="Broken right now")
def test_nonisomorphic_fields_iqas():
    # The fields are different
    mfld1 = ManifoldNT.ManifoldNT("m003(-2,3)")
    mfld2 = ManifoldNT.ManifoldNT("m003(-3,1)")
    while not mfld1._arithmetic_invariants_known():
        mfld1.compute_arithmetic_invariants()
    while not mfld2._arithmetic_invariants_known():
        mfld2.compute_arithmetic_invariants()
    assert not mfld1._isomorphic_quaternion_algebras(mfld2, _invariant_qa=True)
