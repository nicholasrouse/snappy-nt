import pytest
from sage.all import CC, ComplexField, I, NumberField, var

from snappynt import ManifoldNT, misc_functions

# from snappy.snap.find_field import (
#    ApproximateAlgebraicNumber,
#    ExactAlgebraicNumber,
#    ListOfApproximateAlgebraicNumbers,
# )


@pytest.fixture
def fig8():
    mfld = ManifoldNT.ManifoldNT("4_1")
    while not mfld._arithmetic_invariants_known():
        mfld.compute_arithmetic_invariants()
    return mfld


def test_commutator_of_word1():
    word1 = "abAB"
    word2 = "AABB"
    answer = "abABAABBbaBAbbaa"
    assert misc_functions.commutator_of_words(word1, word2) == answer


def test_commutator_of_words2():
    word1 = "abAB"
    word2 = ""
    answer = "abABbaBA"
    assert misc_functions.commutator_of_words(word1, word2) == answer


def test_aan_exception(fig8):
    """
    This *passes* as long as the base ListOfApproximateAlgebraicNumbers class raises
    an exception during iteration.
    """
    gens = fig8._approx_trace_field_gens
    with pytest.raises(IndexError):
        for item in gens:
            item(1000)


def test_aan_iterator(fig8):
    gens = fig8._approx_trace_field_gens
    iterator = misc_functions.aan_iterator(gens)
    s = set()
    for item in iterator:
        s.add(item)
    assert s == set(gens.list())


def test_make_aan():
    x = var("x")
    poly = x ** 2 - 2
    root = CC(1.4)
    aan = misc_functions.make_aan(poly, root)
    min_poly = aan.min_polynomial()
    assert min_poly == poly


def test_make_aan2():
    x = var("x")
    poly = x ** 2 - 2
    root = CC(1.4)
    aan = misc_functions.make_aan(poly, root)
    assert abs(aan(1000) - ComplexField(1000)(1.414)) <= 10 ** (-3)


def test_conjugate_field():
    x = var("x")
    field = NumberField(x ** 2 + 2, "z", embedding=CC(1.4 * I))
    new_field = misc_functions.conjugate_field(field)
    assert abs(CC(-1.414 * I) - new_field.gen_embedding()) <= 10 ** (-3)
