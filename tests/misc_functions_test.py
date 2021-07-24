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


def test_conjugate_EAN(fig8):
    ean = fig8._trace_field_numerical_root
    ean_conjugate = misc_functions.make_aan_conjugate(ean)
    assert abs(ean_conjugate(1000) - ean(1000).conjugate()) <= 10 ** (-1000)


@pytest.mark.parametrize("i", [0, 1, 2])
def test_conjugate_AAN(fig8, i):
    aan = fig8._approx_trace_field_gens[i]
    aan_conjugate = misc_functions.make_aan_conjugate(aan)
    assert abs(aan_conjugate(1000) - aan(1000).conjugate()) <= 10 ** (-1000)


@pytest.mark.parametrize("i", [0, 1, 2])
def test_conjugate_AAN2(i):
    mfld = ManifoldNT.ManifoldNT("7_7")
    aan = mfld._approx_invariant_trace_field_gens[i]
    aan_conjugate = misc_functions.make_aan_conjugate(aan)
    assert abs(aan_conjugate(1000) - aan(1000).conjugate()) <= 10 ** (-1000)


def test_conjugate_LAAN(fig8):
    laan = fig8._approx_trace_field_gens
    laan_conjugate = misc_functions.make_aan_conjugate(laan)
    pairs = zip(
        misc_functions.aan_iterator(laan), misc_functions.aan_iterator(laan_conjugate)
    )
    answers = [
        abs(elt[0](1000) - elt[1](1000).conjugate()) <= 10 ** (-1000) for elt in pairs
    ]
    assert answers == [True] * len(laan_conjugate)
