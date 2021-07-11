import itertools
from collections import Counter

import pytest
from sage.all import NumberField, var

import snappynt.QuaternionAlgebraNF as QuaternionAlgebraNF

# Fixtures: these are mostly particular algebras and fields that we use for tests.


# Fields
@pytest.fixture
def third_cyclo_field():
    x = var("x")
    return NumberField(x ** 2 + x + 1, "z")


@pytest.fixture
def cubic_field():
    # Discriminant equals -31.
    x = var("x")
    return NumberField(x ** 3 + x + 1, "z")


# Quaternion algebras
@pytest.fixture
def div_alg_third_cyclo_field(third_cyclo_field):
    # Ramified at places above 2 and 3.
    field = third_cyclo_field
    z = field.gen()
    entries = field(z - 1), field(2)
    return QuaternionAlgebraNF.QuaternionAlgebraNF(field, *entries)


@pytest.fixture
def mat_alg_third_cyclo_field(third_cyclo_field):
    # Neither -z or -1 are squares, but the algebra is split.
    field = third_cyclo_field
    z = field.gen()
    entries = field(-z), field(-1)
    return QuaternionAlgebraNF.QuaternionAlgebraNF(field, *entries)


@pytest.fixture
def div_alg_cubic_field(cubic_field):
    # Ramified at primes above 2, 3, and 23 and at the real place.
    field = cubic_field
    z = field.gen()
    entries = field(z - 1), field(z ** 2 + 2 * z - 1)
    return QuaternionAlgebraNF.QuaternionAlgebraNF(field, *entries)


@pytest.fixture
def mat_alg_cubic_field(cubic_field):
    # z is not a square, but (z, 1-z) is always split when z != 0,1.
    field = cubic_field
    z = field.gen()
    entries = field(z), field(1 - z)
    return QuaternionAlgebraNF.QuaternionAlgebraNF(field, *entries)


# Fixtures that return iterables. Used for testing behavior that should be uniform
# across instances. Unfortunately, we have to keep these updated by hand if we add more
# in the future.
@pytest.fixture
def matrix_algebras(mat_alg_third_cyclo_field, mat_alg_cubic_field):
    return (mat_alg_third_cyclo_field, mat_alg_cubic_field)


@pytest.fixture
def division_algebras(div_alg_third_cyclo_field, div_alg_cubic_field):
    return (div_alg_third_cyclo_field, div_alg_cubic_field)


@pytest.fixture
def all_algebras(matrix_algebras, division_algebras):
    return itertools.chain(matrix_algebras, division_algebras)


# Tests
def test_div_alg_third_cyclo_field_finite_ramification(div_alg_third_cyclo_field):
    residue_chars = div_alg_third_cyclo_field.ramified_residue_characteristics()
    assert residue_chars == Counter({2: 1, 3: 1})


def test_div_alg_third_cyclo_field_infinite_ramification(div_alg_third_cyclo_field):
    assert div_alg_third_cyclo_field.ramified_real_places() == set()


def test_div_alg_cubic_field_finite_ramification(div_alg_cubic_field):
    residue_chars = div_alg_cubic_field.ramified_residue_characteristics()
    assert residue_chars == Counter({2: 1, 3: 1, 23: 1})


def test_div_alg_cubic_field_infinite_ramification(div_alg_cubic_field):
    assert len(div_alg_cubic_field.ramified_real_places()) == 1


def test_real_ramification_flag(all_algebras):
    failing_algebras = set()
    for algebra in all_algebras:
        algebra.ramified_real_places()
        if not algebra._ramified_real_places_known:
            failing_algebras.add(algebra)
    assert failing_algebras == set()


def test_dyadic_ramification_flag(all_algebras):
    failing_algebras = set()
    for algebra in all_algebras:
        algebra.ramified_dyadic_places()
        if not algebra._ramified_dyadic_places_known:
            failing_algebras.add(algebra)
    assert failing_algebras == set()


def test_nondyadic_ramification_flag(all_algebras):
    failing_algebras = set()
    for algebra in all_algebras:
        algebra.ramified_nondyadic_places()
        if not algebra._ramified_nondyadic_places_known:
            failing_algebras.add(algebra)
    assert failing_algebras == set()


def test_dyadic_dict_right_size(all_algebras):
    failing_algebras = set()
    for algebra in all_algebras:
        dyadic_primes = set(algebra.base_ring().ideal(2).prime_factors())
        algebra.ramified_dyadic_places()
        if dyadic_primes != set(algebra._ramified_dyadic_places_dict.keys()):
            failing_algebras.add(algebra)
    assert failing_algebras == set()


def test_double_computation(all_algebras):
    failing_algebras = set()
    for algebra in all_algebras:
        first_pass = algebra.ramified_real_places() | algebra.ramified_finite_places()
        second_pass = algebra.ramified_real_places() | algebra.ramified_finite_places()
        if first_pass != second_pass:
            failing_algebras.add(algebra)
    assert failing_algebras == set()


def test_different_names(all_algebras):
    failing_algebras = set()
    for algebra in all_algebras:
        field = algebra.base_ring()
        a, b = algebra.invariants()
        new_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(
            field, a, b, names="cat,dog,fish"
        )
        if not new_alg.is_isomorphic(algebra):
            failing_algebras.add(algebra)
    assert failing_algebras == set()
