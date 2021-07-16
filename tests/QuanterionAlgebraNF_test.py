import itertools
from collections import Counter

import pytest
from sage.all import QQ, NumberField, var

import snappynt.QuaternionAlgebraNF as QuaternionAlgebraNF
from snappynt.field_isomorphisms import isomorphisms_between_number_fields

x = var("x")

# Fixtures: these are mostly particular algebras and fields that we use for tests.


# Fields
@pytest.fixture
def third_cyclo_field():
    return NumberField(x ** 2 + x + 1, "z")


@pytest.fixture
def cubic_field():
    # Discriminant equals -31.
    return NumberField(x ** 3 + x + 1, "z")


# Quaternion algebras and relevant ideals
@pytest.fixture
def div_alg_third_cyclo_field(third_cyclo_field):
    # Ramified at places above 2 and 3.
    field = third_cyclo_field
    z = field.gen()
    entries = field(z - 1), field(2)
    return QuaternionAlgebraNF.QuaternionAlgebraNF(field, *entries)


@pytest.fixture
def div_alg_third_cyclo_field_ramified_primes(third_cyclo_field):
    field = third_cyclo_field
    z = field.gen()
    dyadic_ideal = field.ideal(2)
    triadic_ideal = field.ideal(-2 * z - 1)
    return set([dyadic_ideal, triadic_ideal])


@pytest.fixture
def mat_alg_third_cyclo_field(third_cyclo_field):
    # Neither -z or -1 are squares, but the algebra is split.
    field = third_cyclo_field
    z = field.gen()
    entries = field(-z), field(-1)
    return QuaternionAlgebraNF.QuaternionAlgebraNF(field, *entries)


@pytest.fixture
def div_alg_cubic_field(cubic_field):
    # Ramified at a prime above 3 and at the real place.
    field = cubic_field
    z = field.gen()
    entries = field(z - 1), field(z ** 2 + 2 * z - 1)
    return QuaternionAlgebraNF.QuaternionAlgebraNF(field, *entries)


@pytest.fixture
def div_alg_cubic_field_ramified_primes(cubic_field):
    # It's a set even though there's only one element.
    field = cubic_field
    z = field.gen()
    return set([field.ideal(z - 1)])


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
def test_div_alg_third_cyclo_field_finite_residue_chars(div_alg_third_cyclo_field):
    residue_chars = div_alg_third_cyclo_field.ramified_residue_characteristics()
    assert residue_chars == Counter({2: 1, 3: 1})


def test_div_alg_third_cyclo_field_finite_ramification(
    div_alg_third_cyclo_field, div_alg_third_cyclo_field_ramified_primes
):
    # We throw in some other primes for fun (and to get some branch coverage).
    primes = {
        prime
        for prime in div_alg_third_cyclo_field.base_ring().ideal(7).prime_factors()
    }
    primes |= div_alg_third_cyclo_field_ramified_primes
    new_computation = set(
        [prime for prime in primes if div_alg_third_cyclo_field.is_ramified_at(prime)]
    )
    assert new_computation == div_alg_third_cyclo_field.ramified_finite_places()


def test_div_alg_third_cyclo_field_infinite_ramification(div_alg_third_cyclo_field):
    assert div_alg_third_cyclo_field.ramified_real_places() == set()


def test_div_alg_cubic_field_finite_residue_chars(div_alg_cubic_field):
    residue_chars = div_alg_cubic_field.ramified_residue_characteristics()
    assert residue_chars == Counter({3: 1})


def test_div_alg_cubic_field_finite_ramification(
    div_alg_cubic_field, div_alg_cubic_field_ramified_primes
):
    field = div_alg_cubic_field.base_ring()
    dyadic_primes = set(field.ideal(2).prime_factors())
    primes = dyadic_primes | div_alg_cubic_field_ramified_primes
    new_computation = set(
        [prime for prime in primes if div_alg_cubic_field.is_ramified_at(prime)]
    )
    assert new_computation == div_alg_cubic_field.ramified_finite_places()


def test_div_alg_cubic_field_infinite_ramification(div_alg_cubic_field):
    field = div_alg_cubic_field.base_ring()
    place = field.real_places()[0]
    assert div_alg_cubic_field.is_ramified_at(place)


def test_division_algebras(division_algebras):
    failing_algebras = set()
    for algebra in division_algebras:
        if not algebra.is_division_algebra():
            failing_algebras.add(algebra)
    assert failing_algebras == set()


def test_matrix_algebras(matrix_algebras):
    failing_algebras = set()
    for algebra in matrix_algebras:
        if not algebra.is_matrix_ring():
            failing_algebras.add(algebra)
    assert failing_algebras == set()


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


def test_ramification_string_works(all_algebras):
    # This just checks that ramification string runs without errors and returns a
    # string. We do some sorcery with keyword arguments to make coverage's branch
    # coverage tool happy.
    failing_algebras = set()
    for algebra in all_algebras:
        algebra.ramified_places()
        keywords = (
            "full_finite_ramification",
            "full_real_ramification",
            "show_hilbert_symbol",
            "show_field_data",
        )
        bool_tups = itertools.product((True, False), repeat=4)
        keyword_dicts = [dict(zip(keywords, bool_tup)) for bool_tup in bool_tups]
        for keyword_dict in keyword_dicts:
            s = algebra.ramification_string(**keyword_dict)
            if not isinstance(s, str):
                failing_algebras.add(algebra)
    assert failing_algebras == set()


def test_error_raised_for_rationals():
    with pytest.raises(NotImplementedError):
        QuaternionAlgebraNF.QuaternionAlgebraNF(QQ, 1, 1)


# Test for error raised when making new QA.


def test_bad_isomorphism(div_alg_third_cyclo_field):
    old_field = div_alg_third_cyclo_field.base_ring()
    new_field = NumberField(x ** 2 + 3)
    iso = isomorphisms_between_number_fields(new_field, old_field)[0]
    with pytest.raises(ValueError):
        div_alg_third_cyclo_field.new_QA_via_field_isomorphism(iso)


# Tests for is_isomorphic


def test_isomorphic_to_self(all_algebras):
    failing_algebras = set()
    for algebra in all_algebras:
        if not algebra.is_isomorphic(algebra):
            failing_algebras.add(algebra)
    assert failing_algebras == set()


def test_same_field_different_real_ram(div_alg_cubic_field, mat_alg_cubic_field):
    # Testing the branch where the algebras are distinguished by real ramification.
    # This is a priori a bit delicate since it depends on the order that things
    # are tested in is_isomorphic.
    assert not div_alg_cubic_field.is_isomorphic(mat_alg_cubic_field)


def test_same_field_different_nondyadic_ram(div_alg_cubic_field, cubic_field):
    # Have to cook up a new algebra here. Fortunately, (-1,-1)_K will be ramified
    # at the real place and the dyadic place. div_alg_cubic_field has triadic
    # ramification, so they will be distinguished there.
    new_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(cubic_field, -1, -1)
    assert not new_alg.is_isomorphic(div_alg_cubic_field)


def test_same_field_different_dyadic_ram():
    # Need at least two dyadic places for all the other ramification to match up.
    field = NumberField(x ** 2 + 7, "z")
    z = field.gen()
    matrix_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(field, 1, 1)
    division_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(
        field, field((1 / 2)) * z + field((1 / 2)), -1
    )
    assert not matrix_alg.is_isomorphic(division_alg)


def test_same_field_isomorphic_algs(div_alg_cubic_field, cubic_field):
    z = cubic_field.gen()
    new_invariants = [z ** 2 * inv for inv in div_alg_cubic_field.invariants()]
    new_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(cubic_field, *new_invariants)
    assert new_alg.is_isomorphic(div_alg_cubic_field)


# Tests for same_ramification_via_isomorphism


def test_bad_isomorphism_for_via_iso(div_alg_third_cyclo_field):
    old_field = div_alg_third_cyclo_field.base_ring()
    new_field = NumberField(x ** 2 + 3)
    iso = isomorphisms_between_number_fields(new_field, old_field)[0]
    with pytest.raises(ValueError):
        div_alg_third_cyclo_field.same_ramification_via_isomorphism(iso)


def test_different_field_different_real_ram(div_alg_cubic_field):
    # Real ramification should be the first thing checked.
    old_field = div_alg_cubic_field.base_ring()
    new_field = NumberField(x ** 3 + x + 1, "t")
    isos = isomorphisms_between_number_fields(old_field, new_field)
    matrix_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(new_field, 1, 1)
    results = [matrix_alg.same_ramification_via_isomorphism(iso) for iso in isos]
    assert True not in results


"""
def test_different_field_different_nondyadic_ram_computed(div_alg_cubic_field):
    new_field = NumberField(x ** 3 + x + 1, "t")
    new_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(new_field, -1, -1)
    new_alg.ramified_places()
    div_alg_cubic_field.ramified_places()
    assert not new_alg.is_isomorphic(div_alg_cubic_field)


def test_different_field_different_nondyadic_ram(div_alg_cubic_field):
    new_field = NumberField(x ** 3 + x + 1, "t")
    new_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(new_field, -1, -1)
    assert not new_alg.is_isomorphic(div_alg_cubic_field)


def test_different_field_different_dyadic_ram_computed():
    field = NumberField(x ** 2 + 7, "z")
    z = field.gen()
    field2 = NumberField(x ** 2 + 7, "t")
    matrix_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(field2, 1, 1)
    division_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(
        field, field((1 / 2)) * z + field((1 / 2)), -1
    )
    matrix_alg.ramified_places()
    division_alg.ramified_places()
    assert not matrix_alg.is_isomorphic(division_alg)


def test_different_fields_computed(mat_alg_cubic_field, mat_alg_third_cyclo_field):
    mat_alg_cubic_field.ramified_places()
    mat_alg_third_cyclo_field.ramified_places()
    assert not mat_alg_cubic_field.is_isomorphic(mat_alg_third_cyclo_field)


def test_different_field_different_dyadic_ram():
    field = NumberField(x ** 2 + 7, "z")
    z = field.gen()
    field2 = NumberField(x ** 2 + 7, "t")
    matrix_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(field2, 1, 1)
    division_alg = QuaternionAlgebraNF.QuaternionAlgebraNF(
        field, field((1 / 2)) * z + field((1 / 2)), -1
    )
    assert not matrix_alg.is_isomorphic(division_alg)
"""
