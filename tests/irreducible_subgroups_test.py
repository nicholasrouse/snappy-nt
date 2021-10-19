import pytest
import snappy
from sage.all import RR, I

import snappynt.irreducible_subgroups as irreducible_subgroups


@pytest.fixture
def fig8_group():
    G = snappy.Manifold("4_1").polished_holonomy(bits_prec=1000)
    return G


@pytest.fixture
def rank2_free_group_words():
    G = irreducible_subgroups.enumerate_words(2)
    return G


def test_within_epsilon():
    assert not irreducible_subgroups.within_epsilon(1, 2)


def test_within_epsilon2():
    assert not irreducible_subgroups.within_epsilon(1, RR(1.1))


def test_within_epsilon3():
    assert not irreducible_subgroups.within_epsilon(RR(1.1), 1)


def test_within_epsilon4():
    assert not irreducible_subgroups.within_epsilon(I, 1)


def test_within_epsilon5():
    assert irreducible_subgroups.within_epsilon(1, 1)


def test_within_epsilon6():
    a = RR(1)
    b = RR(1) + RR(10 ** (-53) / 2)
    assert irreducible_subgroups.within_epsilon(a, b)


def test_is_parabolic(fig8_group):
    # peripheral_curves() returns a list of tuples, one for each cusp.
    mer, long = (fig8_group(elt) for elt in fig8_group.peripheral_curves()[0])
    assert irreducible_subgroups.is_parabolic(
        mer, epsilon_coefficient=100
    ) and irreducible_subgroups.is_parabolic(long, epsilon_coefficient=100)


def test_is_parabolic2(fig8_group):
    a, b = (fig8_group(elt) for elt in fig8_group.generators())
    assert not (
        irreducible_subgroups.is_parabolic(a) or irreducible_subgroups.is_parabolic(b)
    )


def test_generate_reducible_subgroup(fig8_group):
    a = fig8_group(fig8_group.generators()[0])
    assert irreducible_subgroups.generate_reducible_subgroup(a, a * a)


def test_generate_reducible_subgroup2(fig8_group):
    gens = (fig8_group(elt) for elt in fig8_group.generators())
    assert not irreducible_subgroups.generate_reducible_subgroup(*gens)


def test_enumerate_words(rank2_free_group_words):
    # There are 52 words in Tietze notation with length less than 3.
    elts = list()
    for _ in range(52):
        elts.append(next(rank2_free_group_words))
    assert (
        len(elts) == len(set(tuple(elt) for elt in elts))  # No duplicates
        and False not in (len(item) <= 3 for item in elts)
        and len(next(rank2_free_group_words)) == 4
    )


def test_enumerate_group_elements(fig8_group):
    generator = irreducible_subgroups.enumerate_group_elements(fig8_group, as_word=True)
    elts = list()
    for _ in range(52):
        elts.append(next(generator))
    assert (
        len(elts) == len(set(elts))
        and False not in (len(item) <= 3 for item in elts)
        and len(next(generator)) == 4
    )


def test_enumerate_group_elements2(fig8_group):
    words_gen = irreducible_subgroups.enumerate_group_elements(fig8_group, as_word=True)
    mats_gen = irreducible_subgroups.enumerate_group_elements(fig8_group, as_word=False)
    words = list()
    mats = list()
    for _ in range(52):
        words.append(next(words_gen))
        mats.append(next(mats_gen))
    assert len(set(words)) == len(
        mats
    ) and False not in (  # mutable matrices are unhashable
        fig8_group(word) == mat for word, mat in zip(words, mats)
    )


def test_find_hilbert_symbol_words(fig8_group):
    words = irreducible_subgroups.find_hilbert_symbol_words(fig8_group)
    assert 0 not in words


def test_find_hilbert_symbol_words2(fig8_group):
    def fake_gens():
        return ["a" * i for i in range(1, 10)] + ["b"]

    fig8_group.generators = fake_gens
    words = irreducible_subgroups.find_hilbert_symbol_words(fig8_group)
    assert 0 not in words


def test_find_hilbert_symbol_words3(fig8_group):
    mer, _ = fig8_group.peripheral_curves()[0]

    def fake_gens():
        return [mer, "a"]

    fig8_group.generators = fake_gens
    words = irreducible_subgroups.find_hilbert_symbol_words(fig8_group)
    assert 0 not in words
