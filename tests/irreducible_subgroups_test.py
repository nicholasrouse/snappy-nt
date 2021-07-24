import pytest
import snappy
from sage.all import RR, I

import snappynt.irreducible_subgroups as irreducible_subgroups


@pytest.fixture
def fig8_group():
    G = snappy.Manifold("4_1").polished_holonomy(bits_prec=1000)
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
