import pytest
from sage.all import NumberField, var

import snappynt.QuaternionAlgebraNF as QuaternionAlgebraNF

# Fixtures: these are mostly particular algebras and fields that we use for tests.


@pytest.fixture
def QQ_adjoin_sqrt_minus_two():
    x = var("x")
    field = NumberField(x ** 2 + 2, "z")
    return field


@pytest.fixture
def div_alg_over_QQ_adj_sqrt_minus_two(QQ_adjoin_sqrt_minus_two):
    field = QQ_adjoin_sqrt_minus_two
    z = field.gen()
    return QuaternionAlgebraNF.QuaternionAlgebraNF(field, 1 + z, z)


def test_different_names(QQ_adjoin_sqrt_minus_two):
    field = QQ_adjoin_sqrt_minus_two
    z = field.gen()
    alg1 = QuaternionAlgebraNF.QuaternionAlgebraNF(
        field, 1 + z, z, names="cat,dog,fish"
    )
    alg2 = QuaternionAlgebraNF.QuaternionAlgebraNF(field, 1 + z, z)
    assert alg1.is_isomorphic(alg2)


def test_delayed_ramification():
    pass
