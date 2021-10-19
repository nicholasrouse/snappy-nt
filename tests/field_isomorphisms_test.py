import pytest
from sage.all import CC, I, NumberField, var


@pytest.fixture
def QQ_adjoin_sqrt_minus_two():
    x = var("x")
    return NumberField(x ** 2 + 2, "z")


@pytest.fixture
def QQ_adjoin_sqrt_minus_two_with_embedding():
    x = var("x")
    return NumberField(x ** 2 + 2, "z", embedding=CC(1.4 * I))
