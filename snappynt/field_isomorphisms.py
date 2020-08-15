"""
This module includes functions to express one algebraic number given in terms of powers
of a root of an irreducible polynomial (over the rational numbers) in terms of powers
of roots of another polynomial.

Sage theoretically has this functionality, but as of 9.1, it doesn't actually compute
the isomorphism correctly when one of the defining polynomials is not monic and
integral.  The reason is that sage's implementation converts to PARI using 
pari_polynomial which is always monic and integral.  That said, sage can factor
polynomials over number fields correctly, and those factorizations can be used to 
compute the isomorphisms.
"""

from sage.all import (
    var,
    NumberField,
    NumberFieldElement,
    PolynomialRing,
    localvars,
    factor,
    coerce,
    hom,
    CC,
)
from sage.libs.pari.convert_sage import gen_to_sage


def convert_polmod(polmod, name=None, mod_variable="x"):
    """
    Only implemented for univariate polmods for now. If a variable is passed to the
    name parameter, it needs tohave been initialized as a variable in sage. This
    conforms to the usage of other functions in sage.

    It takes a PARI polmod Mod(f,g) and converts it to the image inside the number
    field defined by g.

    In particular, it doesn't implement any general quotient ring reduction maps and
    will fail if the modulus factors. So it's ill-suited to etale algebra type stuff.

    This could be useful if we need to work with PARI more closely later, but we ended
    up not using it in the functions we really care about.
    Aug-8-2020
    """
    mod_variable = var(mod_variable)
    pari_lift = polmod.lift()
    pari_mod = polmod.mod()
    pari_var = polmod.variable()
    lift_var = pari_var if name == None else name
    sage_mod = gen_to_sage(pari_mod, {str(pari_var): mod_variable})
    sage_lift = gen_to_sage(pari_lift, {str(pari_var): lift_var})
    sage_nf = NumberField(sage_mod, lift_var)
    sage_nf_elt = NumberFieldElement(sage_nf, sage_lift)
    return sage_nf_elt


def isomorphisms_between_number_fields(domain_field, codomain_field):
    """
    Takes in two sage number fields are returns a list of isomorphisms between them.
    The algorithm is simple, but its speed relies on how quickly the factorizations
    can be computed.  There are other algorithms for computing field isomorphisms,
    and they might be added to this module in some form at sometime in the future.
    For fields arising from Kleinian groups (which are often of degree less than 100
    with discriminants that are tractable), this function should be reasonably fast.

    7-Aug-2020
    """
    polynomial_ring_over_codomain_field = PolynomialRing(codomain_field, "x")
    x = polynomial_ring_over_codomain_field.gen()
    domain_min_poly = domain_field.defining_polynomial().change_variable_name("x")
    poly_to_factor = polynomial_ring_over_codomain_field.coerce(domain_min_poly)
    factorization = factor(poly_to_factor)
    iso_list = []
    for factor_with_multiplicity in factorization:
        if factor_with_multiplicity[0].degree() == 1:
            iso_list.append(domain_field.hom([-factor_with_multiplicity[0](0)]))
    return iso_list


def transfer_embedding(isomorphism):
    """
    This function takes an isomorphism whose domain is a number field with a specified
    embedding and codomain is a number field (with or without an embedding). It returns
    a complex number corresponding to an embedding of the generator for the codomain.
    Under this embedding, the image of the generator of the domain field will map
    to the same complex number as it did under the specified embedding of the domain.

    As a side note, the output of this function should be independent of the actual
    choice of isomorphisms. That is, if there are multiple isomorphisms between the
    domain and codomain, the numerical value of the generator (and hence choice of
    embedding) should be the same.
    """
    domain = isomorphism.domain()
    codomain = isomorphism.codomain()
    domain_generator = domain.gen()
    domain_numerical_root = domain.gen_embedding()
    if domain_numerical_root is None:
        raise AttributeError("There is no specified embedding for the number field.")
    codomain_embeddings = [
        embedding
        for embedding in codomain.real_embeddings() + codomain.complex_embeddings()
    ]
    domain_generator_image = isomorphism(domain_generator)
    special_embedding = min(
        codomain_embeddings,
        key=lambda embedding: abs(
            domain_numerical_root - embedding(domain_generator_image)
        ),
    )
    return special_embedding(codomain.gen())