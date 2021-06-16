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
    RR,
    I
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

def canonical_embedding(field_with_embedding):
    """
    It seems sage doesn't have a built-in way to get this map.
    """
    return min(field_with_embedding.complex_embeddings(), key= lambda embedding : abs(CC(field_with_embedding.gen_embedding())-CC(embedding.im_gens()[0])))

def same_subfield_of_CC(field1, field2):
    try:
        iso = isomorphisms_between_number_fields(field1,field2)[0]
    except IndexError:
        return False
    gen1 = field1.gen()
    transfered_gen = iso(gen1)
    automorphisms = field2.automorphisms()
    orbit = [automorphism(transfered_gen) for automorphism in automorphisms]
    embedding2 = canonical_embedding(field2)
    embedded_orbit = [embedding2(elt) for elt in orbit]
    all_im_gens = [embedding.im_gens()[0] for embedding in field1.complex_embeddings()]
    found_elts = [min(all_im_gens, key=lambda x : abs(CC(x)-CC(elt))) for elt in embedded_orbit]
    leftovers = [elt for elt in all_im_gens if elt not in found_elts]
    coerced_elt = min(all_im_gens, key=lambda x : abs(CC(x)-CC(field1.gen_embedding())))
    if coerced_elt in leftovers:
        return False
    elif coerced_elt in found_elts:
        return True
    else:
        raise
    """
    difference_between_embeddings = min(a-b for a in all_im_gens for b in all_im_gens if a != b)
    if min(embedding1(gen1)-elt for elt in embedded_orbit) < difference_between_embeddings/2:
        return True
    else:
        return False
    """

def transfer_embedding(isomorphism):
    """
    This function takes an isomorphism whose domain is a number field with a specified
    embedding and codomain is a number field (with or without an embedding). It returns
    a complex number corresponding to an embedding of the generator for the codomain.
    Under this embedding, the image of the generator of the domain field will map
    to the same complex number as it did under the specified embedding of the domain.

    The basic logic here is to take a generator for the domain with a specified
    embedding into CC. This amounts to some numerical value for this generator. Then we
    compare the image of the generator under the various embeddings of the codomain to
    see which one gets closest to the original numerical value. In terms of the 
    variables this means we compare domain_numerical_root and 
    embedding(domain_generator_image) under the various embeddings of the codomain.
    """
    domain = isomorphism.domain()
    codomain = isomorphism.codomain()
    domain_numerical_root = domain.gen_embedding()
    domain_generator_image = isomorphism(domain.gen())
    if domain_numerical_root is None:
        raise AttributeError("There is no specified embedding for the number field.")
    # Sage's complex_embeddings() gives the real ones as well.
    codomain_embeddings = [embedding for embedding in codomain.complex_embeddings()]
    special_embedding = min(
        codomain_embeddings,
        key=lambda embedding: abs(
            CC(domain_numerical_root) - CC(embedding(domain_generator_image))
        ),
    )
    return CC(special_embedding(codomain.gen()))
"""
def compare_embeddings(field, first_numerical_root, second_numerical_root=None):
    
    Tests whether the two numerical roots define the same embedding. This is to sidestep
    issues of numerical precision. Sage might also have a way to do this, but I couldn't
    find it. Basically the two numerical roots passed should live in some kind of real
    or complex field in sage, but there are many of those (RIFs, ones with specified
    precision, lazy fields, etc.), and the numerical root may come to us in various
    ways.

    One can pass in only a field and one numerical root if the field comes with an
    embedding already attached to it.
    
    generator = field.gen()  # Assumes field is given by a single generator I guess.
    second_numerical_root = (
        field.gen_embedding()
        if second_numerical_root is None
        else second_numerical_root
    )
    if second_numerical_root is None:
        raise AttributeError("Got too few embeddings.")
    embeddings = [embedding for embedding in field.complex_embeddings()]
    first_embedding = numerical_root_to_embedding(field, first_numerical_root)
    second_embedding = numerical_root_to_embedding(field, second_numerical_root)
    second_embedding_conjugate = numerical_root_to_embedding(field, second_numerical_root, conjugate_embedding=True)
    return (first_embedding == second_embedding or first_embedding == second_embedding_conjugate)
"""

"""
def isomorphic_respecting_embeddings(first_field, second_field):
    
    This compares two number fields with distinguished places and checks whether they're
    isomorphic and that there is an isomorphism such that composing the isomorphism with
    the embedding (or its conjugate) of the second field gives the specified embedding
    of the first field. This is equivalent to checking whether the two embedded fields
    are the same subfield of the complex numbers (or conjugates of one another).

    This is a little too implicit. Needs some more documentation eventually.
    
    iso_list = isomorphisms_between_number_fields(first_field, second_field)
    for isomorphism in iso_list:
        transfered_root = transfer_embedding(isomorphism)
        if compare_embeddings(second_field, transfered_root):
            return True
    return False
"""
def run_tests():
    """
    A test bench for the various functions in this module. Probably one day add better
    names for everything and add more tests. The convention is True should always mean
    the test ran correctly.
    """
    # Comparing Embeddings
    """
    x = var("x")
    log_dict = dict()
    Field1 = NumberField(x**2+1, "i", embedding=I)
    Field2 = NumberField(x**2+1, "minusi", embedding=-I)
    # Checks that conjugate embeddings are considered the same.
    log_dict['Conjugate embeddings for QQ(i)'] = isomorphic_respecting_embeddings(Field1, Field2)
    # Checks that giving a nonintegral minimal polynomial doesn't mess anything up.
    # The issue would be finding the isomorphism in the first place.
    Field3 = NumberField(x**2+2*x+(5/4), "a", embedding=-1+(1/2)*I)
    log_dict['Integral and nonintegral minimal polynomials'] = isomorphic_respecting_embeddings(Field1, Field3)
    # These next two fields are the trace fields of the knots 6_1 and 7_7 respectively.
    # We test that we find an isomorphism between them and that they are distinguished
    # by their embeddings into the complex numbers.
    FieldSixOne = NumberField(x**4+x**2-x+1, "b", embedding=CC(0.547423794586058 + 0.585651979689573*I))
    FieldSevenSeven = NumberField(x**4+x**2-x+1, "c", embedding=CC(-0.547423794586059 - 1.12087348993706*I))
    log_dict['Distinguishes the trace fields of 6_1 and 7_7'] = (bool(isomorphisms_between_number_fields(FieldSixOne, FieldSevenSeven)) and not isomorphic_respecting_embeddings(FieldSixOne, FieldSevenSeven))
    return log_dict
    """
    x = var("x")
    log_dict = dict()
    Field1 = NumberField(x**2+1, "i", embedding=I)
    Field2 = NumberField(x**2+1, "minusi", embedding=-I)
    # Checks that conjugate embeddings are considered the same.
    log_dict['Conjugate embeddings for QQ(i)'] = same_subfield_of_CC(Field1, Field2)
    # Checks that giving a nonintegral minimal polynomial doesn't mess anything up.
    # The issue would be finding the isomorphism in the first place.
    Field3 = NumberField(x**2+2*x+(5/4), "a", embedding=-1+(1/2)*I)
    log_dict['Integral and nonintegral minimal polynomials'] = same_subfield_of_CC(Field1, Field3)
    # These next two fields are the trace fields of the knots 6_1 and 7_7 respectively.
    # We test that we find an isomorphism between them and that they are distinguished
    # by their embeddings into the complex numbers.
    FieldSixOne = NumberField(x**4+x**2-x+1, "b", embedding=CC(0.547423794586058 + 0.585651979689573*I))
    FieldSevenSeven = NumberField(x**4+x**2-x+1, "c", embedding=CC(-0.547423794586059 - 1.12087348993706*I))
    log_dict['Distinguishes the trace fields of 6_1 and 7_7'] = (bool(isomorphisms_between_number_fields(FieldSixOne, FieldSevenSeven)) and not same_subfield_of_CC(FieldSixOne, FieldSevenSeven))
    return log_dict