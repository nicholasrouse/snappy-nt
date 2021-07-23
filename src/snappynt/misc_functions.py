"""
Module for odds and ends functions that thematically didn't fit into another module.
"""
from sage.all import CC, ZZ, ComplexField, NumberField, PolynomialRing
from snappy.snap.find_field import (
    ApproximateAlgebraicNumber,
    ExactAlgebraicNumber,
    ListOfApproximateAlgebraicNumbers,
)


def commutator_of_words(word1, word2):
    """
    This takes two strings that represent words in a quotient of a free group and
    forms their commutator as another word. This is probably a bit of band-aid
    function since we probably shouldn't be passing around group elements as
    strings and instead actually make them elements of some class. Anyway, this one
    will work as long as the words are words in a single letter and the other case
    represents the inverse. I.e. the group is a quotient of the free group on the 26
    letters a,b,...,z with inverses A,B,...,Z respectively.
    """
    word1_inverse = word1[::-1].swapcase()
    word2_inverse = word2[::-1].swapcase()
    word = word1 + word2 + word1_inverse + word2_inverse
    return word


def aan_iterator(list_of_aans):
    """
    There's an issue with iterating over a ListOfApproximateAlgebraicNumbers.
    Can use their list() method as well.
    """
    for i, item in enumerate(list_of_aans):
        yield item
        if i >= len(list_of_aans) - 1:
            break


def make_aan(poly, root):
    name = str(poly.variables()[0])
    field = NumberField(poly, name, embedding=root)
    elt = field.gen_embedding()

    def defining_func(prec):
        return ComplexField(prec)(elt)

    aan = ApproximateAlgebraicNumber(defining_func)
    aan._min_poly = PolynomialRing(ZZ, name)(poly)
    return aan


def conjugate_field(field):
    """
    field should be a NumberField with an embedding.
    """
    poly = field.defining_polynomial()
    root = field.gen_embedding()
    name = str(field.gen())
    return NumberField(poly, name, CC(root.conjugate()))


def make_aan_conjugate(aan):
    if isinstance(aan, ExactAlgebraicNumber):
        poly = aan._min_poly
        root = aan._approx_root
        new_aan = ExactAlgebraicNumber(poly, root.conjugate())
        return new_aan
    if isinstance(aan, ApproximateAlgebraicNumber):

        def new_defining_function(prec):
            return aan.f(prec).conjugate()

        new_aan = ApproximateAlgebraicNumber(new_defining_function)
        new_aan._min_poly = aan._min_poly  # might be None
        if hasattr(aan, "_default_precision"):
            new_aan._default_precision = aan._default_precision
        if hasattr(aan, "_approx_root"):
            new_aan._approx_root = aan._approx_root.conjugate()
        return new_aan
    if isinstance(aan, ListOfApproximateAlgebraicNumbers):

        def new_defining_function(prec):
            return [elt.conjugate() for elt in aan.f(prec)]

        new_aan = ListOfApproximateAlgebraicNumbers(new_defining_function)
        for key in aan._field:
            if aan._field[key] is not None:
                old_field = aan._field[key][0]
                new_field = conjugate_field(old_field)
                # new_aan_elt should be an ExactAlgebraicNumber.
                new_aan_elt = make_aan_conjugate(aan._field[key][1])
                # The expressions (the last entry) should be the same, but the parent
                # ring will be different.
                exact_expressions = [new_field(str(elt)) for elt in aan._field[key][2]]
                new_aan._field[key] = (new_field, new_aan_elt, exact_expressions)
        return new_aan
