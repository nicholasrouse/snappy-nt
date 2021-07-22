"""
Module for odds and ends functions that thematically didn't fit into another module.
"""


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


def ramified_real_places(quaternion_algebra):
    """
    Takes in a quaternion algebra over a number field and returns a list of ramified
    places as maps. Obviously this could basically be a somewhat long one-liner, but
    I think it's a bit nicer this way.

    To do: remove this as it is superceded in the QuaternionAlgebraNF module.
    """
    a, b = quaternion_algebra.invariants()
    field = quaternion_algebra.base_ring()
    # Could use real_embeddings() as well.
    real_places = field.real_places()
    ramified_places = [place for place in real_places if place(a) < 0 and place(b) < 0]
    return ramified_places


def aan_iterator(list_of_aans):
    """
    There's an issue with iterating over a ListOfApproximateAlgebraicNumbers.
    """
    for i, item in enumerate(list_of_aans):
        yield item
        if i >= len(list_of_aans) - 1:
            break
