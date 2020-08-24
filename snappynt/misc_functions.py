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
    word = word1+word2+word1_inverse+word2_inverse
    return word