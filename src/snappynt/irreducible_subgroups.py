from sage.all import QQ, ZZ


def within_epsilon(a, b, epsilon_coefficient=1):
    """
    I feel like Sage should probably have a way to do this anyway, but I was unable to
    find it. What this does is take in two number from fixed precision fields and test
    whether they're within the smallest nonzero positive number representatble in both
    fields times the epsilon_coefficient. One can tweak the epsilon_coefficient to be
    more or less conservative (depending on the application) with equality testing.
    Note that setting the epsilon_coefficient to less than 1 might introduce some
    weirdness since you'll have a quantity smaller than the theoretically smallest
    representable quantity in the fixed precision fields. The idea, by the way, is that
    by increasing the epsilon_coefficient, one can be more sure that two numbers are in
    fact different. Of course the very notion of epsilon-nearness relies on some
    underlying metric, so if the elements lack an absolute value function, the function
    will fail.

    Also, this is used extensively in the irreducible_subgroups module, but is perhaps
    thematically better suited to living in another module.

    Comment last updated: Aug-15 2020
    """
    if isinstance(a, int):
        base_field_a = QQ
    else:
        base_field_a = a.parent().base()
    if isinstance(b, int):
        base_field_b = QQ
    else:
        base_field_b = b.parent().base()
    if base_field_a.is_exact() and base_field_b.is_exact():
        epsilon = 0
    else:
        epsilon = max(base_field_a.epsilon(), base_field_b.epsilon())
    if abs(a - b) <= epsilon_coefficient * epsilon or (epsilon == 0 and a == b):
        return True
    else:
        return False


def is_parabolic(element, epsilon_coefficient=10):
    """
    Tests whether an element has trace 2 using the within_epsilon method.  Actually it
    tests whether the element is within epsilon of 2 (see within_epsilon and
    generate_reducible_subgroup functions for more details). By default the function
    is conservative about whether certifying that an element is parabolic, so the
    default use should be to verify that an element is NOT parabolic.

    My plan right now is to use the cusp info from SnapPy to determine whether a Kleinian
    group contains parabolics. Realistically, this is probably pretty safe, but there's no
    way that I know of to exclude the possibility that the trace differs from 2 by an amout
    smaller than what can be detected by for a given fixed precision.
    """
    if within_epsilon(element.trace(), ZZ(2), epsilon_coefficient):
        return True
    else:
        return False


def generate_reducible_subgroup(g, h, epsilon_coefficient=10):
    """
    Tests whether the subgroup generated by g and h is irreducible. The function
    actually tests whether the trace of the commutator of g and h is less than
    epsilon_coeffient times epsilon from 2 where epsilon is the smallest nonzero quantity
    that the base field with specified precision can recognize. This implies that the
    function might incorrectly say that two elements generate a reducible subgroup when
    their commutator has trace closer (but not equal) to 2 that the ambient precision.
    On the other hand, it will never certify that a subgroup is irreducible when it is
    not.

    I'm not actually sure if there's a slicker way to test for equality or exactly
    the implementation details of equality testing in Sage's various rings.

    Also, this function as currently written is a bittle amusingly short, but I'm going
    to leave it because it will be called elsewhere, and it lets me change how
    irreducibility is certified in the future.

    Comment last updated: Aug-13 2020
    """
    if within_epsilon(g.commutator(h).trace(), ZZ(2), epsilon_coefficient):
        return True
    else:
        return False


def enumerate_words(rank, power=1):
    """
    This just enumerates words from a free group in ``Tietze notation". E.g. for two
    generators a and b, [1,2,1] is the element aba; [2,-1,-2] is ba^(-1)b^(-1). The
    rank variable is the rank of the free group, i.e. the number of generators.
    This can probably be rewritten recursively to avoid writing an explicit
    algorithm for computing the next element from the previous. The upshot is that this
    should play nicely with most groups given as quotients of free groups in sage.

    The power argument allows for enumeration of powers of words. E.g. if rank=2 and
    power=2, then the generator will yield [1,1], [2,2], [-1,-1], [-2,-2]. This option
    is used to find group elements used in the construction of invariant trace fields
    and quaternion algebras.
    """

    def next_integer(integer):
        if integer > 0 and integer < rank:
            return integer + 1
        elif integer == rank:
            return -1
        elif integer < 0:
            return integer - 1

    def next_element(element):
        index_list = [(integer == -rank) for integer in element][::-1]
        index_to_change = len(element) - index_list.index(False) - 1
        new_element = element[:]
        new_element[index_to_change] = next_integer(new_element[index_to_change])
        if index_to_change != len(new_element):
            new_element[index_to_change + 1 :] = (
                len(new_element) - index_to_change - 1
            ) * [1]
        return new_element

    def has_simplification(element):
        for i in range(len(element) - 1):
            if element[i] == -element[i + 1]:
                return True
        return False

    previous_element = [1]
    yield previous_element * power
    while True:
        if previous_element == [-rank] * len(previous_element):
            previous_element = [1] * (len(previous_element) + 1)
            yield previous_element * power
        else:
            previous_element = next_element(previous_element)
            if not has_simplification(previous_element):
                yield previous_element * power


def enumerate_group_elements(group, as_word=False, power=1, verbosity=False):
    """
    This is not especially robust at the moment. It basically only works with groups
    coming from polished_holonomy in SnapPy. In particul.ar it makes use of a __call__
    interface with the group, i.e. to get the element a*b^2*a^(-1) of the group G by
    calling G('abbA'). Unfortunately this pattern won't work with most groups in sage
    or (nonpolished) holonomy groups from SnapPy. Later we can go back and make the
    interface more robust if needed. The use case of this function is to get group
    elements that may be used as generators for various fields or quaternion algebras.

    If as_word=True, then the result is returned as a word in the generators.
    Last updated: Aug-26 2020
    """
    gens = (
        group.generators()
    )  # These need to all be lower case letters right now, sadly.
    inverses = [gen.upper() for gen in gens]
    rank = len(gens)
    # The tietze_alphabet below just looks like [1,2,...,rank,-1,-2,...-rank]
    # tietze_alphabet = list(range(1, rank + 1)) + [-i for i in range(1, rank + 1)]

    def tietze_conversion(integer):
        if integer > 0:
            return gens[integer - 1]
        elif integer < 0:
            return inverses[-integer - 1]

    tietze_words_generator = enumerate_words(rank, power=power)
    while True:
        tietze_word = next(tietze_words_generator)
        list_word = [tietze_conversion(integer) for integer in tietze_word]
        word = "".join(list_word)
        if verbosity:
            print(word)
        if as_word:
            yield word
        else:
            yield group(word)


def find_hilbert_symbol_words(group, power=1, epsilon_coefficient=100):
    """
    Given a group, enumerates pairs of words until two, g and h, are found such that
    g and [g,h] (the commutator) are not parabolic. The function is so named because
    these are the conditions for tr^2(g)-4 and tr[g,h]-2 to be entries of the
    quaternion algebra spanned by elements of the group over the field generated by
    traces of the group. Again, the interface needs to more or less conform to that of
    SnapPy's polished holonomy groups. E.g. G('abA') needs to make sense.

    Caution: The returned values are a pair of words (g,h), not the actual Hilbert
    symbol or approximations thereof.

    We had to power up the epsilon_coefficients because we were missing some parabolics.
    I should probably rethink the implementation a bit. Maybe a try--except block where
    we can catch when we convert an approximate algebraic integer to 0 or 2 via LLL.
    """
    gen = enumerate_group_elements(group, as_word=True, power=power)
    while True:
        word1 = next(gen)
        word2 = next(gen)
        numerical_elt1 = group(word1)
        numerical_elt2 = group(word2)
        if not is_parabolic(
            numerical_elt1, epsilon_coefficient=epsilon_coefficient
        ) and not generate_reducible_subgroup(
            numerical_elt1, numerical_elt2, epsilon_coefficient=epsilon_coefficient
        ):
            return (word1, word2)
        if not is_parabolic(
            numerical_elt2, epsilon_coefficient=epsilon_coefficient
        ) and not generate_reducible_subgroup(
            numerical_elt1, numerical_elt2, epsilon_coefficient=epsilon_coefficient
        ):
            return (word2, word1)
