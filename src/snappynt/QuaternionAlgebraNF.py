"""
A module for quaternion algebras over number fields. The class will inherit from Sage's
quaternion algebras. The key method will be an is_isomorphic method which leverages the
Albert-Brauer-Hasse-Noether theorem to check for isomorphism. The primary application
will be for quaternion algebras associted to Kleinian groups, so there might be a few
other methods which are attuned to this purpose.

Our first pass will mostly be focused on
1. Computing ramification data such as real and finite places and residue
characteristics.
2. Providing a convenient __str__ representation for printing in other modules.
3. A fairly general method for testing for isomorphism between two quaternion algebras.

For now, the class QuaternionAlgebraNF won't incorporate the full and rich structure of
quaternion algebras over number fields. E.g. we don't provide methods for base changing
or finding splitting fields. In principle this information can be read off from the
ramification data, and maybe one day we'll add it. To that end one should take caution
in using certain sage methods that this class inherits. For example, sage provides a way
to base extend, but doing so then asking for ramification data will currently produce an
incorrect answer because base changing doesn't remove previously computed ramification
information.

This class basically doesn't work when the base ring is QQ amusingly. The issue is many
number field methods in sage are not implemented for QQ. Maybe there's a clean fix for
this, but for now I'm going to leave it since the functionality already pretty much
exists in Sage for quaternion algebras over the rationals.
"""
from collections import Counter

from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra_ab as _QA_ab
from sage.all import QQ, radical


class QuaternionAlgebraNF(_QA_ab):
    def __init__(
        self,
        base_ring,
        a,
        b,
        names="i,j,k",
    ):
        """
        On initialization we by default compute the ramification set. We don't actually
        type-test that the base_ring is a number field, but we do print a warning if it
        doesn't look like one of Sage's NumberField.
        """
        if base_ring == QQ:
            raise NotImplementedError(
                "Use QuaternionAlgebra when the base field is the rational numbers."
            )
        a, b = base_ring(a), base_ring(b)  # Do the coercion explictly.
        self._ramified_dyadic_residue_chars = Counter()
        self._ramified_real_places = set()
        self._ramified_real_places_known = False
        self._ramified_nondyadic_places = set()
        self._ramified_nondyadic_residue_chars = Counter()
        self._ramified_nondyadic_places_known = False
        self._ramified_dyadic_places = set()
        self._ramified_dyadic_residue_chars = Counter()
        self._ramified_dyadic_places_known = False
        self._ramified_dyadic_places_dict = dict()
        _QA_ab.__init__(self, base_ring, a, b, names)

    def is_ramified_at(self, place):
        """
        Returns True or False depending on whether the algebra is ramified at the place.
        The place parameter can be either a prime ideal of the base_ring, or a real
        place of the base ring as contained in the output of base_ring.real_places().
        """
        if self.base_ring().hilbert_symbol(*self.invariants(), place) == -1:
            try:
                if place.absolute_norm() % 2 != 0:
                    self._ramified_nondyadic_places.add(place)
                else:
                    self._ramified_dyadic_places_dict[place] = True
                    self._ramified_dyadic_places.add(place)
            except AttributeError:
                self._ramified_real_places.add(place)
            return True
        else:
            if place.absolute_norm() % 2 == 0:
                self._ramified_dyadic_places_dict[place] = False
            return False

    def ramified_real_places(self):
        """
        Takes in a quaternion algebra over a number field and returns a list of ramified
        places as maps. Obviously this could basically be a somewhat long one-liner, but
        I think it's a bit nicer this way. The output by the way is a set as the places
        are hashable in Sage.
        """
        if self._ramified_real_places_known:
            return self._ramified_real_places
        field = self.base_ring()
        real_places = set(field.real_places()) - self._ramified_real_places
        ramified_places = set(
            [
                place
                for place in real_places
                if field.hilbert_symbol(*self.invariants(), place) == -1
            ]
        )
        self._ramified_real_places |= ramified_places
        self._ramified_real_places_known = True
        return self._ramified_real_places

    def ramified_nondyadic_places(self):
        if self._ramified_nondyadic_places_known:
            return self._ramified_nondyadic_places
        a, b = self.invariants()
        primes_dividing_a = {
            prime
            for (prime, multiplicity) in self.base_ring().ideal(a).factor()
            if multiplicity % 2 != 0
            and prime.absolute_norm() % 2 != 0
            and prime not in self._ramified_nondyadic_places
        }
        primes_dividing_b = {
            prime
            for (prime, multiplicity) in self.base_ring().ideal(b).factor()
            if multiplicity % 2 != 0
            and prime.absolute_norm() % 2 != 0
            and prime not in self._ramified_nondyadic_places
        }
        ramified_primes = {
            prime
            for prime in primes_dividing_a | primes_dividing_b
            if self.base_ring().hilbert_symbol(a, b, prime) == -1
        }
        self._ramified_nondyadic_places |= ramified_primes
        self._ramified_nondyadic_places_known = True
        return self._ramified_nondyadic_places

    def ramified_dyadic_places(self):
        if self._ramified_dyadic_places_known:
            return self._ramified_dyadic_places
        dyadic_primes = sorted(
            (
                prime
                for prime in self.base_ring().ideal(2).prime_factors()
                if prime not in self._ramified_dyadic_places_dict
            ),
            key=lambda prime: prime.residue_class_degree(),
        )
        if self._ramified_nondyadic_places_known and self._ramified_real_places_known:
            dyadic_primes, last_prime = dyadic_primes[:-1], dyadic_primes[-1]
        else:
            last_prime = None
        for prime in dyadic_primes:
            if self.base_ring().hilbert_symbol(*self.invariants(), prime) == -1:
                self._ramified_dyadic_places_dict[prime] = True
                self._ramified_dyadic_places.add(prime)
            else:
                self._ramified_dyadic_places_dict[prime] = False
        if last_prime is not None:
            if (
                len(
                    self._ramified_real_places
                    | self._ramified_dyadic_places
                    | self._ramified_nondyadic_places
                )
                % 2
                == 1
            ):
                self._ramified_dyadic_places_dict[last_prime] = True
                self._ramified_dyadic_places.add(last_prime)
            else:
                self._ramified_dyadic_places_dict[last_prime] = False
        self._ramified_dyadic_places_known = True
        return self._ramified_dyadic_places

    def ramified_finite_places(self):
        """
        Computes the ramified finite places as a set. We note that the
        parent class has a method ramified_primes() that returns a list. We think that a
        set is better suited and this method also gives us some flexibility to compute
        things like the residue characteristics along the way.

        We try to be a little bit clever with Hilbert reciprocity. In practice the
        dyadic places can take a long time to compute, so we try to avoid finding one of
        them explictly. The point is the total number of ramified real and finite places
        is even.

        We check that the multiplicity of the primes that a and b belong to are odd, but
        this might be unnecessary if sage is smart enough to quickly compute their local
        symbols. On the other hand it probably cuts down slightly on the total number of
        function calls which should be a win. On the third hand having primes appear to
        powers higher than 1 might be so rare in practice that this isn't worth it.
        """
        return self.ramified_nondyadic_places() | self.ramified_dyadic_places()

    def ramified_places(self):
        return self.ramified_real_places() | self.ramified_finite_places()

    def ramified_nondyadic_residue_characteristics(self):
        self._ramified_nondyadic_residue_chars = Counter(
            [
                radical(place.absolute_norm())
                for place in self.ramified_nondyadic_places()
            ]
        )
        return self._ramified_nondyadic_residue_chars

    def ramified_dyadic_residue_characteristics(self):
        self._ramified_dyadic_residue_chars = Counter(
            [radical(place.absolute_norm()) for place in self.ramified_dyadic_places()]
        )
        return self._ramified_dyadic_residue_chars

    def ramified_residue_characteristics(self):
        """
        Find the residue characteristics of the ramified places. It will attempt to
        compute the ramified places if they're not known. The residue characteristics
        are a Counter (morally a multiset) to keep track of multiplicity.
        """
        return (
            self.ramified_dyadic_residue_characteristics()
            | self.ramified_nondyadic_residue_characteristics()
        )

    def is_division_algebra(self):
        """
        Overrides that from Sage's QuaternionAlgebra_ab class. Just whether there is any
        ramification.
        """
        return bool(self.ramified_finite_places()) or bool(self.ramified_real_places())

    def is_matrix_ring(self):
        """
        Overrides the method from the parent class. This (obviously) exploits the
        dichotomy of matrix and division algebra for quaternion algebras.
        """
        return not self.is_division_algebra()

    def new_QA_via_field_isomorphism(self, isomorphism):
        """
        This will use the isomorphism and self to produce a new quaternion algebra over
        isomorphism's codomain.

        That is, if isomorphism is f: K->L and self is (a,b)_K, then this method returns
        (f(a), f(b))_L.

        The domain of isomorphism must be self.base_ring() and its codomain must be a
        NumberField. In particular, this is not appropriate for changing to local
        fields. This method will not recompute the ramification, but will transfer it
        via the isomorphism, remembering whether it has been computed.
        """
        if self.base_ring() != isomorphism.domain():
            raise ValueError("The isomorphism domain must be the base ring of self.")
        a, b = [isomorphism(gen) for gen in self.invariants()]
        new_QA = QuaternionAlgebraNF(isomorphism.codomain(), a, b)
        new_QA._ramified_real_places_known = self._ramified_real_places_known
        new_QA._ramified_nondyadic_places_known = self._ramified_nondyadic_places_known
        new_QA._ramified_dyadic_places_known = self._ramified_dyadic_places_known
        new_QA._ramified_nondyadic_residue_chars = (
            self._ramified_nondyadic_residue_chars
        )
        new_QA._ramified_dyadic_residue_chars = self._ramified_nondyadic_residue_chars
        new_QA._ramified_real_places = set(
            isomorphism(place) for place in self._ramified_real_places
        )
        new_QA._ramified_nondyadic_places = set(
            isomorphism(place) for place in self._ramified_nondyadic_places
        )
        new_QA._ramified_dyadic_places = set(
            isomorphism(place) for place in self._ramified_dyadic_places
        )
        new_QA._ramified_dyadic_places_dict = {
            isomorphism(place): self._ramified_dyadic_places_dict[place]
            for place in self._ramified_dyadic_places_dict
        }
        return new_QA

    def is_isomorphic(self, other):
        """
        Takes in two QuaternionAlgebraNFs. If their base fields do not compare as equal,
        using == in sage, then a ValueError is raised. See the method
        same_ramification_via_isomorphism for ways to compare algebras defined over
        different but isomorphic number fields. The method uses the
        Albert-Brauer-Hasse-Noether theorem for quaternion algebras. The net effect is
        to check whether the ramification sets are the same. The method tries to
        distinguish the algebras via cheaper invariants before more expensive ones.

        Note that this method does not produce an isomorphism between the algebras.
        """
        if self == other:
            return True
        self_field = self.base_ring()
        other_field = other.base_ring()
        if self_field != other_field:
            raise ValueError(
                "The quaternion algebras do not appear to be defined over a common number field. They must compare as equal (==) in sage."
            )
        if self.ramified_real_places() != other.ramified_real_places():
            return False
        if self.ramified_nondyadic_places() != other.ramified_nondyadic_places():
            return False
        if self.ramified_dyadic_places() != other.ramified_dyadic_places():
            return False
        return True

    def same_ramification_via_isomorphism(self, other, isomorphism):
        """
        Given two QuaternionAlgebraNFs and an isomorphism between their base fields,
        this function returns whether the two quaternion algebras are isomorphic once
        they're thought of over a common field. That is, if A is a K-algebra and B is an
        L-algebra, then an isomorphism K->L lets us think of A as an L-algebra and makes
        isomorphism a well-posed question.

        The isomorphism should go from self.base_ring() to other.base_ring().

        The simplest approach is to just create a new quaternion algebra from the
        isomorphism. This will work, but it will cause all the ramification to be
        recomputed, including the potentially expensive dyadic places, so instead we
        use the field isomorphism to move the ramification. We also check to see whether
        there are any obvious obstructions to becoming isomorphic such as different
        numbers of ramified places above some rational or infinite place.
        """
        self_field = self.base_ring()
        other_field = other.base_ring()
        if isomorphism.domain() != self_field or isomorphism.codomain() != other_field:
            raise ValueError(
                f"""The isomorphism does not appear to have correct domain and codomain.
                Domain: expected {self.base_ring()} got {isomorphism.domain()}.
                Codomain: expected {other.base_ring()} got {isomorphism.codomain()}."""
            )
        new_QA = self.new_QA_via_field_isomorphism(isomorphism)
        return new_QA.is_isomorphic(other)

    def ramification_string(
        self,
        full_finite_ramification=True,
        full_real_ramification=True,
        show_hilbert_symbol=True,
        show_field_data=False,
        leading_char="",
    ):
        """
        This returns a somewhat large string that contains the ramification data in a
        human readable format. Its intended use is to display in console sessions or to
        write to some output file. It doesn't read off the field data by default. The
        reason is that for Kleinian groups, the fields are important invariants in their
        own right, so are reproduced elsewhere. Passing in full_ramification=False will
        only give the residue characteristics rather than the ideals. Similarly the
        full_real_ramification keyword argument determines whether the numerical values
        of the the generator for the field at the ramified real places are specified. If
        not, then only the number of ramified real places will be given. Finally the
        leading_char keyword argument allows one to prepend each line of information
        with a particular string. The use case I have in mind is a tab to display with
        other information.

        One should generally not try to access information by parsing the output of this
        method: all the data it contains is much more easily accessible by other methods
        in this class.
        """
        data_strings = list()
        if show_field_data:
            field_data = "Base field: " + str(self.base_ring())
            data_strings.append(field_data)
        if show_hilbert_symbol:
            hilbert_symbol = "Hilbert symbol: " + str(self.invariants())
            data_strings.append(hilbert_symbol)
        if full_finite_ramification:
            finite_ramification_data = "Finite ramification: " + str(
                self.ramified_finite_places()
            )
            data_strings.append(finite_ramification_data)
        ramified_residue_chars_data = (
            "Finite ramification residue characteristics: "
            + str(self.ramified_residue_characteristics())
        )
        data_strings.append(ramified_residue_chars_data)
        number_of_ramified_real_places = "Number of ramified real places: " + str(
            len(self.ramified_real_places())
        )
        data_strings.append(number_of_ramified_real_places)
        if full_real_ramification:
            var = str(self.base_ring().gen())
            small_maps = list()
            for embedding in self.ramified_real_places():
                numerical_value = str(embedding.im_gens()[0])
                small_maps = [var + " |--> " + numerical_value]
            ramified_real_places_data = "Ramified real places: " + str(small_maps)
            data_strings.append(ramified_real_places_data)
        output_string = str()
        for element in data_strings:
            output_string += leading_char + element + "\n"
        return output_string.rstrip()
