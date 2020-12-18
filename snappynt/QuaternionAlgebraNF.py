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
from sage.all import QuaternionAlgebra, QQ, radical, pari
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra_ab
from sage.rings.number_field.number_field import is_NumberField
import field_isomorphisms
from collections import Counter

def convert_QA_toQANF(quaternion_algebra, delay_computations=False, suppress_warnings=False):
    """
    This converts a Sage quaternion algebra over a number field to an instance of the
    subclass QuaternionAlgebraNF. It's possible we can support passing in one of Sage's
    QuaternionAlgebra objects to __init__ for QuaternionAlgebraNF, but it's going
    to take some contortions, so I think this is better.

    We don't actually type-test that the base_ring is a NumberField just like we don't
    in QuaternionAlgebraNF. There is a warning printed though coming from the __init__
    method unless suppress_warnings=True.
    """
    field = quaternion_algebra.base_ring()
    a, b = quaternion_algebra.invariants()
    names = [str(gen) for gen in quaternion_algebra.gens()]
    return QuaternionAlgebraNF(field, a, b, names=names, delay_computations=delay_computations, suppress_warnings=suppress_warnings)

def pari_local_symbol(a,b,prime):
    """
    This takes sage types for number field elements a,b and a prime ideal of the same
    number field. However the calculation is outsourced to PARI.

    This might not be totally necessary. It seems to be more or less how Sage does it
    anyway.
    """
    field = pari(a.parent())
    a = pari(a)
    b = pari(b)
    prime = prime.pari_prime()
    return pari.nfhilbert(field, a, b, prime)


class QuaternionAlgebraNF(QuaternionAlgebra_ab):

    def __init__(self, base_ring, a, b, names="i,j,k", compute_ramification=True, suppress_warnings=False):
        """
        On initialization we by default compute the ramification set. We don't actually
        type-test that the base_ring is a number field, but we do print a warning if it
        doesn't look like one of Sage's NumberField.
        """
        if suppress_warnings and not is_NumberField(base_ring):
            print('The base ring does not appear to be a number field. Proceed with caution.')
        self._ramified_residue_characteristics = None
        self._ramified_real_places = None
        self._ramified_finite_places = None
        self._ramified_dyadic_places_computed = False
        QuaternionAlgebra_ab.__init__(self, base_ring, a, b, names)
        if compute_ramification:
            self.ramified_real_places()
            self.ramified_finite_places()

    def ramified_real_places(self, force_compute=False):
        """
        Takes in a quaternion algebra over a number field and returns a list of ramified
        places as maps. Obviously this could basically be a somewhat long one-liner, but
        I think it's a bit nicer this way. The output by the way is a set as the places
        are hashable in Sage.
        """
        if self._ramified_real_places and not force_compute:
            return self._ramified_real_places
        else:
            a,b = self.invariants()
            field = self.base_ring()
            real_places = field.real_places()
            ramified_places = set([place for place in real_places if place(a) < 0 and place(b) < 0])
            self._ramified_real_places = ramified_places
            return ramified_places
    
    def ramified_residue_characteristics(self, force_compute=False, ignore_dyadic=False):
        """
        Find the residue characteristics of the ramified places. It will attempt to
        compute the ramified places if they're not known. The residue characteristics
        are a Counter (morally a multiset) to keep track of multiplicity. The
        force_compute option will be passed forward to ramified_finite_places if the
        latter method needs to be invoked.
        """
        if self._ramified_finite_places is None:
            if ignore_dyadic:
                self.ramified_nondyadic_places(force_compute=force_compute)
            else:
                self.ramified_finite_places(force_compute=force_compute)
        if not self._ramified_residue_characteristics or force_compute:
            self._ramified_residue_characteristics = Counter([radical(place.absolute_norm()) for place in self._ramified_finite_places])
        return self._ramified_residue_characteristics
    
    def ramified_nondyadic_places(self, force_compute=False):
        if not self._ramified_finite_places or force_compute:
            a, b = self.invariants()
            ramified_finite_primes = None
            primes_dividing_a = {prime for (prime, multiplicity) in self.base_ring().ideal(a).factor() if multiplicity%2 != 0 and prime.absolute_norm()%2 != 0}
            primes_dividing_b = {prime for (prime, multiplicity) in self.base_ring().ideal(b).factor() if multiplicity%2 != 0 and prime.absolute_norm()%2 != 0}
            ramified_finite_primes = {
                prime for prime in primes_dividing_a | primes_dividing_b if self.base_ring().hilbert_symbol(a, b, prime) == -1
            }
            if self._ramified_finite_places is not None:
                self._ramified_finite_places = self._ramified_finite_places | ramified_finite_primes
            else:
                self._ramified_finite_places = ramified_finite_primes
        answer = {prime for prime in self._ramified_finite_places if prime.absolute_norm()%2 != 0}
        return answer
        
    
    def ramified_finite_places(self, force_compute=False):
        """
        Computes the ramified finite places as a set. Passing in force_compute=True will
        recompute the finite places and the residue characteristics. We note that the
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
        if not self._ramified_finite_places or self._ramified_dyadic_places_computed or force_compute:
            a, b = self.invariants()
            ramified_finite_primes = None
            primes_dividing_a = {prime for (prime, multiplicity) in self.base_ring().ideal(a).factor() if multiplicity%2 != 0 and prime.absolute_norm()%2 != 0}
            primes_dividing_b = {prime for (prime, multiplicity) in self.base_ring().ideal(b).factor() if multiplicity%2 != 0 and prime.absolute_norm()%2 != 0}
            ramified_finite_primes = {
                prime for prime in primes_dividing_a | primes_dividing_b if self.base_ring().hilbert_symbol(a, b, prime) == -1
            }
            total_number_ramified_places = len(self.ramified_real_places(force_compute=force_compute) | ramified_finite_primes)
            # Experimentally we most want to avoid computing dyadic places with large residue class degree.
            dyadic_primes = sorted(self.base_ring().ideal(2).prime_factors(), key=lambda prime : prime.residue_class_degree())
            for prime in dyadic_primes:
                if prime == dyadic_primes[-1]:
                    if total_number_ramified_places%2 != 0:
                        ramified_finite_primes.add(prime)
                else:
                    if self.base_ring().hilbert_symbol(a,b,prime) == -1:
                        ramified_finite_primes.add(prime)
                        total_number_ramified_places += 1
            self._ramified_finite_places = ramified_finite_primes
        self.ramified_residue_characteristics(force_compute=force_compute)
        return self._ramified_finite_places
    
    def is_division_algebra(self):
        """
        Overrides that from Sage's QuaternionAlgebra_ab class. Just whether there is any
        ramification.
        """
        return bool(self._ramified_finite_places) or bool(self._ramified_finite_places)
    
    def is_matrix_ring(self):
        """
        Overrides the method from the parent class. This (obviously) exploits the
        dichotomy of matrix and division algebra for quaternion algebras.
        """
        return not self.is_division_algebra()

    def is_isomorphic(self, other, field_isomorphism=None):
        """
        Given two quaternion algebras over number fields, this function tests for
        isomorphism. The first check is that their base fields are (abstractly) isomorphic. Assuming
        that they are, the function uses the Albert-Brauer-Hasse-Noether theorem for
        quaternion algebras. The net effect is to check whether the ramification sets
        are the same. This can be a bit delicate in total generality, e.g. if self and
        other have isomorphic base fields but they are given by different minimal
        polynomials.

        Fortunately for us, given an isomorphism of fields f: K -> L and an ideal I of
        the ring of integers of K, Sage will correctly interpret f(I) as an ideal of 
        (the ring of integers) of L.

        This might not be especially efficient, but our first pass algorithm is just
        going to be to make a QuaternionAlgebraNF out of other which has the same base
        field as self. We do first try to distinguish by some obvious invariants.

        The expensive portion of this method is actually most likely computing an
        isomorphism of the base number fields. For this reason we allow for an
        isomorphism to be passed in.
        """
        self_field = self.base_ring()
        other_field = other.base_ring()
        if field_isomorphism is None:
            try:
                field_isomorphism = field_isomorphisms.isomorphisms_between_number_fields(self_field, other_field)[0]
            except IndexError:
                return False
        # Trying if other already has some ramification data computed.
        try:
            same_number_of_real_ramification = (len(self._ramified_real_places) == len(other._ramified_real_places))
            same_residue_characteristics = (self._ramified_residue_characteristics == other._ramified_residue_characteristics)
            if not (same_number_of_real_ramification and same_residue_characteristics):
                return False
        except AttributeError: pass
        a, b = [field_isomorphism(gen) for gen in other.invariants()]
        new_quaternion_algebra = QuaternionAlgebraNF(self_field, a, b)
        same_real_ramification = (self._ramified_real_places == new_quaternion_algebra._ramified_real_places)
        same_finite_ramification = (self._ramified_finite_places == new_quaternion_algebra._ramified_finite_places)
        return (same_real_ramification and same_finite_ramification)

    def ramification_string(self, full_finite_ramification=True, full_real_ramification=True, show_hilbert_symbol=True, show_field_data=False, leading_char=''):
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
            field_data = 'Base field: ' + str(self.base_ring())
            data_strings.append(field_data)
        hilbert_symbol = 'Hilbert symbol: ' + str(self.invariants())
        data_strings.append(hilbert_symbol)
        if full_finite_ramification:
            finite_ramification_data = 'Finite ramification: ' + str(self.ramified_finite_places())
            data_strings.append(finite_ramification_data)
        ramified_residue_chars_data = 'Finite ramification residue characteristics: ' + str(self.ramified_residue_characteristics())
        data_strings.append(ramified_residue_chars_data)
        number_of_ramified_real_places = 'Number of ramified real places: ' + str(len(self.ramified_real_places()))
        data_strings.append(number_of_ramified_real_places)
        if full_real_ramification:
            var = str(self.base_ring().gen())
            small_maps = list()
            for embedding in self.ramified_real_places():
                numerical_value = str(embedding.im_gens()[0])
                small_maps = [var + ' |--> ' + numerical_value]
            ramified_real_places_data = 'Ramified real places: ' + str(small_maps)
            data_strings.append(ramified_real_places_data)
        output_string = str()
        for element in data_strings:
            output_string += leading_char + element + '\n'
        return output_string.rstrip()