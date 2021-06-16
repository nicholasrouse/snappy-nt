"""
This is the module that contains the class for arbitrary precision Manifolds. 

Things to consider:
    1. Building some kind of database that can be loaded to avoid repeating expensive
    computations. I haven't decided how exactly to build and access such a database
    yet.

    2. Perhaps importing functions that actually compute the various invariants from
    numerical input. I.e. make another module and put all the ugly implementation for
    computations there.
"""


#from testing import compare_against_database
import snappy, denominatorsforsnappy
from sage.all import factor, NumberField, QuaternionAlgebra, radical, cached_function
import functools
import irreducible_subgroups
import misc_functions
import field_isomorphisms
import QuaternionAlgebraNF
from collections import namedtuple


def try_various_precision(func, iterable, fail_value=None):
    """
    This should maybe be in a separate module.
    """
    for item in iterable:
        return_value = func(item)
        if return_value != fail_value:
            break
    return return_value 

PrecDegreeTuple = namedtuple("PrecDegreeTuple", ["prec", "degree"])

def fix_names(name):
    name = name.lower()
    if name == "tf" or name == "trace field":
        return "trace field"
    elif name == "itf" or name == "invariant trace field":
        return "invariant trace field"
    elif name == "qa" or name == "quaternion algebra":
        return "quaternion algebra"
    elif name == "iqa" or name =="invariant quaternion algebra":
        return "invariant quaternion algebra"


class ManifoldAP(snappy.Manifold):

    def __new__(cls, spec=None, *pargs, **kwargs):
        """
        This is to get around some Cython curiosities. Basically, even though we
        override the __init__, a superclass's __cinit__ acts like a __new__ and can mess
        up the call signatures.
        """
        return snappy.Manifold.__new__(cls, spec)
    def __init__(
        self,
        spec=None,
        delay_computations  = False,
        default_starting_prec = 1000,
        default_starting_degree = 10,
        default_max_prec = 5 * 10 ** 5,
        default_max_degree = 100,
        default_prec_increment = 5000,
        default_degree_increment = 5,
    ):
        """
        It's worth noting that we store a lot of attributes here. The reason is that a
        lot of them are somewhat expensive to compute. We could alternatively use
        sage's caching capability, but having all these class attributes should make
        them easier to store and reconstruct later. An alternative would be subclass
        things like Sage's number fields and quaternion algebras and save those
        objects.

        Unless delay_computations=True, we try to compute the main arithmetic
        arithmetic invariants with a precision that should only take at most a few
        seconds to succeed or fail. In case one needs to create a lot of these objects
        and the computations will actually meaningfully slow down whatever is being
        done, one may pass in delay_computations=True.
        """
        snappy.Manifold.__init__(self, spec)
        # We store the fields as a sage NumberField objects with a *monic* generating polynomial.
        # Perhaps subclass Sage's NumberField to store all this info?
        # The prec_record variables record whether a given prec and degree was sucessful.
        self.default_starting_prec = default_starting_prec
        self.default_starting_degree = default_starting_degree
        self.default_max_prec = default_max_prec
        self.default_max_degree = default_max_degree
        self.default_prec_increment = default_prec_increment
        self.default_degree_increment = default_degree_increment
        self._trace_field = None
        self._trace_field_numerical_root = None
        self._trace_field_generators = None
        self._trace_field_prec_record = dict()
        self._invariant_trace_field = None
        self._invariant_trace_field_numerical_root = None
        self._invariant_trace_field_generators = None
        self._invariant_trace_field_prec_record = dict()
        self._quaternion_algebra = None
        self._quaternion_algebra_prec_record = dict()
        self._invariant_quaternion_algebra = None
        self._invariant_quaternion_algebra_prec_record = dict()
        self._dict_of_prec_records = {
            "trace field": self._trace_field_prec_record,
            "invariant trace field": self._invariant_trace_field_prec_record,
            "quaternion algebra": self._quaternion_algebra_prec_record,
            "invariant quaternion algebra": self._invariant_quaternion_algebra_prec_record,
        }
        # denominators will be the empty set (i.e. set()) if there are no denominators.
        self._denominators = None
        self._denominator_residue_characteristics = None
        # This sometimes raises exceptions, but it happens in SnapPy itself.
        self._approx_trace_field_gens = self.trace_field_gens()
        self._approx_invariant_trace_field_gens = self.invariant_trace_field_gens()
        if not delay_computations:
            func = functools.partial(self.compute_arithmetic_invariants, _return_flag=True)
            try_various_precision(func, self.make_prec_degree_generator(), fail_value=False)
            self.denominators()
            

    def next_prec_and_degree(self, invariant):
        """
        This method allows us to move the logic of what the next degree and precision
        to attempt is. This might be an evolving function over time. As of Sept-24
        2020, it just takes the highest failed precision and degree pair and increases
        them by the default increment. It's returned as a named 2-tuple (prec, degree)
        for fields a named 1-tuple (prec) for the algebras. An unattractive aspect of
        the way this is currently written is that for fields, we use a named 2-tuple,
        but just an integer for the algebras because the degree has no meaning there.
        I think using a named 1-tuple is a little insane, but it makes the interfaces
        different.

        One day... I hope to know how long a particular precision, degree pair should
        take and have this return the "lowest hanging fruit" subject to some lower
        bound. But when that happens we can rewrite this method and it should be
        reflected everywhere it's used. It's possible that knowledge about the trace
        or invariant trace field could be used here for degree reasons, but right now,
        this logic lives in the respective methods.

        An assumption I'm going to make is that we want to try to compute the algebras
        with at least as much pecision as we needed for the fields. We have to return
        some degree, but it will get fixed to be that of the field in the algebra
        computation method. This is perhaps slightly undesirable, but it should be
        fine. In fact it is somewhat important that it work this way because the
        methods that compute the algebras try to compute the fields if they're not
        known using whatever precision and degree is passed in.

        The invariant needs to be passed in a string. The acceptable options are:
        "trace field", "invariant trace field", "quaternion algebra", "invariant
        quaternion algebra". The short versions "tf", "itf", "qa", "iqa" are also
        acceptable. See the fix_names function at the top level of the module for more
        information.
        """
        invariant = fix_names(invariant)
        record = self._dict_of_prec_records[invariant]
        if invariant == 'trace field' or invariant == 'invariant trace field':
            if not record:
                return PrecDegreeTuple(self.default_starting_prec, self.default_starting_degree)
            if True in record.values():
                smallest_successful_prec = min(
                    [pair.prec for pair in record if record[pair]]
                )
                smallest_successful_degree = min(
                    [pair.degree for pair in record if record[pair]]
                )
                newpair = PrecDegreeTuple(
                    smallest_successful_prec, smallest_successful_degree
                )
            else:
                # Given the logic so far, checking if not record[pair] is superfluous here, but
                # I think it's better to leave it as is in case things change later.
                largest_failed_prec = max(
                    [pair.prec for pair in record if not record[pair]]
                )
                largest_failed_degree = max(
                    [pair.degree for pair in record if not record[pair]]
                )
                newpair = PrecDegreeTuple(
                largest_failed_prec + self.default_prec_increment,
                largest_failed_degree + self.default_degree_increment,
            )
            return newpair
        if invariant == 'quaternion algebra' or invariant == 'invariant quaternion algebra':
            field = (
                "trace field"
                if invariant == "quaternion algebra"
                else "invariant trace field"
            )
            field_prec_deg = self.next_prec_and_degree(field)
            field_prec = field_prec_deg.prec
            if not record:
                return field_prec
            else:
                if True in record.values(): 
                    return min([prec for prec in record if record[prec]])
                else:
                    largest_failed_prec = max([prec for prec in record if not record[prec]])
                    return max(largest_failed_prec + self.default_prec_increment, field_prec)
                
            

    def has_two_torsion_in_homology(self):
        """
        Returns True if there is two-torsion in homology and False if not. This doesn't
        really need arbitrary precision, but it hopefully makes some other code
        cleaner.

        Obviously this basically factors an integer to check if it's even, which is not
        really optimal, so we should do something about this at some point.
        """
        homology = self.homology()
        elementary_divisors = homology.elementary_divisors()
        elementary_divisors = [divisor for divisor in elementary_divisors if divisor != 0]
        for divisor in elementary_divisors:
            if divisor % 2 == 0:
                return True
        return False

    def defining_function(self, prec):
        return snappy.snap.polished_holonomy(self, bits_prec=prec)

    def trace_field(
        self,
        prec=None,
        degree=None,
        be_smart=True,
        verbosity=False,
        _force_compute=False,
    ):
        """
        If be_smart is False, the function just tries to compute the trace field using
        whatever precision and degree is passed in. In particular when be_smart=False,
        the _force_compute parameter is superfluous.

        When be_smart is True, we use our records for the trace field to pick a degree
        and precision to try. I.e. whatever prec and degrees might have been passed in
        are ignored. The degree might change based on the next paragraph though.

        When be_smart is True and the invariant trace field is known and the orbifold
        has no 2-torsion in homology, then the trace field is equal to the invariant
        trace field. However, in this case, exact generators for the trace field will
        not be computed, and these generators are necessary (as far as I know) to
        compute the primes at which the Kleinian group has nonintegral trace. When
        _force_compute=True, these generators will be found assuming the method
        succeeds in finding the trace field at all. If there is 2-torsion in homology,
        it's possible that trace field is a proper extension of degree of the
        invariant trace field. In this case, we do still get some knowledge about the
        degree, namely that it's at least as large as that of the invariant trace field.

        The verbosity argument prints some information as the method proceeds. This
        can be useful for large calculations.

        If _force_compute is True, then the method won't stop when it finds the trace
        field; it will be sure to find the generators as well.

        Last updated: Sept-24 2020
        """
        if prec == None: prec = self.default_starting_prec
        if degree == None: degree = self.default_starting_degree
        if self._trace_field and not _force_compute:
            return self._trace_field
        if be_smart:
            prec, degree = self.next_prec_and_degree("trace field")
            if self._invariant_trace_field:
                itf_deg = self._invariant_trace_field.degree()
                if not self.has_two_torsion_in_homology():
                    self._trace_field_numerical_root = self._invariant_trace_field_numerical_root
                    self._trace_field = self._invariant_trace_field
                    self._trace_field_generators = self._invariant_trace_field_generators
                    if verbosity:
                        print(
                            "Found invariant trace field and no 2-torsion in homology."
                        )
                    if not _force_compute:
                        return self._trace_field
                    degree = self._invariant_trace_field.degree()
                else:
                    if degree < 2*itf_deg: degree = 2*itf_deg
        if verbosity:
            print(
                f"Trying to compute trace field with precision={prec} and degree={degree}."
            )
        exact_field_data = self._approx_trace_field_gens.find_field(
            prec=prec, degree=degree, optimize=True
        )
        # This will override previous calculations with same prec and degree.
        # It's unclear if we want this behavior.
        self._trace_field_prec_record[PrecDegreeTuple(prec, degree)] = bool(
            exact_field_data
        )
        if exact_field_data is not None:
            if verbosity:
                print("\tTrace field found.")
            self._trace_field = exact_field_data[0]
            self._trace_field_numerical_root = exact_field_data[1]  # An AAN
            self._trace_field_generators = exact_field_data[2]
        else:
            if verbosity:
                print("\tTrace field not found.")
        return self._trace_field

    def invariant_trace_field(
        self,
        prec=None,
        degree=None,
        be_smart=True,
        verbosity=False,
        _force_compute=False,
    ):
        """
        This now should work similarly to compute_trace_field method. There is of
        course slightly different logic for using knowledge of 2-torsion in homology.
        In particular, if the (noninvariant) trace field is known and there is no
        2-torsion in homology, then the fields are the same.
        Last updated: Aug-29 2020
        """
        if prec == None: prec = self.default_starting_prec
        if degree == None: degree = self.default_starting_degree 
        if self._invariant_trace_field and not _force_compute:
            return self._invariant_trace_field
        if be_smart:
            prec, degree = self.next_prec_and_degree("invariant trace field")
            if self._trace_field:
                tf_deg = self._trace_field.degree()
                if (
                    not self.has_two_torsion_in_homology()
                    or self._trace_field.degree() % 2 == 1
                ):
                    self._invariant_trace_field = self._trace_field
                    self._invariant_trace_field_numerical_root = self._trace_field_numerical_root
                    self._invariant_trace_field_generators = self._trace_field_generators
                    if verbosity:
                        print("Found trace field, and it conincides with invariant trace field.")
                    if not _force_compute:
                        return self._invariant_trace_field
                    degree = self._invariant_trace_field.degree()
                else:
                    if degree >= tf_deg:
                        degree = tf_deg
        if verbosity:
            print(
                f"Trying to compute invariant trace field with precision={prec} and degree={degree}."
            )
        exact_field_data = self._approx_invariant_trace_field_gens.find_field(
            prec=prec, degree=degree, optimize=True
        )
        self._invariant_trace_field_prec_record[PrecDegreeTuple(prec, degree)] = bool(
            exact_field_data
        )
        if exact_field_data is not None:
            if verbosity:
                print("\tInvariant trace field found.")
            self._invariant_trace_field = exact_field_data[0]
            self._invariant_trace_field_numerical_root = exact_field_data[1]  # An AAN
            self._invariant_trace_field_generators = exact_field_data[2]
        else:
            if verbosity:
                print("\tInvariant trace field not found.")
        return self._invariant_trace_field

    def approximate_trace(self, word):
        """
        Given a word in the generators for the fundamental group, returns an
        ApproximateAlgebraicNumber which is the trace of that element in SL_2(CC). This
        perhaps shouldn't really be a method of the manifold but rather of the group,
        but we can perhaps change this later.
        """

        def trace_defining_func(prec):
            approximate_group = self.defining_function(prec=prec)
            approximate_matrix = approximate_group(word)
            return approximate_matrix.trace()

        return snappy.snap.find_field.ApproximateAlgebraicNumber(trace_defining_func)

    def compute_approximate_hilbert_symbol(self, power=1, epsilon_coefficient=10):
        """
        Somewhat cumbersomly computes a Hilbert symbol as a pair of
        ApproximateAlgebraicNumbers. It's possible I should use the class
        ListOfApproximateAlgebraicNumbers instead.

        power=1 is for noninvariant quaternion algebra, power=2 is for invariant
        quaternion algebra. More generally power=n will give approximate algebraic
        numbers for traces of elements g,h where g is not parabolic and <g,h> is
        irreducible where g,h belong to G^{(n)}, where G is the fundamental group
        (or really the Kleinian group that gives rise to the orbifold). I know of no
        use for any power beyond 2, though.

        Possible improvement: try using the exact expression that are known to trace
        field before getting new elements as AANs. Might not work though since probably
        need a commutator somewhere.

        Last updated: Aug-29 2020
        """
        (word1, word2) = irreducible_subgroups.find_hilbert_symbol_words(
            self.defining_function(prec=self.default_starting_prec), power=power, epsilon_coefficient=epsilon_coefficient
        )
        first_entry = self.approximate_trace(
            word1
        ) ** 2 - snappy.snap.find_field.ApproximateAlgebraicNumber(4)
        commutator_word = misc_functions.commutator_of_words(word1, word2)
        second_entry = self.approximate_trace(
            commutator_word
        ) - snappy.snap.find_field.ApproximateAlgebraicNumber(2)
        return (first_entry, second_entry)

    def quaternion_algebra(
        self, prec=None, be_smart=True, verbosity=False, _force_compute=False, compute_ramification=True, **kwargs
    ):
        """
        This method won't try to compute the trace field if it isn't known. The
        precision variable is used for expressing the Hilbert symbol entries in terms
        of the elements of the trace field.

        The be_smart parameter is used to get a smarter precision based on previous
        computations.

        There are no special kwargs; it just allows for passing in the same argument
        signature as the other invarinats (e.g. compute_trace_field).

        Possible refactor: Just have one method for computing quaternion algebras from
        ApproximateAlgebraicNumbers. In that case, probably easiest to make another
        module wherein we subclass ApproximateAlgebraicNumber. This could simplify
        other things as well. On the other hand, from a mathemtical perspective, one
        can more or less interchangably work with the Kleinian group or the orbifold,
        so there's not a compelling mathematical reason for say the quaternion algebra
        to be associated to the orbifold rather than the group. From a coding
        perspective, I think I would prefer to hide things like finding Hilbert symbols
        in terms of approximate algebraic integers in another module (probably the
        one that deals with the group more directly) because these are not really of
        great mathematical interest to the orbifolds in the way that the trace fields
        and quaternion algebras are.

        For now though most everything is a ManifoldAP method.
        """
        if prec == None: prec = self.default_starting_prec
        if self._quaternion_algebra and not _force_compute:
            return self._quaternion_algebra
        if not self._trace_field:
            if verbosity:
                failure_message = (
                    "Trace field not known. It can be "
                    "computed with the compute_trace_field method. It's "
                    "possible that one needs high precision or degree to find the field, "
                    "though."
                )
                print(failure_message)
            return None
        if be_smart:
            prec = self.next_prec_and_degree("quaternion algebra")
        primitive_element = self._trace_field_numerical_root  # An AAN
        epsilon_coefficient = 10
        while True:
            
            (
                approx_first_entry,
                approx_second_entry,
            ) = self.compute_approximate_hilbert_symbol(power=1, epsilon_coefficient=epsilon_coefficient)
            first_entry = primitive_element.express(approx_first_entry, prec=prec)
            second_entry = primitive_element.express(approx_second_entry, prec=prec)
            if first_entry == 0 or second_entry == 0:
                epsilon_coefficient *= 10
            else:
                break
        self._quaternion_algebra_prec_record[prec] = bool(first_entry and second_entry)
        if first_entry == None or second_entry == None:
            if verbosity:
                print('Failed to find quaternion algebra.')
            return None
        else:
            if verbosity: print('Found quaternion algebra.')
            first_entry, second_entry = first_entry(self._trace_field.gen()), second_entry(self._trace_field.gen())
            self._quaternion_algebra = QuaternionAlgebraNF.QuaternionAlgebraNF(
                self._trace_field, first_entry, second_entry, compute_ramification=compute_ramification
            )
        return self._quaternion_algebra

    def invariant_quaternion_algebra(
        self, prec=None, be_smart=True, verbosity=False, _force_compute=False, compute_ramification=True, **kwargs
    ):
        """
        See docstring for compute_quaterion_algebra_fixed_prec. Should try to refactor this
        somehow since it's so similar to the one for the noninvariant quaternion algebra.

        Last updated: Aug-29 2020
        """
        if prec == None: prec = self.default_starting_prec
        if self._invariant_quaternion_algebra and not _force_compute:
            return self._invariant_quaternion_algebra
        if not self._invariant_trace_field:
            if verbosity:
                failure_message = (
                    "Invariant trace field not known. It can be "
                    "computed with the compute_invariant_trace_field method. It's "
                    "possible that one needs high precision or degree to find the field, "
                    "though."
                )
                print(failure_message)
            return None
        if be_smart:
            prec = self.next_prec_and_degree('invariant quaternion algebra')
        #degree = self._invariant_trace_field.degree()
        primitive_element = self._invariant_trace_field_numerical_root  # An AAN
        epsilon_coefficient = 10
        while True:
            (
                approx_first_entry,
                approx_second_entry,
            ) = self.compute_approximate_hilbert_symbol(power=2, epsilon_coefficient=epsilon_coefficient)
            first_entry = primitive_element.express(approx_first_entry, prec=prec)
            second_entry = primitive_element.express(approx_second_entry, prec=prec)
            if verbosity:
                print("epsilon_coefficient=", epsilon_coefficient, "\nfirst_entry=", first_entry, "\nsecond_entry=", second_entry)
                print(first_entry == self.invariant_trace_field(0))
            if first_entry == 0 or second_entry == 0:
                epsilon_coefficient *= 10
            else:
                break
        self._invariant_quaternion_algebra_prec_record[prec] = bool(first_entry and second_entry)
        if first_entry == None or second_entry == None:
            if verbosity: print('Failed to find invariant quaternion algebra.')
            return None
        else:
            if verbosity: print('Found invariant quaternion algebra.')
            first_entry, second_entry = first_entry(self._invariant_trace_field.gen()), second_entry(self._invariant_trace_field.gen())
            self._invariant_quaternion_algebra = QuaternionAlgebraNF.QuaternionAlgebraNF(
                self._invariant_trace_field, first_entry, second_entry, compute_ramification=compute_ramification
            )
        if compute_ramification:
            if self._invariant_quaternion_algebra:
                discriminant_list = list(
                    self._invariant_quaternion_algebra.discriminant().factor()
                )
                self._invariant_quaternion_algebra_ramified_places = [
                    ideal for (ideal, multiplicity) in discriminant_list
                ]
                self._invariant_quaternion_algebra_ramified_places_residue_characteristics = list(
                    {
                        radical(place.absolute_norm())
                        for place in self._invariant_quaternion_algebra_ramified_places
                    }
                )
                self._invariant_quaternion_algebra_ramified_places_residue_characteristics.sort()
        return self._invariant_quaternion_algebra
    
    def denominators(self, verbosity=False, **kwargs):
        """
        This function incidentally computes the residue characteristics of the
        denominators for easy access later.

        It's also worth pointing out that the denominators are returned as a set of
        ideals of a number field. This is different from the behavior in the
        denominatorsforsnappy module that just returns the residue characteristics. We
        could add this as an optional argument at some point though.

        Recall the convention that self._denominators is None if they haven't been
        computed but set() if they have been to computed to be the empty set.
        """
        if self._denominators or self._denominators == set():
            return self._denominators
        if not self._trace_field_generators:
            if verbosity:
                failure_message = (
                    "Generators for the (noninvariant) trace field are not known. "
                    "They can be computed using the trace_field method, but "
                    "make sure to use _force_compute=True to make sure the generators "
                    "are found."
                )
                print(failure_message)
            return None
        denominator_ideals = {
            element.denominator_ideal() for element in self._trace_field_generators
        }
        prime_ideals = set()
        for ideal in denominator_ideals:
            factorization = ideal.factor()
            for element in factorization:
                prime_ideals.add(element[0])
        self._denominators = prime_ideals
        norms = {ideal.absolute_norm() for ideal in prime_ideals}
        self._denominator_residue_characteristics = denominatorsforsnappy.find_prime_factors_in_a_set(
            norms
        )
        return prime_ideals
    
    def denominator_residue_characteristics(self, verbosity=False, **kwargs):
        if self._denominators is None:
            self.denominators()
        if self._denominator_residue_characteristics is None:
            prime_ideals = self._denominators
            norms = {ideal.absolute_norm() for ideal in prime_ideals}
            self._denominator_residue_characteristics = denominatorsforsnappy.find_prime_factors_in_a_set(norms)
        return self._denominator_residue_characteristics

    def make_prec_degree_generator(
        self,
        starting_prec=None,
        starting_degree=None,
        prec_increment=None,
        degree_increment=None,
        max_prec=None,
        max_degree=None,
        be_smart=True,
        verbosity=False,
        _force_compute=False
    ):
        """
        This makes a generator the output of which can be passed into some functions to
        compute invariants. Most importantly, the generator itself can
        be passed into the try_various_precision top level function to try several to
        compute the invariants using different precision and degree parameters.
        Everything should be robust enough that even passing in degrees to methods that
        don't need them (e.g. to computing quaternion algbras) shouldn't break anything.

        This function is also used together with try_various_precision at object
        initialization time unless delay_computations=True.
        """
        if starting_prec == None: starting_prec = self.default_starting_prec
        if starting_degree == None: starting_degree = self.default_starting_degree
        if prec_increment == None: prec_increment = self.default_prec_increment
        if degree_increment == None: degree_increment = self.default_degree_increment
        if max_prec == None: max_prec = self.default_max_prec
        if max_degree == None: max_degree = self.default_max_degree
        def gen():
            """
            This can be neater perhaps?
            """
            prec, degree = starting_prec, starting_degree
            yield {'prec':prec, 'degree':degree, 'be_smart': be_smart, 'verbosity': verbosity, "_force_compute": _force_compute}
            while prec < max_prec and degree <= max_degree:
                prec = min(prec + prec_increment, max_prec)
                degree = min(degree + degree_increment, max_degree)
                yield {'prec':prec, 'degree':degree, 'be_smart': be_smart, 'verbosity': verbosity, "_force_compute": _force_compute}
                if prec == max_prec and degree == max_degree: break
        return gen()


    def compute_arithmetic_invariants(
        self,
        prec=None,
        degree=None,
        be_smart=True,
        verbosity=False,
        _force_compute=False,
        _return_flag=False,
        print_results=False
    ):
        """
        This tries to compute the four basic arithmetic invariants: the two trace
        fields and the two quaternion algebras.

        It will also try to compute the other invariants to fill out all the attributes
        of the instance. Right now it's called upon creation of a ManifoldAP instance,
        but this can be disabled with a keyword argument when a ManifoldAP object is
        initialized.

        Warning: if you pass _force_compute in as True, then this method will try to
        recompute the arithmetic invariants. If you combine this with vary_precision,
        there is a good chance the computation will spend a long time recomputing known
        invariants.

        This function returns True if all four of the main arithmetic invariants are
        found and False if not. The reason for this is to use with try_various_precision.
        This might be overcoupled though, so we should think about ways to deal with
        this. This dodges the need for needing to create some partial functions, but
        maybe that's just the answer.
        """
        if prec == None: prec = self.default_starting_prec
        if degree == None: degree = self.default_starting_degree
        arguments = {'prec': prec, 'degree': degree, 'be_smart': be_smart, 'verbosity': verbosity}
        invariant_method_pairs = [
            (self._trace_field, ManifoldAP.trace_field),
            (self._quaternion_algebra, ManifoldAP.quaternion_algebra),
            (self._invariant_trace_field, ManifoldAP.invariant_trace_field),
            (
                self._invariant_quaternion_algebra,
                ManifoldAP.invariant_quaternion_algebra,
            ),
        ]
        for (invariant, method) in invariant_method_pairs:
            method(self,**arguments)
        if self._trace_field_generators:
            ManifoldAP.denominators(self, **arguments)
        if _return_flag:
            return bool(self._trace_field and self._quaternion_algebra and self._invariant_trace_field and self._invariant_quaternion_algebra)
        if print_results:
            self.p_arith()
            
    
    def is_arithmetic(self):
        """
        This checks whether the manifold (really the Kleinian group) is arithmetic.
        It doesn't itself explicitly compute the necessary invariants if they aren't
        already known.
        For why this works, see MR Theorem 8.3.2 pp.261-262.
        This could be a one-liner, but I think it's clearer this way.
        """
        (
            number_of_real_places,
            number_of_complex_places,
        ) = self._invariant_trace_field.signature()
        number_of_ramified_real_places = len(
            misc_functions.ramified_real_places(self._invariant_quaternion_algebra)
        )
        return (
            number_of_ramified_real_places == number_of_real_places
            and number_of_complex_places == 1
            and self._denominators == set()
        )

    def p_arith(self):
        """
        This is so named for the common use case of typing p arith in snap to get the
        arithmetic invariants.

        This function is probably a good argument for subclassing a lot of our objects
        and defining __str__ methods in those classes.
        """
        print("Orbifold name:", self)
        print("Volume:", self.volume())
        if self._trace_field:
            print("Trace field:", self._trace_field)
            print("\t Signature:", self._trace_field.signature())
            print("\t Discriminant:", self._trace_field.discriminant())
            if self._quaternion_algebra:
                print('Quaternion algebra:')
                print(self._quaternion_algebra.ramification_string(leading_char='\t'))
            else:
                print("Quaternion algebra not found.")
        else:
            print("Trace field not found.")
        if self._invariant_trace_field:
            print("Invariant Trace field:", self._invariant_trace_field)
            print("\t Signature:", self._invariant_trace_field.signature())
            print("\t Discriminant:", self._invariant_trace_field.discriminant())
            if self._invariant_quaternion_algebra:
                print('Invariant quaternion algebra:')
                print(self._invariant_quaternion_algebra.ramification_string(leading_char='\t'))
            else:
                print("Invariant quaternion algebra not found.")
        else:
            print("Invariant trace field not found.")
        if self._denominators is None:
            print("Denominators not found (trace field probably not computed)")
        else:
            print("Integer traces:", not bool(self._denominators))
            if len(self._denominators) >= 1:
                print("\t Denominator ideals:", self._denominators)
                print(
                    "\t Denominator Residue Characteristics:",
                    self._denominator_residue_characteristics,
                )
        if self._trace_field and self._invariant_quaternion_algebra:
            print("Arithmetic:", bool(self.is_arithmetic()))

    def delete_arithmetic_invariants(self):
        """
        This method sets all the invariants back to their starting values (usually
        None). I don't love having a method like this, but when we Dehn fill we need
        to forget all this data since other methods try to use known information about
        the various invariants to compute other ones. I'll try to find a better way to
        do this in a later version.

        Last updated: Sept-9 2020
        """
        self._trace_field = None
        self._trace_field_numerical_root = None
        self._trace_field_generators = None
        self._trace_field_prec_record = dict()
        self._invariant_trace_field = None
        self._invariant_trace_field_numerical_root = None
        self._invariant_trace_field_generators = None
        self._invariant_trace_field_prec_record = dict()
        self._quaternion_algebra = None
        self._quaternion_algebra_prec_record = dict()
        self._invariant_quaternion_algebra = None
        self._invariant_quaternion_algebra_prec_record = dict()
        self._dict_of_prec_records = {
            "trace field": self._trace_field_prec_record,
            "invariant trace field": self._invariant_trace_field_prec_record,
            "quaternion algebra": self._quaternion_algebra_prec_record,
            "invariant quaternion algebra": self._invariant_quaternion_algebra_prec_record,
        }
        self._denominators = None
        self._denominator_residue_characteristics = None
        self._approx_trace_field_gens = self.trace_field_gens()
        self._approx_invariant_trace_field_gens = self.invariant_trace_field_gens()

    def dehn_fill(self, filling_data, which_cusp=None, delay_computations=False):
        """
        This performs dehn surgery "in place," i.e. it doesn't return a new manifold
        but rather changes self. This, of course, necessitates recomputing all the
        arithmetic invariants, so we have to override SnapPy's method to make sure
        we clean out the invariants.
        """
        snappy.Manifold.dehn_fill(self, filling_data, which_cusp=which_cusp)
        self.delete_arithmetic_invariants()
        if not delay_computations:
            func = functools.partial(self.compute_arithmetic_invariants, _return_flag=True)
            try_various_precision(func, self.make_prec_degree_generator(), fail_value=False)
    
    def compare_arithmetic_invariants(self, other):
        """
        This takes two ManifoldAP's and computes whether they have isomorphic trace
        fields, invariant trace fields, quaternion algebras, invariant quaternion
        algebras, and denominators. It does check whether the numerical roots of fields
        agree. This function is primarily useful for testing various revisions of this
        package against known correct computations. To that end, this function might get
        moved out of this module into some module designed for testing. Note that
        checking whether number fields are isomorphic can be an expensive calculation if
        the fields are not distinguished by simple invariants like degree or
        discriminant.

        We use a custom function in the field_isomorphism module to check that two
        number fields are isomorphic and that their distinguished complex places
        conincide.

        One day other might be able to be a normal Manifold from SnapPy whence we
        compute its arithmetic invariants and compare them to self.
        """
        arith_dict = dict()
        arith_dict['trace field'] = field_isomorphisms.same_subfield_of_CC(self._trace_field, other._trace_field)
        arith_dict['invariant trace field'] = field_isomorphisms.same_subfield_of_CC(self._invariant_trace_field, other._invariant_trace_field)
        if arith_dict['trace field']:
            arith_dict['quaternion algebra'] = self._quaternion_algebra.is_isomorphic(other._quaternion_algebra)
            # This should be computed and saved when we find whether the fields are isomorphic in the first place.
            iso = field_isomorphisms.isomorphisms_between_number_fields(self._trace_field, other._trace_field)[0]
            other_denominators = {iso(ideal) for ideal in other._denominators}
            arith_dict['denominators'] = (self._denominators == other_denominators)
        else: 
            arith_dict['quaternion algebra'] = False
            # When the trace fields differ, the convention we take is that the
            # denominators are the same if and only both orbifolds have integral traces.
            arith_dict['denominators'] = True if (self._denominators == set() and other._denominators == set()) else False
        if arith_dict['invariant trace field']:
            arith_dict['invariant quaternion algebra'] = self._invariant_quaternion_algebra.is_isomorphic(other._invariant_quaternion_algebra)
        else: arith_dict['invariant quaternion algebra'] = False
        return arith_dict
        
    def has_same_arithmetic_invariants(self, other):
        arith_dict = self.compare_arithmetic_invariants(other)
        return not (False in arith_dict.values())
    
    def identify(self):
        """
        I'm overriding this for now because for some reason I get an error that the
        triangulation is empty when I call it on ManifoldAP objects. It seems to work
        fine for snappy ones. My guess is that it's related to the __new__ Cython issue,
        but I'll try to fix it later.
        """
        temp_mfld = snappy.Manifold(str(self))
        return temp_mfld.identify()