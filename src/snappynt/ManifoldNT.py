"""
This is the module that contains the class for arbitrary precision Manifolds.

Things to consider:

    1. Perhaps importing functions that actually compute the various invariants from
    numerical input. I.e. make another module and put all the ugly implementation for
    computations there.
"""


from collections import namedtuple

import snappy
from sage.all import ZZ, log

from . import (
    QuaternionAlgebraNF,
    denominatorsforsnappy,
    field_isomorphisms,
    irreducible_subgroups,
    misc_functions,
)


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
    elif name == "iqa" or name == "invariant quaternion algebra":
        return "invariant quaternion algebra"


class ManifoldNT:
    def __init__(
        self,
        spec=None,
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
        self._snappy_mfld = snappy.Manifold(spec)
        # The fields are sage NumberField objects with a *monic* generating polynomial.
        # Perhaps subclass Sage's NumberField to store all this info?
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
        self._approx_trace_field_gens = self._snappy_mfld.trace_field_gens()
        if self.homology_two_rank() == 0:
            self._approx_invariant_trace_field_gens = self._approx_trace_field_gens
        else:
            self._approx_invariant_trace_field_gens = (
                self._snappy_mfld.invariant_trace_field_gens()
            )
        """
        if not delay_computations:
            func = functools.partial(
                self.compute_arithmetic_invariants, _return_flag=True
            )
            try_various_precision(
                func, self.make_prec_degree_generator(), fail_value=False
            )
            self.denominators()
        """

    def __getattr__(self, attr):
        return getattr(self._snappy_mfld, attr)

    def __str__(self):
        return str(self._snappy_mfld)

    def __repr__(self):
        return self._snappy_mfld.__repr__()

    # The _approx_trace_field_gens are unpicklable, so we have these
    def __getstate__(self):
        # Maybe move this, so we keep track of unpicklable object elsewhere.
        unpicklable = ["_approx_trace_field_gens", "_approx_invariant_trace_field_gens"]
        state = self.__dict__.copy()
        for attr in unpicklable:
            del state[attr]
        # SnapPy has some issues with pickle.
        # E.g. m003(-2,3) becomes m003(254,3) upon unpickling.
        state["snappy_mfld_name"] = str(self._snappy_mfld)
        del state["_snappy_mfld"]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        # We reconstruct the wrapped snappy manifold each time since it's cheap and has
        # some issues being pickled.
        self._snappy_mfld = snappy.Manifold(state["snappy_mfld_name"])
        self._approx_trace_field_gens = self._snappy_mfld.trace_field_gens()
        self._approx_invariant_trace_field_gens = (
            self._snappy_mfld.invariant_trace_field_gens()
        )

    def _arithmetic_invariants_known(self):
        basic_invariants_known = all(
            [
                self._trace_field,
                self._invariant_trace_field,
                self._quaternion_algebra,
                self._invariant_quaternion_algebra,
            ]
        )
        denominators_known = (
            self._denominators is not None
        )  # Denominators are allowed to be set(), and bool(set()) == False.
        return basic_invariants_known and denominators_known

    def defining_function(self, prec):
        return snappy.snap.polished_holonomy(self, bits_prec=prec)

    def next_prec_and_degree(
        self,
        invariant,
        starting_prec=1000,
        starting_degree=20,
        prec_increment=5000,
        degree_increment=5,
    ):
        """
        This method allows us to move the logic of what the next degree and precision
        to attempt is. This might be an evolving function over time. As of Sept-24
        2020, it just takes the highest failed precision and degree pair and increases
        them by the default increment. It's returned as a named 2-tuple (prec, degree)
        for fields an integer (i.e. not a namedtuple) for the algebras.

        An assumption I'm going to make is that we want to try to compute the algebras
        with at least as much pecision as we needed for the fields. We have to return
        some degree, but it will get fixed to be that of the field in the algebra
        computation method. This is perhaps slightly undesirable, but it should be
        fine. In fact it is somewhat important that it work this way because the
        methods that compute the algebras try to compute the fields if they're not
        known using whatever precision and degree is passed in.

        The invariant needs to be passed in as a string. The acceptable options are:
        "trace field", "invariant trace field", "quaternion algebra", "invariant
        quaternion algebra". The short versions "tf", "itf", "qa", "iqa" are also
        acceptable. See the fix_names function at the top level of the module for more
        information.
        """
        invariant = fix_names(invariant)
        record = self._dict_of_prec_records[invariant]
        if invariant == "trace field" or invariant == "invariant trace field":
            if record == dict():
                return PrecDegreeTuple(starting_prec, starting_degree)
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
                largest_failed_prec = max(
                    [pair.prec for pair in record if not record[pair]]
                )
                largest_failed_degree = max(
                    [pair.degree for pair in record if not record[pair]]
                )
                newpair = PrecDegreeTuple(
                    largest_failed_prec + prec_increment,
                    largest_failed_degree + degree_increment,
                )
                if (
                    invariant == "trace field"
                    and self._invariant_trace_field is not None
                ):
                    v = ZZ.valuation(2)
                    biggest_valuation = min(
                        v(self.homology_two_rank()),
                        v(self._invariant_trace_field.degree()),
                    )
                    possible_degrees = [
                        self._invariant_trace_field.degree() / k
                        for k in (2 ** j for j in range(biggest_valuation + 1))
                    ]
                    newpair.degree = max(
                        deg for deg in possible_degrees if deg < newpair.degree
                    )
                if (
                    invariant == "invariant trace field"
                    and self._trace_field is not None
                ):
                    v = ZZ.valuation(2)
                    possible_degrees = [
                        self._trace_field.degree() * k
                        for k in (2 ** j for j in range(self.homology_two_rank() + 1))
                    ]
                    newpair.degree = max(
                        deg for deg in possible_degrees if deg < newpair.degree
                    )
            return newpair
        if (
            invariant == "quaternion algebra"
            or invariant == "invariant quaternion algebra"
        ):
            field = (
                "trace field"
                if invariant == "quaternion algebra"
                else "invariant trace field"
            )
            field_prec = self.next_prec_and_degree(field).prec
            if not record:
                return field_prec
            else:
                if True in record.values():
                    return min([prec for prec in record if record[prec]])
                else:
                    largest_failed_prec = max(
                        [prec for prec in record if not record[prec]]
                    )
                    return max(largest_failed_prec + prec_increment, field_prec)

    def homology_two_rank(self):
        """
        Returns the number of Z/2Z factors in the first homology group.
        """
        homology = self.homology()
        elementary_divisors = homology.elementary_divisors()
        two_factor = [
            divisor
            for divisor in elementary_divisors
            if divisor != 0 and divisor % 2 == 0
        ].pop()
        return log(two_factor, 2)

    def trace_field(
        self,
        prec=None,
        degree=None,
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
        if self._trace_field and prec is None and degree is None:
            return self._trace_field
        if prec is None:
            prec = self.next_prec_and_degree("trace field").prec
        if degree is None:
            degree = self.next_prec_and_degree("trace field").degree
        exact_field_data = self._approx_trace_field_gens.find_field(
            prec=prec, degree=degree, optimize=True
        )
        # This will override previous calculations with same prec and degree.
        # It's unclear if we want this behavior.
        self._trace_field_prec_record[PrecDegreeTuple(prec, degree)] = bool(
            exact_field_data
        )
        if exact_field_data is not None:
            self._trace_field = exact_field_data[0]
            self._trace_field_numerical_root = exact_field_data[1]  # An AAN
            self._trace_field_generators = exact_field_data[2]
            if self._invariant_trace_field is None and self.homology_two_rank() == 0:
                self._invariant_trace_field_prec_record[
                    PrecDegreeTuple(prec, degree)
                ] = True
                self._invariant_trace_field = exact_field_data[0]
                self._invariant_trace_field_numerical_root = exact_field_data[1]
                self._invariant_trace_field_generators = exact_field_data[2]
        else:
            return None
        return self._trace_field

    def invariant_trace_field(
        self,
        prec=None,
        degree=None,
    ):
        """
        This now should work similarly to compute_trace_field method. There is of
        course slightly different logic for using knowledge of 2-torsion in homology.
        In particular, if the (noninvariant) trace field is known and there is no
        2-torsion in homology, then the fields are the same.
        Last updated: Aug-29 2020
        """
        if self._invariant_trace_field and prec is None and degree is None:
            return self._invariant_trace_field
        if prec is None:
            prec = self.next_prec_and_degree("invariant trace field").prec
        if degree is None:
            degree = self.next_prec_and_degree("invariant trace field").degree
        exact_field_data = self._approx_invariant_trace_field_gens.find_field(
            prec=prec, degree=degree, optimize=True
        )
        self._invariant_trace_field_prec_record[PrecDegreeTuple(prec, degree)] = bool(
            exact_field_data
        )
        if exact_field_data is not None:
            self._invariant_trace_field = exact_field_data[0]
            self._invariant_trace_field_numerical_root = exact_field_data[1]  # An AAN
            self._invariant_trace_field_generators = exact_field_data[2]
            if self._trace_field is None and self.homology_two_rank() == 0:
                self._trace_field_prec_record[PrecDegreeTuple(prec, degree)] = True
                self._trace_field = exact_field_data[0]
                self._trace_field_numerical_root = exact_field_data[1]
                self._trace_field_generators = exact_field_data[2]
        else:
            return None
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
            self.defining_function(prec=self.default_starting_prec),
            power=power,
            epsilon_coefficient=epsilon_coefficient,
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
        self,
        prec=None,
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

        For now though most everything is a ManifoldNT method.
        """
        if self._quaternion_algebra and prec is None:
            return self._quaternion_algebra
        if prec is None:
            prec = self.next_prec_and_degree("quaternion algebra")
        if not self._trace_field:
            self._trace_field = self.trace_field(prec=prec)
            if self._trace_field is None:
                return None
        primitive_element = self._trace_field_numerical_root  # An AAN
        epsilon_coefficient = 10
        while True:
            (
                approx_first_entry,
                approx_second_entry,
            ) = self.compute_approximate_hilbert_symbol(
                power=1, epsilon_coefficient=epsilon_coefficient
            )
            first_entry = primitive_element.express(approx_first_entry, prec=prec)
            second_entry = primitive_element.express(approx_second_entry, prec=prec)
            if first_entry == 0 or second_entry == 0:
                epsilon_coefficient *= 10
            else:
                break
        self._quaternion_algebra_prec_record[prec] = bool(first_entry and second_entry)
        if first_entry is None or second_entry is None:
            return None
        else:
            first_entry, second_entry = first_entry(
                self._trace_field.gen()
            ), second_entry(self._trace_field.gen())
            self._quaternion_algebra = QuaternionAlgebraNF.QuaternionAlgebraNF(
                self._trace_field,
                first_entry,
                second_entry,
            )
        return self._quaternion_algebra

    def invariant_quaternion_algebra(self, prec=None):
        """
        See docstring for compute_quaterion_algebra_fixed_prec. Should try to refactor this
        somehow since it's so similar to the one for the noninvariant quaternion algebra.

        Last updated: Aug-29 2020
        """
        if self._invariant_quaternion_algebra and prec is None:
            return self._invariant_quaternion_algebra
        if prec is None:
            prec = self.next_prec_and_degree("invariant quaternion algebra")
        if not self._invariant_trace_field:
            self._invariant_trace_field = self.invariant_trace_field(prec=prec)
            if self._invariant_trace_field is None:
                return None
        primitive_element = self._invariant_trace_field_numerical_root  # An AAN
        epsilon_coefficient = 10
        while True:
            (
                approx_first_entry,
                approx_second_entry,
            ) = self.compute_approximate_hilbert_symbol(
                power=2, epsilon_coefficient=epsilon_coefficient
            )
            first_entry = primitive_element.express(approx_first_entry, prec=prec)
            second_entry = primitive_element.express(approx_second_entry, prec=prec)
            if first_entry == 0 or second_entry == 0:
                epsilon_coefficient *= 10
            else:
                break
        self._invariant_quaternion_algebra_prec_record[prec] = bool(
            first_entry and second_entry
        )
        if first_entry is None or second_entry is None:
            return None
        else:
            first_entry, second_entry = first_entry(
                self._invariant_trace_field.gen()
            ), second_entry(self._invariant_trace_field.gen())
            self._invariant_quaternion_algebra = (
                QuaternionAlgebraNF.QuaternionAlgebraNF(
                    self._invariant_trace_field,
                    first_entry,
                    second_entry,
                )
            )
        return self._invariant_quaternion_algebra

    def denominators(self):
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
        if self.trace_field() is None:
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
        self._denominator_residue_characteristics = (
            denominatorsforsnappy.find_prime_factors_in_a_set(norms)
        )
        return prime_ideals

    def denominator_residue_characteristics(self):
        if self._denominators is None:
            self.denominators()
        if self._denominator_residue_characteristics is None:
            prime_ideals = self._denominators
            norms = {ideal.absolute_norm() for ideal in prime_ideals}
            self._denominator_residue_characteristics = (
                denominatorsforsnappy.find_prime_factors_in_a_set(norms)
            )
        return self._denominator_residue_characteristics

    def make_prec_degree_generator(
        self,
        starting_prec=None,
        starting_degree=None,
        prec_increment=None,
        degree_increment=None,
        max_prec=None,
        max_degree=None,
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
        if starting_prec is None:
            starting_prec = self.default_starting_prec
        if starting_degree is None:
            starting_degree = self.default_starting_degree
        if prec_increment is None:
            prec_increment = self.default_prec_increment
        if degree_increment is None:
            degree_increment = self.default_degree_increment
        if max_prec is None:
            max_prec = self.default_max_prec
        if max_degree is None:
            max_degree = self.default_max_degree

        def gen():
            """
            This can be neater perhaps?
            """
            prec, degree = starting_prec, starting_degree
            yield {
                "prec": prec,
                "degree": degree,
            }
            while prec < max_prec and degree <= max_degree:
                prec = min(prec + prec_increment, max_prec)
                degree = min(degree + degree_increment, max_degree)
                yield {
                    "prec": prec,
                    "degree": degree,
                }
                if prec == max_prec and degree == max_degree:
                    break

        return gen()

    def compute_arithmetic_invariants(
        self,
        prec=None,
        degree=None,
    ):
        """
        This tries to compute the four basic arithmetic invariants: the two trace
        fields and the two quaternion algebras.

        It will also try to compute the other invariants to fill out all the attributes
        of the instance. Right now it's called upon creation of a ManifoldNT instance,
        but this can be disabled with a keyword argument when a ManifoldNT object is
        initialized.
        """
        tf_prec = (
            self.next_prec_and_degree("trace field").prec if prec is None else prec
        )
        tf_degree = (
            self.next_prec_and_degree("degree").degree if degree is None else degree
        )
        itf_prec = (
            self.next_prec_and_degree("invariant trace field").prec
            if prec is None
            else prec
        )
        itf_degree = (
            self.next_prec_and_degree("invariant trace field").degree
            if degree is None
            else degree
        )
        qa_prec = (
            self.next_prec_and_degree("quaternion algebra") if prec is None else prec
        )
        iqa_prec = (
            self.next_prec_and_degree("invariant quaternion algebra")
            if prec is None
            else prec
        )
        self.trace_field(prec=tf_prec, degree=tf_degree)
        self.invariant_trace_field(prec=itf_prec, degree=itf_degree)
        self.quaternion_algebra(prec=qa_prec)
        self.invariant_quaternion_algebra(prec=iqa_prec)
        if self._trace_field_generators:
            self.denominators()

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
                print("Quaternion algebra:")
                print(self._quaternion_algebra.ramification_string(leading_char="\t"))
            else:
                print("Quaternion algebra not found.")
        else:
            print("Trace field not found.")
        if self._invariant_trace_field:
            print("Invariant Trace field:", self._invariant_trace_field)
            print("\t Signature:", self._invariant_trace_field.signature())
            print("\t Discriminant:", self._invariant_trace_field.discriminant())
            if self._invariant_quaternion_algebra:
                print("Invariant quaternion algebra:")
                print(
                    self._invariant_quaternion_algebra.ramification_string(
                        leading_char="\t"
                    )
                )
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
        self._approx_trace_field_gens = self._snappy_mfld.trace_field_gens()
        self._approx_invariant_trace_field_gens = (
            self._snappy_mfld.invariant_trace_field_gens()
        )

    def dehn_fill(self, filling_data, which_cusp=None):
        """
        This performs dehn surgery "in place," i.e. it doesn't return a new manifold
        but rather changes self. This, of course, necessitates recomputing all the
        arithmetic invariants, so we have to override SnapPy's method to make sure
        we clean out the invariants.
        """
        self._snappy_mfld.dehn_fill(filling_data, which_cusp=which_cusp)
        self.delete_arithmetic_invariants()

    def _isomorphic_quaternion_algebras(self, other, _invariant_qa=False):
        self_field = (
            self.trace_field() if not _invariant_qa else self.invariant_trace_field()
        )
        other_field = (
            other.trace_field() if not _invariant_qa else other.invariant_trace_field()
        )
        self_qa = (
            self.quaternion_algebra()
            if not _invariant_qa
            else self.invariant_quaternion_algebra()
        )
        other_qa = (
            other.quaternion_algebra()
            if not _invariant_qa
            else other.invariant_quaternion_algebra()
        )
        if not all((self_field, other_field, self_qa, other_qa)):
            raise RuntimeError("Trace fields or quaternion algebras not known.")
        self_tf_wo_embedding = self_field.absolute_field(str(self_field.gen()))
        other_tf_wo_embedding = other_field.absolute_field(str(other_field.gen()))
        if self_tf_wo_embedding == other_tf_wo_embedding:
            return self._quaternion_algebra.is_isomorphic(other._quaternion_algebra)
        elif not self_field.is_isomorphic(other_field):
            return False
        else:
            if len(self_qa.ramified_real_places() != other_qa.ramified_real_places()):
                return False
            elif (
                self_qa.ramified_residue_characteristics()
                != other_qa.ramified_residue_characteristics()
            ):
                return False
            else:
                primitive_element = (
                    self._trace_field_numerical_root
                    if not _invariant_qa
                    else self._invariant_trace_field_numerical_root
                )
                old_anchor = (
                    other._trace_field_gens
                    if not _invariant_qa
                    else other._invariant_trace_field_gens
                )
                approx_gens = (
                    other._approx_trace_field_gens
                    if not _invariant_qa
                    else other._approx_invariant_trace_field_gens
                )
                new_anchor = [primitive_element.express(gen) for gen in approx_gens]
                special_iso = field_isomorphisms.special_isomorphism(
                    self_field,
                    other_field,
                    old_anchor,
                    new_anchor,
                )
                new_QA = other_qa.new_QA_via_field_isomorphism(special_iso)
                return self_qa.is_isomorphic(new_QA)

    def _same_denominators(self, other):
        """
        Returns whether self and other have the same non-integral primes in their trace
        fields. If the trace fields are isomorphic, it first has to force them into a
        common field.
        """
        if not all(
            (
                self.denominators(),
                other.denominators(),
                self.trace_field(),
                other.trace_field(),
            )
        ):
            raise RuntimeError("Denominators not known.")
        if (
            self._denominator_residue_characteristics
            != other._denominator_residue_characteristics
        ):
            return False
        if self._trace_field == other._trace_field:
            return self._denominators == other._denominators
        elif not self._trace_field.is_isomorphic(other._trace_field):
            return False
        else:
            primitive_element = self._trace_field_numerical_root
            old_anchor = other._trace_field_gens
            new_anchor = [
                primitive_element.express(gen) for gen in other._approx_trace_field_gens
            ]
            special_iso = field_isomorphisms.special_isomorphism(
                self._trace_field,
                other._trace_field,
                old_anchor,
                new_anchor,
            )
            new_denoms = {special_iso(ideal) for ideal in other._denominators}
            return self._denominators == new_denoms

    def compare_arithmetic_invariants(self, other):
        """
        This takes two ManifoldNTs and computes whether they have isomorphic trace
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

        Implementation detail: This function somewhat depends on what trace field
        generators are given because it distinguishes an between number fields by
        checking which one takes one set of generators to another, so it depends on the
        set of generators not e.g. containing an entire Galois orbit. In the
        future there should probably be a more robust way of making sure the generating
        set is not going to mess things up.
        """
        arith_dict = dict()
        arith_dict["trace field"] = field_isomorphisms.same_subfield_of_CC(
            self._trace_field, other._trace_field
        )
        arith_dict["invariant trace field"] = field_isomorphisms.same_subfield_of_CC(
            self._invariant_trace_field, other._invariant_trace_field
        )
        # We take the convention that quaternion algebras can be the same when the trace
        # trace fields are isomorphic even if the fields come with different embeddings.
        arith_dict["quaternion algebra"] = self._isomorphic_quaternion_algebras(other)
        arith_dict[
            "invariant quaternion algebra"
        ] = self._isomorphic_quaternion_algebras(other, _invariant_qa=True)
        arith_dict["denominators"] = self._same_denominators(other)
        return arith_dict

    def has_same_arithmetic_invariants(self, other):
        arith_dict = self.compare_arithmetic_invariants(other)
        return not (False in arith_dict.values())
