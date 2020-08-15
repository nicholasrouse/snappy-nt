import snappy, denominatorsforsnappy
from sage.all import factor, NumberField
import math
import functools


class ManifoldAP(snappy.Manifold):
    default_starting_prec = 1000
    default_starting_degree = 10
    default_max_prec = 10 ** 5
    default_max_degree = 100
    default_prec_increment = 2000
    default_degree_increment = 5

    def __init__(self):
        snappy.Manifold.__init__(self)
        # We store the fields as a sage NumberField objects with a *monic* generating polynomial.
        self.trace_field = None  
        self.trace_field_numerical_root = None 
        self.trace_field_generators = None
        self.invariant_trace_field = None
        self.invariant_trace_field_generators = None
        self.quaternion_algebra = None
        self.invariant_trace_field = None

    def defining_function(self, prec):
        return snappy.snap.polished_holonomy(bits_prec=prec)

    def compute_trace_field_fixed_prec(
        self, prec=default_starting_prec, degree=default_starting_degree
    ):
        approx_trace_field = snappy.snap.trace_field_gens(self)
        exact_field_data = approx_trace_field.find_field(
            prec=prec, degree=degree, optimize=True
        )
        if exact_field is not None:
            exact_field[0] = self.trace_field
            exact_field[1] = self.trace_field_numerical_root
            exact_field[2] = self.trace_field_generators
        return exact_field

    def compute_trace_field(
        self,
        starting_prec=ManifoldAP.default_starting_prec,
        starting_degree=ManifoldAP.default_starting_degree,
        prec_increment=default_prec_increment,
        degree_increment=default_degree_increment,
        max_prec=None,
        max_degree=None,
        verbosity=False,
    ):
        """
        This is the exact field, returned as a sage NumberField.
        The exact generators for the field are not returned by this method
        to allow for easier interface with sage. They are however computed and
        stored as a class attribute for later use.
        The method has a semicomplicated interface to allow for 
        multiple attempts to find the field.
        If only one attempt is required (e.g. because the requisite precision and degree are known),
        the method compute_trace_field_fixed_prec is probably better and will store the result if successful.
        Aug-1-2020
        """
        if self.trace_field:
            return self.trace_field
        else:
            prec = starting_prec
            degree = starting_degree
            approx_trace_field = snappy.snap.trace(self)
            exact_field = None
            while exact_field == None:
                if verbosity:
                    print(
                        str(manifold) + ":",
                        f"Trying with precision={prec} and degree={degree}",
                    )
                exact_field = compute_trace_field_fixed_prec(
                    manifold, prec=prec, degree=degree
                )
                if prec == max_prec and degree == max_degree:
                    return None
                if prec + prec_increment <= max_prec:
                    prec = prec + prec_increment
                else:
                    prec = max_prec
                if degree + degree_increment <= max_degree:
                    degree = degree + degree_increment
                else:
                    degree = max_degree
            return self.trace_field
