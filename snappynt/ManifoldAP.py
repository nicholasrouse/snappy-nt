"""
This is the module that contains the class for arbitrary precision Manifolds. As of
Aug-23 2020, the only objects that are actually computed to arbitrary precision are
reps into SL_2(CC) and the arithmetic invariants (e.g. trace fields and quaternion
algebras).

Things to consider:
    1. Building some kind of database that can be loaded to avoid repeating expensive
    computations. I haven't decided how exactly to build and access such a database
    yet.

    2. Perhaps importing functions that actually compute the various invariants from
    numerical input. I.e. make another module and put all the ugly implementation for
    computations there.

"""


import snappy, denominatorsforsnappy
from sage.all import factor, NumberField, QuaterionAlgebra
import math
import functools



class ManifoldAP(snappy.Manifold):
    # Probably make these changeable via a class method at some point.
    # Not sure if this a great pattern to have these as class level attributes.
    # It's possible it should be managed on a per instance level.
    default_starting_prec = 1000
    default_starting_degree = 10
    default_max_prec = 10 ** 6
    default_max_degree = 100
    default_prec_increment = 5000
    default_degree_increment = 5

    def __init__(self, spec=None):
        snappy.Manifold.__init__(self, spec)
        # We store the fields as a sage NumberField objects with a *monic* generating polynomial.
        # Perhaps subclass Sage's NumberField to store all this info?
        # The prec_record variables record whether a given prec and degree wsa sucessful.
        self.trace_field = None  
        self.trace_field_numerical_root = None 
        self.trace_field_generators = None
        self.trace_field_prec_record = dict()
        self.invariant_trace_field = None
        self.invariant_trace_field_numerical_root = None
        self.invariant_trace_field_generators = None
        self.invariant_trace_field_prec_record = dict()
        self.quaternion_algebra = None
    
    def has_two_torsion_in_homology(self):
        """
        Returns True if there is two-torsion in homology and False if not. This doesn't
        really need arbitrary precision, but it hopefully makes some other code
        cleaner.
        """
        homology = self.homology()
        elementary_divisors = homology.elementary_divisors()
        for divisor in elementary_divisors:
            if divisor%2 == 0:
                return True
        return False

    def defining_function(self, prec):
        return snappy.snap.polished_holonomy(bits_prec=prec)

    def compute_trace_field_fixed_prec(
        self, prec=default_starting_prec, degree=default_starting_degree
    ):
        approx_trace_field = snappy.snap.trace_field_gens(self)
        exact_field_data = approx_trace_field.find_field(
            prec=prec, degree=degree, optimize=True
        )
            # This will override previous calculations with same prec and degree.
            # It's unclear if we want this behavior.
            self.trace_field_prec_record[(prec, degree)] = bool(exact_field)
        if exact_field_data is not None:
            self.trace_field = exact_field_data[0]
            self.trace_field_numerical_root = exact_field_data[1] 
            self.trace_field_generators = exact_field_data[2] 
        return self.trace_field

    def compute_trace_field(
        self,
        starting_prec=default_starting_prec,
        starting_degree=default_starting_degree,
        prec_increment=default_prec_increment,
        degree_increment=default_degree_increment,
        max_prec=default_max_prec,
        max_degree=default_max_degree,
        verbosity=False,
        use_last_known_failed=False
    ):
        """
        This is the exact field, returned as a sage NumberField. The exact generators for the
        field are not returned by this method to allow for easier interface with sage. They are
        however computed and stored as an  attribute for later use. The method has a
        semicomplicated interface to allow for multiple attempts to find the field. If only one
        attempt is required (e.g. because the requisite precision and degree are known), the
        method compute_trace_field_fixed_prec is probably better and will store the result if
        successful.

        I think perhaps I should put all the code for computing this in another module
        and just import for this one. I do need to decide on an interface for this one
        though.
         
        Aug-1-2020
        """
        if self.trace_field:
            return self.trace_field
        elif self.invariant_trace_field and not self.has_two_torsion_in_homology():
            return self.invariant_trace_field
        else:
            prec = starting_prec
            degree = starting_degree
            exact_field = None
            while exact_field == None:
                if verbosity:
                    print(
                        str(manifold) + ":",
                        f"Trying with precision={prec} and degree={degree}",
                    )
                exact_field = ManifoldAP.compute_trace_field_fixed_prec(
                    self, prec=prec, degree=degree
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
    
    def compute_quaternion_algebra_fixed_prec(self, prec=default_starting_prec, degree=default_starting_degree):
        """
        If the trace field isn't known, whatever precision and degree are passed here
        are used to try to compute it. If it fails to do so, the entire function fails,
        and will return None.
        """
        if not self.trace_field:
            self.compute_trace_field_fixed_prec(prec=prec, degree=degree)
        if not self.trace_field: return None
        
        
        
