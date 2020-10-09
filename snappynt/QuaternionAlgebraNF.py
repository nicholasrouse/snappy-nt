"""
A module for quaternion algebras over number fields. The class will inherit from Sage's
quaternion algebras. The key method will be an is_isomorphic method which leverages the
Albert-Brauer-Hasse-Noether theorem to check for isomorphism. The primary application
will be for quaternion algebras associted to Kleinian groups, so there might be a few
other methods which are attuned to this purpose.
"""
from sage.all import QuaternionAlgebra
from sage.algebras.quatalg.quaternion_algebra import QuaternionAlgebra_ab as QuaternionAlgebra
import field_isomorphisms

class QuaternionAlgebraNF(QuaternionAlgebra):
    """
    __init__ method should be the same as Sage's so no need to overload it
    """
    def is_isomorphic(self, other):
        """
        Given two quaternion algebras over number fields, this function tests for
        isomorphism. The first check is that their base fields are isomorphic. Assuming
        that they are, the function uses the Albert-Brauer-Hasse-Noether theorem for
        quaternion algebras. The net effect is to check whether the ramification sets
        are the same. This can be a bit delicate in total generality, e.g. if self and
        other have isomorphic base fields but they are given by different minimal
        polynomials.

        Fortunately for us, given an isomorphism of fields f: K -> L and an ideal I of
        the ring of integers of K, Sage will correctly interpret f(I) as an ideal of 
        (the ring of integers) of L.
        """
