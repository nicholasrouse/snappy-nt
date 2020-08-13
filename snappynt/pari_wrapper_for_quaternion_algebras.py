"""
This file is to embedd a quaternion algebra coming from  PARI.
Doing so is actually pretty necessary because sage's functionality (as of 9.1)
seems to be fairly limited for quaternion algebras defined over number fields other than QQ.
However, PARI's tools are much more extensive.
The big one is computing ramified places. One can more or less compute ramified places
in sage using tools from NumberField, but it does so without actually creating a
quaternion algebra. Using a wrapper for PARI allows us to use all of PARI's functions
for quaternion algebras inside sage.
"""
import cypari2  # This is necessary for some type testing.
from sage.libs.pari.convert_sage import gen_to_sage
from sage.all import (NumberField, )  # Can I do this somewhere else?


def is_probably_pari_nf(nf):
    """
    This tests somewhat crudely inside sage whether the object looks
    like the ouput of a PARI nfinit. Note that this doesn't work for bnfinit etc.
    """
    if not isinstance(
        nf, cypari2.gen.Gen
    ):  # All pari objects in sage have this python type.
        return False
    if len(nf) != 9:  # Maybe come back and try to catch TypeErrors for this one.
        return False
    if nf.type() != "t_VEC":
        return False
    if nf[0].type() != "t_POL":
        return False
    return True


def convert_parinf_to_sage(nf, variable="z"):
    """
    Converts the output of a pari_nf() call on a sage number field back to a sage number field.
    It sounds dumb when I say it like that, but it allows us to have functions that can easily accept both sage
    and PARI number fields.

    Note that pari_poly will always be monic. So if the original sage number field had a nonmonic generator,
    the output of this function will have a different minimal polynomial.
    """
    x = var("x")
    if not is_probably_pari_nf(nf):
        raise TypeError("This does not appear to be the output of a pari_nf() call")
    pari_poly = gen_to_sage(
        nf[0], {"y": x}
    )  # This assumes the default 'y' variable for pari_nf().
    sage_number_field = NumberField(
        pari_poly, variable
    )  # Doesn't actually define the variable outside this function.
    return sage_number_field


class PARIQuaternionAlgebraOverNumberField:
    """
    This class wraps a PARI quaternion algebra over a number field.
    The idea is to make a sage accessible object so that conversions to and from PARI
    are taken care of.
    """

    def __init__(self, number_field, a, b, variable="z"):
        """
        As of Jul-15-2020, this isn't as robust as I would like to be, but the current goal
        is to get it working with SnapPy. Eventually I would like it to be a completely general
        tool for working with quaternion algebras over number fields inside sage.

        One potential obstacle is that PARI's number fields always have monic defining polynomial.
        The problem is that we need to make sure we express the Hilbert symbol entries in terms
        of PARI's generator since it actually does all the computations. We use the module 

        This initalizes a quaternion algebra described by the Hilbert symbol (a,b)
        over the number field number_field.
        Ideally this will be robust enough to accept number fields and elements that are either the correct type in sage
        or PARI objects inside sage or even polynomials to create the associated number fields.
        Let us restrain ourselves for now to avoid creating field extensions out of arbitrary
        symbolic expressions like ``sqrt(2)".

        In usage I expect the number_field argument to be passed in as a sage number field;
        however, in building the functionality for hyperbolic 3-manifolds, it will be useful to feed PARI
        number fields in with converting them to sage number fields.

        It's also worth noting that sage number fields always carry the pari nf object with them somehow.
        Therefore, it might somehow be more memory efficient to rewrite this to not save the pari number field.
        On the other hand, I'm not sure how expensive the function calls to pari_nf() are.

        We don't systematically keep track of the embeddings of the number field.
        Aug-12 2020
        """
        if isinstance(
            number_field, sage.rings.number_field.number_field.NumberField_absolute
        ):
            self.sage_number_field = number_field
            self.pari_number_field = (
                self.number_field.pari_nf()
            )  # This makes PARI variable "y". We can change this if necessary.
        elif isinstance(
            number_field, sage.rings.number_field.number_field.NumberField_relative
        ):  # Relative number fields are converted to absolute ones.
            self.sage_number_field = number_field.absolute_field(
                names=variable
            )  # E.g. if variable='z' (the default). This makes a sage number field with the variable z.
            self.pari_number_field = self.sage_number_field.pari_nf()
        elif is_probably_pari_nf(number_field):
            self.pari_number_field = number_field
            self.sage_number_field = convert_parinf_to_sage(number_field)
        if a in self.sage_number_field and b in self.sage_number_field:
            self.a = a
            self.b = b
        else:
            raise TypeError(
                "The Hilbert symbol entries are not recognized by Sage as belonging to the number field."
            )
