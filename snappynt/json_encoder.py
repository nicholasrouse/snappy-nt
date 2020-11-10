"""
This module gives custom JSON encoding and decoding for ManifoldAP objects. The purpose
is to eventually have a test bench that can be verified for accuracy then quickly loaded
to compare future versions of this package against.

The core invariants computed in the ManifoldAP module are the two trace fields, and the
two quaternion algebras. Another invariant are the denominators, though these are often
computed alongside the (noninvariant) trace field. We want this encoding and decoding to
be somewhat flexible to allow for later information to be added. 

We should add that we don't intend for the JSON representations to be totally minimal:
since JSON is human readable, we include some extra information for quick reference.
"""

import json
import ManifoldAP
import QuaternionAlgebraNF
from sage.rings.number_field.number_field import is_NumberField
from sage.all import CC, radical


class FieldEncoder(json.JSONEncoder):
    """
    Returns a JSON dict for fields as they come to us from Kleinian groups. This will
    intentionally raise an exception if there's no specified embedding when gen_embedding
    gets called.
    """

    def default(self, field):
        if is_NumberField(field):
            if field.gen_embedding() is None:
                raise AttributeError("There is no specified embedding")
            d = {
                "defining polynomial": str(field.defining_polynomial()),
                "numerical root": {
                    "real part": float(CC(field.gen_embedding()).real_part()),
                    "imaginary part": float(CC(field.gen_embedding()).imag_part()),
                },
                "discriminant": int(field.discriminant()),
                "signature": str(field.signature()),
            }
            return d
        else:
            return json.JSONEncoder(self, field)


class QuaternionAlgebraEncoder(json.JSONEncoder):
    """
    Returns a JSON dictionary for Quaternion Algebras over number fields. This one will
    be a little narrower than what one would expect for a totally general such object.
    For example, we don't repeat the base field because this will be recorded as the
    trace field (or invariant trace field), but the ground field is of course totally
    indispensable information for a quaternion algebra in general. Since computing
    ramification sets can actually be computationally expensive in general, we have a
    special field to record the dyadic ramification (which tends to be the hard ones).
    """

    def default(self, algebra):
        if isinstance(algebra, QuaternionAlgebraNF.QuaternionAlgebraNF):
            a, b = algebra.invariants()
            d = {
                "Hilbert symbol": {"first entry": str(a), "second entry": str(b),},
                "ramification": {
                    "dyadic ramification": [
                        str(element.gens())
                        for element in sorted(
                            [
                                ideal
                                for ideal in algebra.ramified_finite_places()
                                if radical(ideal.absolute_norm()) == 2
                            ],
                            key=lambda ideal: ideal.residue_class_degree(),
                        )
                    ],
                    "nondyadic ramification": [
                        str(element.gens())
                        for element in sorted(
                            [
                                ideal
                                for ideal in algebra.ramified_finite_places()
                                if radical(ideal.absolute_norm()) != 2
                            ],
                            key=lambda x: radical(x.absolute_norm()),
                        )
                    ],
                    "residue characteristics": [int(element) for element in algebra.ramified_residue_characteristics().elements()],
                    "real ramification": [
                        {
                            "real part": float(CC(place.im_gens()[0]).real_part()),
                            "imaginary part": float(CC(place.im_gens()[0]).imag_part()),
                        }
                        for place in algebra.ramified_real_places()
                    ],
                },
            }
            return d
        else:
            return json.JSONEncoder(self, algebra)


class ManifoldAP_Encoder(json.JSONEncoder):
    def default(self, mfld):
        d = {
            "name": str(mfld),
            "volume": float(mfld.volume()),
            "trace field": FieldEncoder().default(mfld.trace_field()),
            "quaternion algebra": QuaternionAlgebraEncoder().default(
                mfld.quaternion_algebra()
            ),
            "invariant trace field": FieldEncoder().default(
                mfld.invariant_trace_field()
            ),
            "invariant quaternion algebra": QuaternionAlgebraEncoder().default(
                mfld.invariant_quaternion_algebra()
            ),
            "denominators": {
                "ideals": [
                    str(element)
                    for element in sorted(
                        [ideal.gens() for ideal in mfld.denominators()],
                        key=lambda ideal: radical(ideal.absolute_norm()),
                    )
                ],
                "residue characteristics": sorted(
                    [
                        int(residue_char)
                        for residue_char in mfld.denominator_residue_characteristics()
                    ]
                ),
            },
            "arithmetic" : mfld.is_arithmetic()
        }
        return d