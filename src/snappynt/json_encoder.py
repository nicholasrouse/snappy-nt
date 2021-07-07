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
from collections import Counter

from sage.all import CC, QQ, ZZ, NumberField, PolynomialRing, RealField, radical
from sage.rings.number_field.number_field import is_NumberField

from . import ManifoldAP, QuaternionAlgebraNF


def nested_decoder(func):
    """
    A decorator that tries to fix decoding of special objects inside containers. For
    example, if I have a JSON array whose elements are encoded ManifoldAP objects, I'd
    like to be able to call json.loads with cls=ManifoldAP_Decoder. As written without
    this decorator, there will be an error because ManifoldAP_Decoder is coded to expect
    a dictionary from the default JSON decoder, which is then pieces back into a
    ManifoldAP object. There is maybe a better way to code my decoders that makes
    this decorator unnecessary.

    The basic logic is to try to decode the object, if it's a nested object, then there
    will be some error thrown. In particular, all the encoding in this module converts
    the objects to a dictionary all of whose keys are strings that is serializable. So
    if the object fed in is actually a list, we'll get a TypeError from something like
    L["name"] where L is a list. Similarly, if the object we feed in is a dictionary but
    not one of the form that gets encoded by JSON, then we'll need to catch a KeyError.
    """

    def wrapper(instance, text):
        """
        The instance should be an instance of the decoder class. We need it for the
        call signature to be correct.
        """
        try:
            return func(instance, text)
        except (TypeError, AttributeError):
            # Think that text actually gets decoded to a Python list.
            # What we need to do is load the text, which return a Python list of Python
            # dictionaries. We can't just call func on these dictionaries because func
            # expects a string, but they will be dumpable via the default JSON encoder.
            # So after we dump them, we get strings that func can handle.
            raw_list = json.JSONDecoder().decode(text)
            list_of_encoded_text = [json.dumps(elt) for elt in raw_list]
            list_of_decoded_objects = [
                func(instance, elt) for elt in list_of_encoded_text
            ]
            return list_of_decoded_objects
        except KeyError:
            # Think that text is a dictionary whose values are the dicts that func is a
            # special decoder for.
            raw_dict = json.JSONDecoder().decode(text)
            dict_of_encoded_text = {key: json.dumps(raw_dict[key]) for key in raw_dict}
            dict_of_decoded_objects = {
                key: func(instance, dict_of_encoded_text[key])
                for key in dict_of_encoded_text
            }
            return dict_of_decoded_objects

    return wrapper


class FieldEncoder(json.JSONEncoder):
    """
    Returns a JSON dict for fields as they come to us from Kleinian groups. This will
    intentionally raise an exception if there's no specified embedding when gen_embedding
    gets called. This has some shortcomings for use with totally general number field.
    In particular, the fields encoded should be generated by a single element. Of course
    every number field is generated by a primitive element, but they're not always
    presented that way. Also, relative extensions will probably break things. However,
    as far as I'm aware, there's no reason trace fields should be relative extensions.
    """

    def default(self, field):
        if is_NumberField(field):
            if field.gen_embedding() is None:
                raise AttributeError("There is no specified embedding")
            d = {
                "defining polynomial": str(field.defining_polynomial()),
                "generator name": str(field.gen()),
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


def string_to_poly(s, variable="x"):
    """
    Takes a string representation s of a polynomial and makes it a polynomial that sage
    can recognize. The polynomial should have rational coefficients. If the coefficients
    are (rational) integers, then the returned polynomial will be recognized by sage as
    an polynomial in ZZ[x]. The variable keyword argument needs to match the variable
    that appears in the string representation. This will sometimes work with
    """
    x = variable
    try:
        ring = PolynomialRing(ZZ, x)
        poly = ring(s)
    except TypeError:
        ring = PolynomialRing(QQ, x)
        poly = ring(s)
    return poly


def fix_ideal_string(s):
    """
    Takes a string representing an ideal in a number field and returns a list of strings
    that sage can parse correctly. Sometimes sage won't interpret the strings correctly.
    For example, in K=Q(sqrt2), we have
    sage: K.ideal("(2,sqrt2)")
    Fractional ideal (4)

    But of course we really want the ideal generated by sqrt2. For that example, this
    function returns ["2", "sqrt2"], which can be parsed by sage correctly:
    sage:K.ideal(["sqrt2", "2"])
    Fractional ideal (sqrt2)
    """
    if s[0] == "(":
        s = s[1:]
    if s[-1] == ")":
        s = s[0:-1]
    return s.split(",")


def dict_to_field(d):
    """
    Takes a Python/Sage dictionary to a number field object. The keys should be
    "defining polynomial"
    "generator name"
    "numerical root", which has
        "real part"
        "imaginary part"
    There can be other keys, but these must be present. The value for "defining
    polynomial" can possibly be a sage polynomial, but this is an "off label" use case.
    The point of having this be a separate function rather than just in the FieldDecoder
    class is to better handle nested JSON objects. E.g. ManifoldAP objects have fields
    associated to them, but when we decode ManifoldAP database objects, the field
    information becomes a Python dictionary rather than encoded JSON.
    """
    d = {key.lower(): d[key] for key in d}
    defining_poly = string_to_poly(d["defining polynomial"])
    real = d["numerical root"]["real part"]
    imag = d["numerical root"]["imaginary part"]
    numerical_root = CC(real, imag)
    field = NumberField(defining_poly, d["generator name"], embedding=numerical_root)
    return field


class FieldDecoder(json.JSONDecoder):
    @nested_decoder
    def decode(self, text):
        raw_dict = json.JSONDecoder().decode(text)
        # defining_poly = string_to_poly(raw_dict["defining polynomial"])
        # real= raw_dict["numerical root"]["real part"]
        # imag = raw_dict["numerical root"]["imaginary part"]
        # numerical_root = CC(real, imag)
        # field = NumberField(defining_poly, raw_dict["generator name"], embedding=numerical_root)
        return dict_to_field(raw_dict)


class QuaternionAlgebraEncoder(json.JSONEncoder):
    """
    Returns a JSON dictionary for Quaternion Algebras over number fields. Since computing
    ramification sets can actually be computationally expensive in general, we have a
    special field to record the dyadic ramification (which tends to be the hard ones).
    """

    def default(self, algebra):
        if isinstance(algebra, QuaternionAlgebraNF.QuaternionAlgebraNF):
            a, b = algebra.invariants()
            d = {
                "base field": FieldEncoder().default(algebra.base_ring()),
                "Hilbert symbol": {
                    "first entry": str(a),
                    "second entry": str(b),
                },
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
                    "residue characteristics": {
                        str(element): int(
                            algebra._ramified_residue_characteristics[element]
                        )
                        for element in algebra._ramified_residue_characteristics
                    },
                    "real ramification": [
                        str(place.im_gens()[0])
                        for place in algebra.ramified_real_places()
                    ],
                },
            }
            return d
        else:
            return json.JSONEncoder(self, algebra)


def real_place_finder(field, s):
    """
    For a field field and a string s representing a real place, finds the real place of
    the field closest to the string. The string should be a string that looks like a
    real number.
    """
    number = RealField(106)(s)
    places = field.real_places()
    return min(places, key=lambda place: abs(number - place.im_gens()[0]))


def dict_to_quaternion_algebra(d):
    """
    Converts a dictionary to a QuaternionAlgebraNF object. The main use case is decoding
    serialized JSON objects. However, it's good to have it a separate function for
    nested decodings. The dictionary needs to conform to the following interface. It
    should have keys and subkeys
    "base field"
        (see dict_to_field function for interface)
    "ramification"
        "real ramification"
        "nondyadic ramification"
        "dyadic ramification"
        "residue characteristics"
    "Hilbert symbol"
        "first entry"
        "second entry"
    The interface looks a bit complicated, but the primary usecase is decoding encoded
    quaternion algebras. If they are encoded using the QuaternionAlgebraEncoder class
    and then decoded using the default JSON decoder, the original quaternion algebra
    will be recovered. E.g.
    F = NumberField(x^2-2,"z", embedding=1.4)
    QA = QuaternionAlgebraNF.QuaternionAlgebraNF(F, F.gen(), F.gen() - 5)
    s = json.dumps(QA, cls=json_encoder.QuaternionAlgebraEncoder)
    d = json.loads(s)
    QAPrime = dict_to_quaternion_algebra(d)
    QA.is_isomorphic(QAPrime) # Should return True.
    """
    ramification_dict = d["ramification"]
    field = dict_to_field(d["base field"])
    a = field(d["Hilbert symbol"]["first entry"])
    b = field(d["Hilbert symbol"]["second entry"])
    algebra = QuaternionAlgebraNF.QuaternionAlgebraNF(
        field, a, b, compute_ramification=False
    )
    ramified_real_places = {
        real_place_finder(field, s) for s in ramification_dict["real ramification"]
    }
    ramified_nondyadic_places = {
        field.ideal(fix_ideal_string(s))
        for s in ramification_dict["nondyadic ramification"]
    }
    ramified_dyadic_places = {
        field.ideal(fix_ideal_string(s))
        for s in ramification_dict["dyadic ramification"]
    }
    algebra._ramified_real_places = ramified_real_places
    algebra._ramified_finite_places = ramified_nondyadic_places | ramified_dyadic_places
    algebra._ramified_residue_characteristics = Counter(
        {
            int(key): ramification_dict["residue characteristics"][key]
            for key in ramification_dict["residue characteristics"]
        }
    )
    algebra._ramified_dyadic_places_computed = True
    return algebra


class QuaternionAlgebraDecoder(json.JSONDecoder):
    """
    Decoding a JSON encoded QuaternionAlgebraNF. See the dict_to_quaternion_algebra
    function for details.
    """

    @nested_decoder
    def decode(self, text):
        text_dict = json.JSONDecoder().decode(text)
        algebra = dict_to_quaternion_algebra(text_dict)
        return algebra


class ManifoldAP_Encoder(json.JSONEncoder):
    def default(self, mfld):
        d = {
            "name": str(mfld),
            "volume": float(mfld.volume()),
            "quaternion algebra": QuaternionAlgebraEncoder().default(
                mfld.quaternion_algebra()
            ),
            "invariant quaternion algebra": QuaternionAlgebraEncoder().default(
                mfld.invariant_quaternion_algebra()
            ),
            "denominators": {
                "ideals": [
                    str(element.gens())
                    for element in sorted(
                        list(mfld.denominators()),
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
            "arithmetic": mfld.is_arithmetic(),
        }
        return d


def dict_to_manifold(d):
    """
    For similar reason as the other dict_to functions, we have a standalone function
    that converts the Python dictionary returned by the default JSON encoding to a
    ManifoldAP object. The necessary interface is
    "name"
    "quaternion algebra"
    "invariant quaternion algebra"
    "denominators"
    It's worth pointing out that the trace fields are part of the quaternion algebras.
    Also, since computing volume is cheap, we just recompute it when we recreate the
    manifold so that we don't lose precision. We also recompute the residue
    characteristics, but this may or may not be faster. The idea is we don't want to
    break things if one day we decide that the residue characteristics don't belong in
    the JSON object.
    """
    mfld = ManifoldAP.ManifoldAP(d["name"], delay_computations=True)
    quaternion_algebra = dict_to_quaternion_algebra(d["quaternion algebra"])
    invariant_quaternion_algebra = dict_to_quaternion_algebra(
        d["invariant quaternion algebra"]
    )
    trace_field = quaternion_algebra.base_ring()
    invariant_trace_field = invariant_quaternion_algebra.base_ring()
    denominators = {
        trace_field.ideal(fix_ideal_string(s)) for s in d["denominators"]["ideals"]
    }
    mfld._trace_field, mfld._invariant_trace_field = trace_field, invariant_trace_field
    mfld._quaternion_algebra, mfld._invariant_quaternion_algebra = (
        quaternion_algebra,
        invariant_quaternion_algebra,
    )
    mfld._denominators = denominators
    mfld._denominator_residue_characteristics = (
        mfld.denominator_residue_characteristics()
    )
    return mfld


class ManifoldAP_Decoder(json.JSONDecoder):
    """
    Converts a serialized JSON object to a ManifoldAP object.
    """

    @nested_decoder
    def decode(self, text):
        decoded_text = json.JSONDecoder().decode(text)
        mfld = dict_to_manifold(decoded_text)
        return mfld


def encode_list_of_manifolds(list_of_manifolds):
    """
    Given a list of manifolds (or some other iterable with a similar interface), returns
    a Python list of encoded manifolds. The returned object should be serializable with
    the default JSON encoder. This ends up being unnecessary and unused.
    """
    return [json.dumps(mfld, cls=ManifoldAP_Encoder) for mfld in list_of_manifolds]


def decode_list_of_manifolds(list_of_manifolds):
    """
    Given a Python list whose elements are encoded ManifoldAP objects, returns a list
    with the decoded ManifoldAP objects.
    """
    return [json.loads(mfld, cls=ManifoldAP_Decoder) for mfld in list_of_manifolds]


class ManifoldAP_List_Encoder(json.JSONEncoder):
    """
    I think this class is basically unneccesary as if we pass a list of ManifoldAPs in
    with the class stipulated as ManifoldAP_Encoder, the list will be encoded correctly.
    However, the dual decode function will not work because of how the ManifoldAP
    encoder is written (which suggests I should rewrite it at some point). I leave this
    for now so that we have ManifoldAP_List_Encoder and ManifoldAP_List_Decoder
    """

    def default(self, list_of_manifolds):
        return [ManifoldAP_Encoder().default(mfld) for mfld in list_of_manifolds]


class ManifoldAP_List_Decoder(json.JSONDecoder):
    def decode(self, text):
        text_list = json.JSONDecodeError().decode(text)
        return decode_list_of_manifolds(text_list)