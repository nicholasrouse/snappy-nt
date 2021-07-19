from paths import convert_rel_to_abs
from sage.all import NumberField, var

import snappynt.database as database
import snappynt.ManifoldNT as ManifoldNT


def compare_against_database(filename):
    """
    Given a filename, we load a ManifoldNTDatabase using that filename, then compare its
    contents to those obtained by computing the invariants from scratch. The returned
    value is a dictionary whose keys are the manifolds and whose values booleans
    indicating whether the computed invariants are expected.
    """
    differences_dict = dict()
    with database.ManifoldNTDatabase(filename) as db:
        for dbmfld_name in db:
            dbmfld = db[dbmfld_name]
            new_mfld = ManifoldNT.ManifoldNT(str(dbmfld))
            for _ in range(4):
                if not new_mfld._arithmetic_invariants_known():
                    new_mfld.compute_arithmetic_invariants()
            differences_dict[str(dbmfld)] = dbmfld.has_same_arithmetic_invariants(
                new_mfld
            )
    return differences_dict


def agrees_with_database(filename):
    return not (False in compare_against_database(filename).values())


def test_fig8_field():
    mfld = ManifoldNT.ManifoldNT("4_1")
    tf = mfld.trace_field()
    x = var(tf.defining_polynomial().variable_name())
    field = NumberField(x ** 2 + 3, "z")
    assert field.is_isomorphic(tf)


def test_knots():
    assert agrees_with_database(convert_rel_to_abs("data/MRKnots.json"))


def test_cusped():
    assert agrees_with_database(convert_rel_to_abs("data/MRCusped.json"))


def test_closed():
    assert agrees_with_database(convert_rel_to_abs("data/MRClosed.json"))
