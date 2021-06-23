"""
Module for general testing of correctness and performance of this package.
"""

import ManifoldAP, database
import unittest
from sage.all import var, NumberField, CC, I
from field_isomorphisms import same_subfield_of_CC, isomorphisms_between_number_fields

def compare_against_database(filename):
    """
    Given a filename, we load a ManifoldAPDatabase using that filename, then compare its
    contents to those obtained by computing the invariants from scratch. The returned
    value is a dictionary whose keys are the manifolds and whose values booleans
    indicating whether the computed invariants are expected.
    """
    differences_dict = dict()
    with database.ManifoldAPDatabase(filename) as db:
        for dbmfld_name in db:
            dbmfld = db[dbmfld_name]
            new_mfld = ManifoldAP.ManifoldAP(str(dbmfld))
            differences_dict[str(dbmfld)] = dbmfld.compare_arithmetic_invariants(new_mfld)
    return differences_dict

def agrees_with_database(filename):
    return not (False in compare_against_database(filename).values())

def MR_test_dictionary():
    database_path = "data/"
    result = dict()
    # Census tests
    result["MR Knots"] = agrees_with_database(database_path + "MRKnots.json")
    result["MR Cusped"] = agrees_with_database(database_path + "MRKnots.json")
    result["MR Closed"] = agrees_with_database(database_path + "MRKnots.json")
    return result

def field_isomorphism_tests_as_dict():
    x = var("x")
    log_dict = dict()
    Field1 = NumberField(x**2+1, "i", embedding=I)
    Field2 = NumberField(x**2+1, "minusi", embedding=-I)
    # Checks that conjugate embeddings are considered the same.
    log_dict['Conjugate embeddings for QQ(i)'] = same_subfield_of_CC(Field1, Field2)
    # Checks that giving a nonintegral minimal polynomial doesn't mess anything up.
    # The issue would be finding the isomorphism in the first place.
    Field3 = NumberField(x**2+2*x+(5/4), "a", embedding=-1+(1/2)*I)
    log_dict['Integral and nonintegral minimal polynomials'] = same_subfield_of_CC(Field1, Field3)
    # These next two fields are the trace fields of the knots 6_1 and 7_7 respectively.
    # We test that we find an isomorphism between them and that they are distinguished
    # by their embeddings into the complex numbers.
    FieldSixOne = NumberField(x**4+x**2-x+1, "b", embedding=CC(0.547423794586058 + 0.585651979689573*I))
    FieldSevenSeven = NumberField(x**4+x**2-x+1, "c", embedding=CC(-0.547423794586059 - 1.12087348993706*I))
    log_dict['Distinguishes the trace fields of 6_1 and 7_7'] = (bool(isomorphisms_between_number_fields(FieldSixOne, FieldSevenSeven)) and not same_subfield_of_CC(FieldSixOne, FieldSevenSeven))
    return log_dict

class Figure8Tester(unittest.TestCase):
    def setUp(self):
        self.variable = var("x")
        self.mfld = ManifoldAP.ManifoldAP('4_1')
    
    def test_trace_field(self):
        self.assertTrue(self.mfld.trace_field().is_isomorphic(NumberField(self.variable**2+3, "z")))