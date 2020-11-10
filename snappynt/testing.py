"""
Module for general testing of correctness and performance of this package.
"""

import ManifoldAP

def compare_against_shelve_object(shelve_object):
    """
    Given a shelve object containing some ManifoldAP objects, we compare whether fresh
    fresh computations give the same arithmetic invariants.
    """
    success = True
    for item in shelve_object:
        name = str(item)
        mfld = ManifoldAP.ManifoldAP(name)
        if not mfld.has_same_arithmetic_invariants(item):
            success = False
            print('Differnce at', name)
    return success