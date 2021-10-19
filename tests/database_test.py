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
