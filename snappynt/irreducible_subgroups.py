from snappy import polished_holonomy

def generate_reducible_subgroup(g, h):
    """
    Tests whether the subgroup of  generated by g and h is irreducible. The function
    actually tests whether the trace of the commutator of g and h is less than epsilon
    from 2 where epsilon is the smallest nonzero quantity that the base field with
    specified precision can recognize.  This has the effect that the function might
    incorrectly say that two elements generate a reducible subgroup when their
    commutator has trace closer (but not equal) to 2 that the ambient precision.
    On the other hand, it will never certify that a subgroup is irreducible when it is
    not.
    Aug-13 2020
    """
    # This makes things a bit more robust, but may also slow things down unnecessarily.
    base_field_g = g.parent().base()
    base_field_h = h.parent().base()
    if base_field_g.is_exact() and base_field_h.is_exact(): epsilon = 0
    else: epsilon = max(base_field_g.epsilon(), base_field_h.epsilon())
    trace = g.commutator(b).trace()
    if abs(trace-2) <= epsilon: return True
    else: return False