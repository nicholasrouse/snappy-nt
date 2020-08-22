"""
This is meant to be imported into sage with snappy already imported.
The only goal is to make an iterator of large knots using the
largeknotsforsnappy text file.

Ok, now (Jul 29 2020) we're expanding the module.

"""
from snappy import Manifold

def large_knots():
    with open('largeknotsforsnappy') as large_knots:
        for knot in large_knots:
            yield Manifold(knot)