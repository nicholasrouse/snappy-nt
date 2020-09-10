"""
This module contains a class for arbitrary precision holonomy groups. It doesn't do
anything that can't be done working with SnapPy's builtin objects, but I think it will
ultimately make our life easier, particularly in making the groups play nicer with Sage.

It's important to note that we can't actually subclass the polished holonomy class
because that one has a particular precision specified at creation.

Last updated: Sept-8 2020
"""