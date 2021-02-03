"""
Module for making and managing databases of ManifoldAP objects. These databases are
fairly rudimentary at this point. Their primary purposes are to have a test bench of
known correct computations to test revisions of the package against and to allow saving
of calculations in a persistent and human readable way. It's worth pointing out that
basically all of this is just translating between some JSON encoding and a shelve
object, but it at least allows building the test bench from a human readable source.

It might be important to note that one shouldn't load untrusted binary data using this
module for the same reasons one shouldn't do so with shelve or pickle.

The JSON encoded object should be a JSON array the entries of which are encoded JSON
encoded ManifoldAP objects.
"""

import shelve
import json_encoder

class ManifoldAPDatabase:
    """
    This class wraps both a quasi-JSON human readable file object and a shelve object,
    and keeps them consistent.
    """
    def __init__(filename):
        """
        The extant object could be either a shelve object or a file with an array of
        encoded ManifoldAP objects. The file can also not yet exist, in which case one
        will be created with the name.
        """
        with shelve.open(filename) as shelve_object:
            pass