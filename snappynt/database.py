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
import shelve, dbm, os.path, time
import json_encoder, json
import collections.abc

def json_array_to_dict(s):
    """
    Takes a Python string that represents an encoded JSON array with elements ManifoldAP
    objects and converts it to a Python dictionary where the keys are the names of the
    manifolds and the values are the ManifoldAP objects. Note that s cannot be a file,
    but we provide another function json_file_to_dict to handle this. It's also
    important that a single encoded ManifoldAP string is not passed in.
    """
    list_of_manifolds = json.loads(s, cls=json_encoder.ManifoldAP_Decoder)
    dict_of_manifolds = {str(mfld) : mfld for mfld in list_of_manifolds}
    return dict_of_manifolds

def json_file_to_dict(fp):
    """
    Takes a file object containing JSON text representing an array of ManifoldAP objects, and
    return a Python dictionary whose keys are the names of the manifolds and whose
    values are the ManifoldAP objects.
    """
    list_of_manifolds = json.load(fp, cls=json_encoder.ManifoldAP_Decoder)
    dict_of_manifolds = {str(mfld) : mfld for mfld in list_of_manifolds}
    return dict_of_manifolds

def looks_like_a_json_file(filename, cls=json_encoder.ManifoldAP_Decoder):
    """
    Takes a filename and tests whether it can load it using the
    ManifoldAP_Decoder (or whatever cls is passed in).
    """
    try:
        with open(filename, 'r') as fp:
            list_of_manifolds = json.load(fp, cls=cls)
            for mfld in list_of_manifolds:
                mfld.trace_field()
            return True
    except (FileNotFoundError, UnicodeDecodeError, json.JSONDecodeError):
        return False

def looks_like_a_shelve_file(filename):
    """
    Takes a filename and tests whether the corresponding file contains a shelve object.
    Basically it just tries to open it and do some rudimentary operations. It catches
    the dbm.error exception in the case that it doesn't look like a shelve object.
    It's important to note that if the file doesn't exist, this will not return true.
    """
    try:
        with shelve.open(filename, flag='r') as shelve_object:
            temp_dict = dict()
            for elt in shelve_object:
                temp_dict[elt] = shelve_object[elt]
                break
        return True
    except dbm.error[0]:
        return False

def change_file_extension(filename, old_extension, new_extension):
    """
    This function rips off old_extension (if it's present) and appends new_extension.
    It tries to be smart enough to handle whether the initial period is passed in. So,
    old_extenion="json" and old_extension=".json" should both work. If new_extension is
    already present at the end of the filename, then this function returns the original
    filename. If neither old_extension nor new_extension ends filename, then
    old_extension is appended without ripping anything off.
    """
    old_extension = old_extension[1:] if old_extension[0] == "." else old_extension
    new_extension = new_extension[1:] if new_extension[0] == "." else new_extension
    if filename[-len(old_extension):] == old_extension:
        new_name = filename[:-len(old_extension)]
        new_name = new_name + new_extension if new_name[-1] == "." else new_name + "." + new_extension
    elif filename[-len(new_extension)] == new_extension:
         pass
    else:
        new_name = filename + new_extension if filename[-1] == "." else filename + "." + new_extension
    return new_name

class ManifoldAPDatabase(collections.abc.MutableMapping):
    """
    This class wraps both a human readable JSON file object and a shelve object,
    and keeps them consistent.

    There is some Manifold-specific modifications here as well, but for the most part
    the actual object behaves like a shelve object, which mostly behaves like a
    dictionary.
    """
    def __init__(self, filename, other_filename=None, timing=False):
        """
        The filename should refer to a file containing either a shelve object or to a
        JSON file with an array of encoded ManifoldAP objects. The file can also not yet
        exist, in which case one will be created with the name. The other_filename
        parameter allows one to pass in a name for the JSON or shelve object to be
        created. By default the newly created file will be named by appending .json or
        .shelve. It will rip off the old extension of .json or .shelve for the new one
        though.
        """
        if timing: time0 = time.perf_counter()
        self._original_filename = filename
        if looks_like_a_json_file(filename):
            self._json_filename = filename
            if other_filename:
                self._shelve_filename = other_filename
            else:
                self._shelve_filename = change_file_extension(filename, "json", "shelve")
            with open(filename, 'r') as fp, shelve.open(self._shelve_filename) as shelve_object:
                temp_dict = json_file_to_dict(fp)
                for key in temp_dict: shelve_object[key] = temp_dict[key]
        elif looks_like_a_shelve_file(filename):
            if timing: shelve_time0 = time.perf_counter()
            print("Took", shelve_time0 - time0, "seconds to evaluate the shelve conditional.")
            self._shelve_filename = filename
            if other_filename:
                self._json_filename = other_filename
            else:
                self._json_filename = change_file_extension(filename, "shelve", "json")
        elif not os.path.isfile(filename):
            self._shelve_filename = filename
            self._json_filename = change_file_extension(filename, "shelve", "json")
        else:
            raise RuntimeError("The database couldn't be created")
        if timing:
            time1 = time.perf_counter()
            print("Took", time1 - time0, "seconds to do startup before opening shelve object")
        self._shelve_object = shelve.open(self._shelve_filename, writeback=True)
        if timing:
            time2 = time.perf_counter()
            print("Took", time2 - time1, "seconds to open the shelve object")
    
    def __getitem__(self, key):
            return self._shelve_object[key]
    
    def __setitem__(self, key, value):
        self._shelve_object[key] = value
    
    def __delitem__(self, key):
        del self._shelve_object[key]
    
    def __len__(self):
        return len(self._shelve_object)
    
    def __iter__(self):
        return iter(self._shelve_object)

    def update_shelve(self):
        self._shelve_object.sync()
    
    def update_json(self):
        # It is perhaps a little nonstandard to have a "JSON stream", but I would like
        # to avoid having all the manifolds in memory before we serialize them.
        # This is probably kind of hacky right now, and there should be a proper
        # way done in the module with the custom json classes.
        with open(self._json_filename, 'w') as fp:
            s = str()
            for elt in self._shelve_object.values():
                s += json.dumps(elt, cls=json_encoder.ManifoldAP_Encoder, indent=4) + ",\n"
            t = str()
            for line in s.splitlines():
                t += '    ' + line + '\n'
            t = "[\n" + t[:-2] + "\n]"
            fp.write(t)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.update_json()
        self.update_shelve()
        self._shelve_object.close()
