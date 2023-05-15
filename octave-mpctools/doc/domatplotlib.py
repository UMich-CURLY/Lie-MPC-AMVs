import numpy as np
import matplotlib.pyplot as plt # analysis:ignore
import scipy.io

def loadmat(filename):
    """
    Loads the given mat file, transforming structs to dicts.
    
    Note that if any structs are nonscalar, you may get unexpected results.
    """
    rawdata = scipy.io.loadmat(filename, squeeze_me=True)
    data = {}
    for (k, v) in rawdata.items():
        if not k.startswith("__"):
            data[k] = _struct_to_dict(v)
    return data
 
def _struct_to_dict(obj):
    """Converts the object to a Python dict if it is a struct."""
    if isinstance(obj, np.ndarray):
        if obj.dtype.fields is not None:
            obj = {k : obj[k][()] for k in obj.dtype.fields}
        elif obj.dtype == np.dtype("O"):
            for x in np.nditer(obj, op_flags=["readwrite"], flags=["refs_ok"]):
                x[...] = _struct_to_dict(x[()])
    return obj
