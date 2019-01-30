import os
from platform import platform

import numpy as np

from ctypes import *

def _find_library(name, version):
    """Find the shared library, checking SEMEL_LIB_PATH first, then system paths."""
    if platform().startswith('macOS'):
        filename = "lib%s.%s.dylib" % (name, version)
    elif platform().startswith('Linux'):
        filename = "lib%s.so.%s" % (name, version)
    else:
        return None

    lib_path = os.environ.get('SEMEL_LIB_PATH')
    if lib_path:
        full = os.path.join(lib_path, filename)
        if os.path.isfile(full):
            return full

    return filename

def lib_load(name, version):
    path = _find_library(name, version)
    if path is None:
        raise RuntimeError("Unsupported platform: %s" % platform())
    lib = cdll.LoadLibrary(path)
    lib.s_init()
    return lib

def lib_unload(lib):
    lib.s_free()
    return None

def dict_keys_str_to_int(x):
    if not isinstance(x, dict):
        return x
    x_ = {}
    for key in x.keys():
        if (len(key) > 0):
            x_[int(key)] = np.array(x[key], dtype=np.float64)
    return x_

def numpy_to_pointer(x):
    [px, nx, mx, tx] = [None, 0, 0, None]
    if isinstance(x, np.ndarray):
        tx = x.dtype

        if (tx == np.float64):
            px = x.flatten().ctypes.data_as(POINTER(c_double))
        elif (tx == np.float32):
            px = x.flatten().ctypes.data_as(POINTER(c_float))
        elif (tx == np.uint64):
            px = x.flatten().ctypes.data_as(POINTER(c_ulong))
        elif (tx == np.uint32):
            px = x.flatten().ctypes.data_as(POINTER(c_uint))
        elif (tx == np.int64):
            px = x.flatten().ctypes.data_as(POINTER(c_long))
        elif (tx == np.int32):
            px = x.flatten().ctypes.data_as(POINTER(c_int))
    return [px, x.shape, tx]

def pointer_to_numpy(px, xshape, tx):
    x = np.ctypeslib.as_array(px, shape=xshape).astype(tx)
    return x

def primitive_value(x):
    x_ = x
    if isinstance(x, (c_double, c_float, c_ulong, c_uint, c_long, c_int)):
        x_ = x.value
    return x_
