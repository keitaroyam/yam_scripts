#!/usr/bin/env cctbx.python
"""
xds_merge_framecbf.py

(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function

import pycbf
import os
import numpy

def load_minicbf_as_numpy(filein, quiet=True): # This can also read XDS special cbf
    assert os.path.isfile(filein)
    if not quiet:
        print("reading", filein, "as minicbf")
    h = pycbf.cbf_handle_struct()
    h.read_file(filein, pycbf.MSG_DIGEST)
    h.require_category("array_data")
    h.find_column("data")
    compression, binary_id, elsize, elsigned, elunsigned, elements, minelement, maxelement, bo, ndimfast, ndimmid, ndimslow, padding = h.get_integerarrayparameters_wdims()
    assert elsize == 4 or elsize == 8
    assert elsigned == 1
    assert ndimslow <= 1
    arr = numpy.fromstring(h.get_integerarray_as_string(), dtype=numpy.int32 if elsize==4 else numpy.int64)
    return arr, ndimfast, ndimmid

# load_minicbf_as_numpy()


def save_numpy_data_as_cbf(data, size1, size2, title, cbfout, pilatus_header=None):
    h = pycbf.cbf_handle_struct()
    h.new_datablock(title)

    h.require_category('array_data')

    if pilatus_header is not None:
        h.require_column('header_convention')
        h.set_value('"PILATUS_1.2"')
        h.require_column('header_contents')
        h.set_value(pilatus_header)


    h.require_category('array_data')
    h.require_column('data')

    elsigned = 1
    if data.dtype in (numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64):
        elsigned = 0

    h.set_integerarray_wdims_fs(pycbf.CBF_BYTE_OFFSET, 1, data.tostring(), data.dtype.itemsize,
                                elsigned, len(data), "little_endian",
                                size1, size2, 1, 0)
    h.write_file(cbfout, pycbf.CBF,
                 pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, pycbf.ENC_NONE)
# save_numpy_data_as_cbf()


def run(files, cbfout):
    merged = None
    for i, f in enumerate(files):
        repl = -10 * (i+1)
        print("%s %d" % (f, repl))

        data, ndimfast, ndimmid = load_minicbf_as_numpy(f)
        if i == 0:
            merged = data.copy()
            continue

        merged[data==-10] = 65540 # repl # for adxv visualization. only assuming two files.

    save_numpy_data_as_cbf(merged, ndimfast, ndimmid, "merged_predictions", cbfout)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 3:
        print("Usage: %s FRAME.1.cbf FRAME.2.cbf .." % sys.argv[0])
        quit()

    run(sys.argv[1:], "FRAME_merged.cbf")
