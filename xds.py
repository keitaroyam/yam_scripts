"""
xds.py

(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import re, os, glob, shutil, random, string

re_xds_kwd = re.compile("([^ =]+)= *((?:(?! [^ =]+=).)*)")

def get_xdsinp_keyword(xdsinp):
    ##
    # Return the list of tuple (keyword, value) in XDS.INP
    #

    #    re_kwd = re.compile("([^ =]+)=(.*)(?: [^ =]+=)")
    #    re_kwd = re.compile("([^ =]+)=(.*)[^ =]+=")
    res = []
    for l in open(xdsinp):
        l = l[:l.find("!")] # Remove comment
        r = re_xds_kwd.findall(l)
        res.extend(r)

    return res
# get_xdsinp_keyword()


def modify_xdsinp(xdsinp, inp_params):
    """
    inp_params = [(key, val), (key, val), ...]
    If val is None, does not insert key=val line, but comment out original speicifications (thus use default).
    """

    kwds = get_xdsinp_keyword(xdsinp)

    ofs = open(xdsinp, "w")
    for k, v in inp_params:
        if v is None:
            continue
        ofs.write(" %s= %s\n" % (k, v))
    ofs.write("\n")

    for kwd, val in kwds:
        if kwd in [k for k,v in inp_params]:
            ofs.write("!")
        ofs.write(" %s= %s\n" % (kwd,val))

    ofs.close()
    #print open(xdsinp).read()
# modify_xdsinp()

def make_backup(backup_needed, bk_prefix=None, wdir=None, quiet=False):
    if wdir is None:
        wdir = os.getcwd()

    if bk_prefix is None:
        while True:
            bk_prefix = "".join(random.choice(string.digits + string.lowercase) for i in xrange(10)) + "_"
            if len(glob.glob(os.path.join(wdir, bk_prefix + "*"))) == 0:
                break

    for f in backup_needed:
        if os.path.isfile(os.path.join(wdir, f)):
            shutil.copyfile(os.path.join(wdir, f), 
                            os.path.join(wdir, bk_prefix+f))
            if not quiet: print "Backup: %s => %s" % (f, bk_prefix+f)

    return bk_prefix
# make_backup()

def remove_backups(backup_needed, bk_prefix, wdir=None):
    if wdir is None:
        wdir = os.getcwd()

    map(lambda x: os.remove(x),
        glob.glob(os.path.join(wdir, bk_prefix+"*")))
# remove_backups()

def revert_files(backup_needed, bk_prefix, wdir=None, quiet=False):
    if wdir is None:
        wdir = os.getcwd()

    for f in backup_needed:
        if os.path.isfile(os.path.join(wdir, bk_prefix+f)):
            os.rename(os.path.join(wdir, bk_prefix+f),
                      os.path.join(wdir, f))
            if not quiet: print "Reverted: %s => %s" % (bk_prefix+f, f)
# revert_files()

def optimal_delphi_by_nproc(xdsinp=None, osc_width=None, nframes=None, nproc=None, min_delphi=5):
    if xdsinp is not None:
        kwds = dict(get_xdsinp_keyword(xdsinp))
        osc_width = float(kwds["OSCILLATION_RANGE"])
        data_range = map(int, kwds["DATA_RANGE"].split())
        nframes = data_range[1] - data_range[0] + 1
    if None in (osc_width, nframes):
        print "osc_width and nframes not specified."
        return None

    if nproc is None:
        from yamtbx import util
        nproc = util.get_number_of_processors()

    # Need to check nframes%nproc

    delphi = nproc * osc_width

    fac = 2
    while delphi < min_delphi:
        delphi = delphi * fac
        fac += 1

    return delphi
# optimal_delphi_by_nproc()
