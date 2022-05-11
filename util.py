"""
util.py

(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function

import os
import sys
import shutil
import subprocess
import glob

def call(cmd, arg="",
         stdin=None, stdout=subprocess.PIPE,
         wdir=None,
         expects_in=[], expects_out=[]):

    ##
    # call the external program using subprocess.
    #
    # @param expects_in expected files before running
    # @param expects_out expected files after running
    #
    # expected_in/out must be written as relative path from wdir or absolute path.
    #

    def check_exist(files):
        is_exist = [os.path.isfile(f) for f in files]

        if sum(is_exist) != len(is_exist):
            not_founds = [ f for f, e in zip(files, is_exist) if not e ]
            raise Exception("Expected file(s) not found: " + " ".join(not_founds))

    # check_exist()

    if wdir is None:
        wdir = os.getcwd()

    # Go to working directory
    cwd = os.getcwd()
    os.chdir(wdir)

    # check before run
    check_exist(expects_in)

    # call the program
    p = subprocess.Popen("%s %s" % (cmd, arg),
                         shell=True,
                         stdin=subprocess.PIPE,
                         stdout=stdout,
                         stderr=stdout,
                         universal_newlines=True)

    if stdin is not None:
        p.stdin.write(stdin)    

    if stdout == subprocess.PIPE:
        out, err = p.communicate()
    else:
        out, err = None, None
        p.stdin.close()
        p.wait()

    if p.returncode < 0:
        print(cmd, ": returncode is", p.returncode, file=sys.stderr)

    # check after run
    check_exist(expects_out)

    # go back to the previous working directory
    os.chdir(cwd)

    return p.returncode, out, err
# call()

def rotate_file(filename, copy=False):
    """
    Rotate file like logrotate.
    If given filename already exists, rename it to "filename".n, n=1...
    Filename with larger n is older one.
    """

    # If not exist,
    if not os.path.isfile(filename):
        return

    # make list [ [filename, number], ... ]
    old_list = []
    dot_files = glob.glob(filename + ".*")
    for f in dot_files:
        suffix = f.replace(filename+".", "")
        try:
            i = int(suffix)
            if str(i) == suffix: # ignore if suffix was such as 003...
                old_list.append([f, i])
        except ValueError as e:
            continue

    old_list.sort(lambda x,y: x[1]-y[1])

    # rotate files
    for f, i in reversed(old_list):
        os.rename(f, "%s.%d" % (f[:f.rfind(".")], i+1))

    if copy:
        shutil.copyfile(filename, filename + ".1")
    else:
        os.rename(filename, filename + ".1")

    return filename + ".1"
# rotate_file()

def get_number_of_processors(default=4):
    nproc = default

    if os.path.isfile("/proc/cpuinfo"):
        nproc = len([x for x in open("/proc/cpuinfo") if x.startswith("processor")])
    else:
        try:
            nproc = int(subprocess.getoutput("sysctl -n hw.ncpu"))
        except:
            pass

    return nproc
# get_number_of_processors()

def safe_float(v):
    try:
        return float(v)
    except ValueError:
        return float("nan")
# safe_float()

def num_th_str(v):
    s = str(v)
    if s[-1] == "1": return s+"st"
    if s[-1] == "2": return s+"nd"
    if s[-1] == "3": return s+"rd"
    return s+"th"
# num_th_str()
