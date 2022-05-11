#!/usr/bin/env cctbx.python
from __future__ import print_function
"""
xds_plot_integrate.py

(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
"""
TODO: plot differences in direct beam and rotation axis
"""

import sys
import re
import collections

from cctbx import sgtbx

class IntegrateLp:
    def __init__(self, lpin):
        if lpin is not None:
            self.parse(lpin)
    # __init__()

    def parse(self, int_lp):
        re_im = re.compile("^ (.....)   0 +([0-9\.]+) +([0-9]+) +([0-9]+) +([0-9]+) +([0-9]+) +([0-9]+) +([0-9\.]+) +([0-9\.]+)")
        re_cell = re.compile("^ UNIT CELL PARAMETERS *([0-9\.]+) *([0-9\.]+) *([0-9\.]+) *([0-9\.]+) *([0-9\.]+) *([0-9\.]+)")
        re_rotation = re.compile("^ CRYSTAL ROTATION OFF FROM INITIAL ORIENTATION *([-0-9\.]+) *([-0-9\.]+) *([-0-9\.]+)") #
        re_mosaicity = re.compile("^ CRYSTAL MOSAICITY \(DEGREES\) *([0-9\.]+)") #
        re_axis = re.compile("^ LAB COORDINATES OF ROTATION AXIS *([-0-9\.]+) *([-0-9\.]+) *([-0-9\.]+)") #
        re_beam = re.compile("^ DIRECT BEAM COORDINATES \(REC\. ANGSTROEM\) *([-0-9\.]+) *([-0-9\.]+) *([-0-9\.]+)") #
        re_dist = re.compile("^ CRYSTAL TO DETECTOR DISTANCE \(mm\) *([-0-9\.]+)")
        re_dev_spot = re.compile("^ STANDARD DEVIATION OF SPOT    POSITION \(PIXELS\) *([0-9\.]+)")
        re_dev_spindle = re.compile("^ STANDARD DEVIATION OF SPINDLE POSITION \(DEGREES\) *([0-9\.]+)")
        re_orig = re.compile("^ DETECTOR ORIGIN \(PIXELS\) AT *([0-9\.]+) *([0-9\.]+)")

        images = [] # as key of params
        self.cell_changes = []
        self.blockparams = collections.OrderedDict()
        clear_flag = False

        self.frames = []
        self.scales, self.overloads, self.rejecteds, self.sigmads, self.sigmars = [], [], [], [], []

        self.space_group = None

        # Read INTEGRATE.LP file
        for l in open(int_lp):
            r_im = re_im.search(l)
            r_cell = re_cell.search(l)
            r_rotation = re_rotation.search(l)
            r_dist = re_dist.search(l)
            r_spot = re_dev_spot.search(l)
            r_spindle = re_dev_spindle.search(l)
            r_orig = re_orig.search(l)

            if l.startswith(" SPACE_GROUP_NUMBER="):
                sgnum = int(l.strip().split()[-1])
                if sgnum > 0:
                    self.space_group = sgtbx.space_group_info(sgnum).group()

            if r_im:
                if clear_flag:
                    images = []
                    clear_flag = False
                image, scale, nbkg, novl, newald, nstrong, nrej, sigmad, sigmar = r_im.groups()
                images.append(int(image))

                # for plot
                self.frames.append(int(image))
                self.scales.append(scale)
                self.overloads.append(int(novl))
                self.rejecteds.append(int(nrej))
                self.sigmads.append(sigmad)
                self.sigmars.append(sigmar)

            if r_cell:
                #a, b, c, alpha, beta, gamma = r_cell.groups()
                self.blockparams.setdefault(tuple(images), {})["cell"] = r_cell.groups()
                self.cell_changes.append((images, r_cell.groups()))
                clear_flag = True

            if r_rotation:
                self.blockparams.setdefault(tuple(images), {})["rotation"] = r_rotation.groups()
                clear_flag = True

            if r_dist:
                self.blockparams.setdefault(tuple(images), {})["dist"] = r_dist.group(1)
                clear_flag = True
            if r_spot:
                self.blockparams.setdefault(tuple(images), {})["spot"] = r_spot.group(1)
                clear_flag = True
            if r_spindle:
                self.blockparams.setdefault(tuple(images), {})["spindle"] = r_spindle.group(1)
                clear_flag = True
            if r_orig:
                self.blockparams.setdefault(tuple(images), {})["orig"] = r_orig.groups()
                clear_flag = True

            if l.startswith(" SIGMAB (degree)"):
                self.blockparams.setdefault(tuple(images), {})["sigmab9"] = l.strip().split()[-9:]
                clear_flag = True

            if l.startswith(" SIGMAR (degree)"):
                self.blockparams.setdefault(tuple(images), {})["sigmar9"] = l.strip().split()[-9:]
                clear_flag = True

    # parse_integrate_lp()
# class IntegrateLp


class CellConstraints:
    def __init__(self, space_group):
        self.cs = space_group.crystal_system()
    # __init__()

    def is_b_equal_a(self): return self.cs in ("Tetragonal", "Hexagonal", "Trigonal", "Cubic")
    def is_c_equal_a_b(self): return self.cs == "Cubic"

    def is_angle_constrained(self, angle):
        assert angle in ("alpha", "beta", "gamma")
        if self.cs == "Triclinic": return False
        if self.cs == "Monoclinic": return angle != "beta"

        return True
    # is_angle_constrained()
# class CellConstraints


def make_plot(lp, log_out):
    ofs = open(log_out, "w")

    ofs.write("$TABLE: Parameters estimated for each frame:\n")
    ofs.write("$GRAPHS\n")
    ofs.write(":scales")
    ofs.write(":A:1,2:\n")
    ofs.write(":number of overloaded reflections")
    ofs.write(":A:1,3:\n")
    ofs.write(":number of unexpected reflections")
    ofs.write(":A:1,4:\n")
    ofs.write(":SIGMAB (beam divergence e.s.d.)")
    ofs.write(":A:1,5:\n")
    ofs.write(":SIGMAR (reflecting range e.s.d.)")
    ofs.write(":A:1,6:\n")
    ofs.write("$$\n")
    ofs.write("Frame scale overlods nrej sigmaD sigmaM $$\n$$\n")
    for f, scale, novl, nrej, sd, sm in zip(lp.frames, lp.scales, lp.overloads, lp.rejecteds, lp.sigmads, lp.sigmars):
        ofs.write("%5d %s %d %d %s %s\n" % (f, scale, novl, nrej, sd, sm))

    ofs.write("$$\n")
    ofs.write("\n\n\n")

    ofs.write("$TABLE: Parameters estimated for each block:\n")
    ofs.write("$GRAPHS\n")
    ofs.write(":unit cell length a")
    ofs.write(":A:1,2:\n")

    cellconstr = CellConstraints(lp.space_group)
    if not cellconstr.is_b_equal_a():
        ofs.write(":unit cell length b")
        ofs.write(":A:1,3:\n")
    if not cellconstr.is_c_equal_a_b():
        ofs.write(":unit cell length c")
        ofs.write(":A:1,4:\n")
    if not cellconstr.is_angle_constrained("alpha"):
        ofs.write(":unit cell angle alpha")
        ofs.write(":A:1,5:\n")
    if not cellconstr.is_angle_constrained("beta"):
        ofs.write(":unit cell angle beta")
        ofs.write(":A:1,6:\n")
    if not cellconstr.is_angle_constrained("gamma"):
        ofs.write(":unit cell angle gamma")
        ofs.write(":A:1,7:\n")
    ofs.write(":rotations off from initial orientation")
    ofs.write(":A:1,8,9,10:\n")
    ofs.write(":distance")
    ofs.write(":A:1,11:\n")
    ofs.write(":deviations from predicted positions")
    ofs.write(":A:1,12,13:\n")
    ofs.write(":beam center")
    ofs.write(":A:1,14,15:\n")
    ofs.write("$$\n")
    ofs.write("#image a b c alpha beta gamma rotx roty rotz dist spot spindle orgx orgy$$\n$$\n")
    for images, param in sorted(lp.blockparams.items()):
        for i in images:
            print("%4d " % i, "  ".join(param.get("cell", ["D"]*6)), " ".join(param.get("rotation", ["D"]*3)), param.get("dist","D"), param.get("spot","D"), param.get("spindle","D"), " ".join(param.get("orig",["D"]*2)), file=ofs)

    ofs.write("$$\n")
    ofs.write("\n\n\n")


    ofs.write("$TABLE: sigmaB and sigmaR on 9 areas for each block:\n")
    ofs.write("$GRAPHS\n")
    ofs.write(":SIGMAB")
    ofs.write(":A:1,2,3,4,5,6,7,8,9,10:\n")
    ofs.write(":SIGMAR")
    ofs.write(":A:1,11,12,13,14,15,16,17,18,19:\n")
    ofs.write("$$\n")
    ofs.write("#image %s %s$$\n$$\n" % (" ".join(["sigmab%d"%x for x in range(1,10)]), " ".join(["sigmar%d"%x for x in range(1,10)])))
    for images, param in sorted(lp.blockparams.items()):
        for i in images:
            print("%4d " % i, " ".join(param["sigmab9"]), " ".join(param["sigmar9"]), file=ofs)

    ofs.write("$$\n")
    ofs.write("\n\n\n")
# make_plot()

def run(int_lp, log_out="plot_integrate.log"):
    lpobj = IntegrateLp(int_lp)
    make_plot(lpobj, log_out)
# run()

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        int_lp = sys.argv[1]
    else:
        int_lp = "INTEGRATE.LP"

    log_out = "plot_integrate.log"

    run(int_lp, log_out)

    print()
    print("Run:")
    print("loggraph", log_out)
