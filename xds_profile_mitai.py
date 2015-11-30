#!/usr/bin/env python
"""
xds_profile_mitai.py

(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import os, sys, math, re

def write_xplor_map(res, NX, NZ, filename):
    ofs = open(filename, "w")

    ofs.write("\n")
    ofs.write("      2 !NTITLE\n")
    ofs.write("REMARKS COMMENT1\n")
    ofs.write("REMARKS COMMENT2\n")
    ofs.write(("%8d%8d%8d" %(NX-1,0,NX-1))*2 + "%8d%8d%8d" %(NZ-1,0,NZ-1) + "\n")
    ofs.write("%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e\n" % (NX,NX,NZ,90,90,90))
    ofs.write("ZYX\n")

    acc = []
    for i in xrange(NZ):
        if i > 0:
                ofs.write("\n")
        ofs.write("%8d\n" % i)
        count = 0
        for j in xrange(NX):
            for k in xrange(NX):
                try:
                    ofs.write("%12.5e" % res[i][j][k])
                    acc.append(res[i][j][k])
                except:
                    print "EXCEPTION!=", i, j, k
                count += 1
                if count > 5:
                    ofs.write("\n")
                    count = 0

    ofs.write("\n")
    ofs.write("%8d\n" % -9999)
    mean = float(sum(acc))/len(acc)
    stddev = math.sqrt(sum(map(lambda x:(x-mean)**2, acc))/len(acc))
    ofs.write("%12.4f %12.4f \n" % (mean, stddev))

# write_xplor_map()

def run(intlp):
    lines = []
    profile_ids = []

    read_flag = False
    prefix = "profile"

    ofs_pml = open("load_%s.pml" % prefix, "w")

    for l in open(intlp):
        if "PROCESSING OF IMAGES " in l:
            from_to = map(lambda x:int(x.strip()), l.replace("PROCESSING OF IMAGES","").split("..."))
            print "images", from_to
            profile_ids.append("image%.3d-%.3d" % tuple(from_to))

        if "***** AVERAGE THREE-DIMENSIONAL PROFILE" in l or "RUN-AVERAGE OF PROFILE #" in l:
            r = re.search("\*\*\*\*\* RUN-AVERAGE OF PROFILE #  ([1-9]) \*\*\*\*\*", l)
            if r:
                print "run-average", r.group(1)
                profile_ids.append("run_average_%s" % r.group(1))
            read_flag = True
            lines.append([])
            continue

        if "REFLECTION INTENSITIES INTEGRATED BY PROFILE FITTING" in l:
            read_flag = False
            continue

        if read_flag and re.search("[^-0-9 ]", l.rstrip()):
            read_flag = False
            continue

        if read_flag:
            lines[-1].append(l)

    print
    for i, (ll,pid) in enumerate(zip(lines, profile_ids)):
        vals = []
        for l in ll:
            if l.strip() == "" and (len(vals) == 0 or len(vals[-1]) != 0):
                vals.append([])
                continue

            sp = re.findall("...", l.rstrip())
            tmp = []
            for s in sp:
                if s.strip() == "":
                    tmp.append([])
                else:
                    tmp[-1].append(int(s))

            if len(tmp) > 0:
                vals[-1].append(tmp)

        if len(vals) > 0 and len(vals[-1]) == 0:
            del vals[-1]

        nx = len(vals[0][0][0])
        nz = sum(map(lambda x:len(x[0]), vals))
        print pid, "nx, nz=", nx, nz
        data = [[[None for x in xrange(nx)] for y in xrange(nx)] for z in xrange(nz)] # [z][y][x]

        z,y,x = 0,0,0
        for iv, v in enumerate(vals):
            y = 0
            for ivv, vv in enumerate(v):
                z = len(vals[0][0]) * iv
                for vvv in vv:
                    assert len(data[z][y]) == len(vvv)
                    data[z][y] = vvv
                    #print z,y,data[z][y]
                    z += 1
                y += 1

        filename = "%s_%s_%.4d.xplor" % (prefix, pid, i+1)
        write_xplor_map(data, nx, nz, filename)
        ofs_pml.write("load %s, %s_%s_%d\n"%(filename, prefix, pid, i+1))
        ofs_pml.write("isomesh msh_%d, %s_%s_%d\n"%(i+1, prefix, pid, i+1))

    print
    print "Start:"
    print "pymol load_profile.pml"
# run()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: %s INTEGRATE.LP" % sys.argv[0]
        sys.exit()

    intlp_in = sys.argv[1]
    run(intlp_in)
