"""
load_xparm.py

(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import os, shutil, tempfile, subprocess
import numpy

from pymol import cmd, stored
from pymol.cgo import *    # get constants


"""
XPARAM.XDS

     1       50.0000    1.0000 -0.999997 -0.001378 -0.002095
       1.071460      -0.001051      -0.000103       0.933305
      3072      3072    0.073242    0.073242
     159.917953    1541.014771    1538.410645
       1.000000       0.000000       0.000000
       0.000000       1.000000       0.000000
       0.000000       0.000000       1.000000
    18     75.7997     97.0630    202.3445  90.000  90.000  90.000
      74.986534      -7.994404       7.661397
      14.166846      72.044411     -63.483360
      -1.222476     133.907318     151.692612

1. Starting image number (STARTING_FRAME=), spindle angle at start (STARTING_ANGLE=), oscillation range, and laboratory coordinates of the rotation axis.
2.  Wavelength (A) and laboratory coordinates of the incident beam wavevector.
3.Number of pixels along the detector X-axis (NX=) and Y-axis (NY=) in a data image and pixel sizes (mm) (QX=, QY=) along X and Y.
4. Signed distance between crystal and detector (mm), detector X-coordinate (pixels) of origin, detector Y-coordinate (pixels) of origin.
5. Laboratory coordinates of the unit vector along the detector X-axis.
6. Laboratory coordinates of the unit vector along the detector Y-axis.
7. Laboratory coordinates of the unit vector along the detector normal.
8. Space group number and unit cell parameters (A and degrees).
9. Laboratory coordinates of the unit cell a-axis of the unrotated crystal.
10. Laboratory coordinates of the unit cell b-axis of the unrotated crystal.
11. Laboratory coordinates of the unit cell c-axis of the unrotated crystal.
"""

class XPARM:
    def __init__(self, xparm_file):
        lines = open(xparm_file).readlines()

        is_new_format = "XPARM.XDS" in lines[0]

        if not is_new_format:
            starting_frame, starting_angle, osc_range, rotx, roty, rotz = lines[0].split()
            wavelength, ibeamx, ibeamy, ibeamz = lines[1].split()
            nx, ny, qx, qy = lines[2].split()
            distance, orgx, orgy = lines[3].split()
            Xx, Xy, Xz = lines[4].split()
            Yx, Yy, Yz = lines[5].split()
            Zx, Zy, Zz = lines[6].split()
            spacegroup, a, b, c, alpha, beta, gamma = lines[7].split()
            ax, ay, az = lines[8].split()
            bx, by, bz = lines[9].split()
            cx, cy, cz = lines[10].split()
        else:
            starting_frame, starting_angle, osc_range, rotx, roty, rotz = lines[1].split()
            wavelength, ibeamx, ibeamy, ibeamz = lines[2].split()
            spacegroup, a, b, c, alpha, beta, gamma = lines[3].split()
            ax, ay, az = lines[4].split()
            bx, by, bz = lines[5].split()
            cx, cy, cz = lines[6].split()
            nseg, nx, ny, qx, qy = lines[7].split()
            orgx, orgy, distance = lines[8].split()
            Xx, Xy, Xz = lines[9].split()
            Yx, Yy, Yz = lines[10].split()
            Zx, Zy, Zz = lines[11].split()

        self.starting_frame = int(starting_frame)
        self.starting_angle = float(starting_angle)
        self.osc_range = float(osc_range)
        self.rotation_axis = numpy.array((float(rotx), float(roty), float(rotz)))
        self.wavelength = float(wavelength)
        self.incident_beam = numpy.array((float(ibeamx), float(ibeamy), float(ibeamz)))
        self.nx = float(nx)
        self.ny = float(ny)
        self.qx = float(qx)
        self.qy = float(qy)
        self.distance = float(distance)
        self.origin = numpy.array((float(orgx), float(orgy)))
        self.X_axis = numpy.array((float(Xx), float(Xy), float(Xz)))
        self.Y_axis = numpy.array((float(Yx), float(Yy), float(Yz)))
        self.Z_axis = numpy.array((float(Zx), float(Zy), float(Zz)))
        self.spacegroup = int(spacegroup)
        self.unit_cell = numpy.array((float(a), float(b), float(c), float(alpha), float(beta), float(gamma)))
        self.a_axis = numpy.array((float(ax), float(ay), float(az)))
        self.b_axis = numpy.array((float(bx), float(by), float(bz)))
        self.c_axis = numpy.array((float(cx), float(cy), float(cz)))
    # __init__()
# class XPARM

def reciprocal(a, b, c):
    ##
    # @param a,b,c numpy.array
    ##

    V = numpy.dot(numpy.cross(a,b),c)

    astar = numpy.cross(b,c) / V
    bstar = numpy.cross(c,a) / V
    cstar = numpy.cross(a,b) / V

    return astar, bstar, cstar
# reciprocal()




def rot_matrix(axis, angle):
    cos = numpy.cos(angle)
    sin = numpy.sin(angle)
    Vx, Vy, Vz = axis

    M = numpy.matrix([[Vx*Vx*(1-cos) + cos,	Vx*Vy*(1-cos) - Vz*sin,	Vz*Vx*(1-cos) + Vy*sin],
                      [Vx*Vy*(1-cos) + Vz*sin,	Vy*Vy*(1-cos) + cos,	Vy*Vz*(1-cos) - Vx*sin],
                      [Vz*Vx*(1-cos) - Vy*sin,	Vy*Vz*(1-cos) + Vx*sin,	Vz*Vz*(1-cos) + cos]]
                     )

    return M

# rot_matrix


def obj_axis(o, v, abc, radius=0.3):
    # arrow offset
    length = numpy.linalg.norm(v-o)
    if length >  1:
        ao = o + (v-o)/length
    else:
        ao = v
        radius = length

    if abc == "a":
        col = [1., .3, .3] # red
    elif abc == "b":
        col = [.3, 1., .3] # green
    else:
        col = [1., 1., .3] # yellow

    obj = [ CYLINDER,
            o[0], o[1], o[2],      # XYZ 1
            v[0]-ao[0], v[1]-ao[1], v[2]-ao[2],      # XYZ 2
            radius,                # Radius
            col[0],col[1],col[2],         # RGB Color 1
            col[0],col[1],col[2],         # RGB Color 2

            CONE,
            v[0]-ao[0], v[1]-ao[1], v[2]-ao[2], # XYZ 1
            v[0], v[1], v[2],                   # XYZ 2
            radius*1.5,                         # Radius 1
            0.0,                                # Radius 2
            col[0],col[1],col[2],               # RGB Color 1
            col[0],col[1],col[2],               # RGB Color 2
            1.0, 1.0,                           # Caps 1 & 2
            ]
    return obj
# obj_axis()

def obj_phigonio(o, rot_ax, scale):
    col = [1., 1., 1.]

    ra = -rot_ax * scale

    obj = [CONE,
           ra[0], ra[1], ra[2], # XYZ 1
           o[0], o[1], o[2],    # XYZ 2
           2,                   # Radius 1
           0,                   # Radius 2
           col[0],col[1],col[2],               # RGB Color 1
           col[0],col[1],col[2],               # RGB Color 2
           1.0, 1.0
           ]

    return obj
# obj_phigonio()

def obj_beam(o, beam, l, scale):
    s = o - beam/numpy.linalg.norm(beam) * 4. / l * scale
    e = o - beam/numpy.linalg.norm(beam) * 1. / l * scale

    obj = [BEGIN, LINES,
           COLOR, 1.0, 1.0, 1.0,

           VERTEX, s[0],s[1],s[2],
           VERTEX, e[0],e[1],e[2],
           END
           ]

    return obj
# obj_beam()

def obj_ewald(o, beam, l, scale):
    ##
    # @param o origin
    # @param beam indicent_beam
    # @param l lambda
    ##

    #l = 1.0

    pos = o - beam/numpy.linalg.norm(beam)/l * scale

    obj = [ALPHA, 0.3,
           SPHERE,
           pos[0], pos[1], pos[2],  # XYZ
           1/l*scale,  # Radius
           ]

    return obj
# obj_ewald()

def load_ewald(o, beam, l, scale, name="Ewald"):
    ##
    # As pseudoatom
    # @param o origin
    # @param beam indicent_beam
    # @param l lambda
    ##

    #l = 1.0

    pos = o - beam/numpy.linalg.norm(beam)/l * scale

    print "lambda=", l
    print "pos=", pos

    cmd.pseudoatom(object=name,
                   pos=list(pos),
                   vdw=1./l*scale,
                   )
    cmd.hide("everything", name)
    cmd.show("spheres", name)
    cmd.set("sphere_transparency", 0.5)


# load_ewald()


def draw(xparmfile, end, offset="0"):
    """
DESCRIPTION
    Visualize XPARM.XDS

USAGE
    load_xparm XPARM.XDS, framenumber [,frame offset]
    end is the end number of frames

AUTHORS
    Keitaro Yamashita, 2011
    """

    end = int(end)
    offset = int(offset)

    if not hasattr(stored, "__xparm_internal__"):
        stored.__xparm_internal__ = 1
    else:
        stored.__xparm_internal__ += 1

    suffix = str(stored.__xparm_internal__)

    param = XPARM(xparmfile)
    origin = numpy.zeros(3)

    gonio  = obj_phigonio(origin, param.rotation_axis, scale=6)
    beam   = obj_beam(origin, param.incident_beam, param.wavelength, scale=10)

    cmd.load_cgo(gonio, "gonio_"+suffix)
    cmd.load_cgo(beam,  "beam_"+suffix)

    load_ewald(origin, param.incident_beam, param.wavelength, scale=10, name="Ewald_"+suffix)

    # Rotate space by phi axis
    for i, d in enumerate([param.starting_angle + (i+0.5+offset) * param.osc_range
                           for i in xrange(end)]
                          ):
        print i, d
        R = rot_matrix(param.rotation_axis, numpy.radians(d))

        A = numpy.array(numpy.dot(R, param.a_axis))[0]
        B = numpy.array(numpy.dot(R, param.b_axis))[0]
        C = numpy.array(numpy.dot(R, param.c_axis))[0]

        rA, rB, rC = reciprocal(A, B, C)

        #print "real vectors:", A, B, C
        #print "reciprocal vectors:", rA, rB, rC

        scale = .5 * .2 * min([numpy.linalg.norm(A), numpy.linalg.norm(B), numpy.linalg.norm(C)]) / max([numpy.linalg.norm(rA), numpy.linalg.norm(rB), numpy.linalg.norm(rC)])

        # Unit vectors
        a_axis = obj_axis(origin, A*.2, "a")
        b_axis = obj_axis(origin, B*.2, "b")
        c_axis = obj_axis(origin, C*.2, "c")

        # Reciprocal unit vectors
        as_axis = obj_axis(origin, rA*scale, "a")
        bs_axis = obj_axis(origin, rB*scale, "b")
        cs_axis = obj_axis(origin, rC*scale, "c")

        cmd.load_cgo(a_axis+b_axis+c_axis, "axes_"+suffix, i+1)
        cmd.load_cgo(as_axis+bs_axis+cs_axis, "reciprocal_axes_"+suffix, i+1)

# draw()

cmd.extend("load_xparm", draw)
