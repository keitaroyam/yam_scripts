"""
xparm.py

(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import numpy
from xds import get_xdsinp_keyword
from util import safe_float

class XPARM:
    def __init__(self, xparm_file=None):
        if xparm_file is not None:
            self.parse_xparm_file(xparm_file)
        else:
            self.starting_frame = 1
            self.starting_angle = 0.
            self.osc_range = -1.
            self.rotation_axis = numpy.array((1.,0.,0.))
            self.wavelength = -1.
            self.incident_beam = numpy.array((0.,0.,1.))
            self.nx, self.ny = 0, 0
            self.qx, self.qy = 0., 0.
            self.distance = 0.
            self.origin = numpy.array((0., 0.))
            self.X_axis = numpy.array((1., 0., 0.))
            self.Y_axis = numpy.array((0., 1., 0.))
            self.Z_axis = numpy.array((0., 0., 1.))
            self.spacegroup = 1
            self.unit_cell = numpy.array((100., 100., 100., 90., 90., 90.))
            self.a_axis = numpy.array((100., 0., 0.))
            self.b_axis = numpy.array((0., 100., 0.))
            self.c_axis = numpy.array((0., 0., 100.))
    # __init__()

    def parse_xparm_file(self, xparm_file):
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
        self.a_axis = numpy.array(map(safe_float, (ax, ay, az)))
        self.b_axis = numpy.array(map(safe_float, (bx, by, bz)))
        self.c_axis = numpy.array(map(safe_float, (cx, cy, cz)))
    # parse_xparm_file()

    def set_info_from_xdsinp(self, xdsinp):
        # XXX x, y, z axes

        table = [("STARTING_FRAME", "starting_frame", lambda x: int(x)),
                 ("STARTING_ANGLE", "starting_angle", lambda x: float(x)),
                 ("OSCILLATION_RANGE", "osc_range", lambda x: float(x)),
                 ("ROTATION_AXIS", "rotation_axis", lambda x: numpy.array(map(lambda y:float(y), x.split()))),
                 ("X-RAY_WAVELENGTH", "wavelength", lambda x: float(x)),
                 ("INCIDENT_BEAM_DIRECTION", "incident_beam", lambda x: numpy.array(map(lambda y:float(y), x.split()))),
                 ("NX", "nx", lambda x: int(x)),
                 ("NY", "ny", lambda x: int(x)),
                 ("QX", "qx", lambda x: float(x)),
                 ("QY", "qy", lambda x: float(x)),
                 ("DETECTOR_DISTANCE", "distance", lambda x: float(x)),
                 ("SPACE_GROUP_NUMBER", "spacegroup", lambda x: int(x)),
                 ("UNIT_CELL_CONSTANTS", "unit_cell", lambda x: numpy.array(map(lambda y:float(y), x.split()))),
                 ("UNIT_CELL_A-AXIS", "a_axis", lambda x: numpy.array(map(lambda y:float(y), x.split()))),
                 ("UNIT_CELL_B-AXIS", "b_axis", lambda x: numpy.array(map(lambda y:float(y), x.split()))),
                 ("UNIT_CELL_C-AXIS", "c_axis", lambda x: numpy.array(map(lambda y:float(y), x.split())))
                 ]
        inp = dict(get_xdsinp_keyword(xdsinp)) # I believe dict() removes duplicated parameters and keeps last.

        for k, at, f in table:
            if k in inp and inp[k].strip() != "":
                setattr(self, at, f(inp[k]))
        if "ORGX" in inp:
            self.origin[0] = float(inp["ORGX"])
        if "ORGY" in inp:
            self.origin[1] = float(inp["ORGY"])
    # set_info_from_xdsinp()

    def xparm_str(self, old_format=False):
        assert not old_format # Currently, only new format is supported!
        xparm_str = """ XPARM.XDS    VERSION March 30, 2013
%6d%14.4f%10.4f%10.6f%10.6f%10.6f
%15.6f%15.6f%15.6f%15.6f
%6d%12.6f%12.6f%12.6f%8.3f%8.3f%8.3f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
%10d%10d%10d%12.6f%12.6f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
%15.6f%15.6f%15.6f
%10d%10d%10d%10d%10d
%8.2f%8.2f%8.2f%9.5f%9.5f%9.5f%9.5f%9.5f%9.5f
""" % (self.starting_frame, self.starting_angle, self.osc_range, self.rotation_axis[0], self.rotation_axis[1], self.rotation_axis[2], 
       self.wavelength, self.incident_beam[0], self.incident_beam[1], self.incident_beam[2], 
       self.spacegroup, self.unit_cell[0], self.unit_cell[1], self.unit_cell[2], self.unit_cell[3], self.unit_cell[4], self.unit_cell[5],
       self.a_axis[0], self.a_axis[1], self.a_axis[2], 
       self.b_axis[0], self.b_axis[1], self.b_axis[2], 
       self.c_axis[0], self.c_axis[1], self.c_axis[2],
       1, self.nx, self.ny, self.qx, self.qy,
       self.origin[0], self.origin[1], self.distance, 
       self.X_axis[0], self.X_axis[1], self.X_axis[2], 
       self.Y_axis[0], self.Y_axis[1], self.Y_axis[2], 
       self.Z_axis[0], self.Z_axis[1], self.Z_axis[2], 
       1, 1, self.nx, 1, self.ny,
       0., 0., 0., 1., 0., 0., 0., 1., 0.,
       )
        return xparm_str
    # xparm_str()

    def crystal_symmetry(self):
        from cctbx import crystal
        return crystal.symmetry(tuple(self.unit_cell),
                                self.spacegroup)
    # crystal_symmetry()
# class XPARM


def get_xparm_from_integrate_lp(lpfile, frame):
    assert 0 < frame

    keys = {"beam direction": "DIRECT BEAM COORDINATES (REC. ANGSTROEM)",
            "beam center": "DETECTOR ORIGIN (PIXELS) AT",
            "distance": "CRYSTAL TO DETECTOR DISTANCE (mm)",
            "rotation axis": "LAB COORDINATES OF ROTATION AXIS",
            "a axis": "COORDINATES OF UNIT CELL A-AXIS",
            "b axis": "COORDINATES OF UNIT CELL B-AXIS",
            "c axis": "COORDINATES OF UNIT CELL C-AXIS",
            "cell": "UNIT CELL PARAMETERS",
            "spacegroup": "SPACE GROUP NUMBER"
            }
    data = {}
    
    flag_read = False
    for l in open(lpfile):
        if "PROCESSING OF IMAGES" in l:
            flag_read = False
            l = l.strip()
            first, last = map(lambda x:int(x.strip()), l[l.index("PROCESSING OF IMAGES")+len("PROCESSING OF IMAGES"):].split("..."))
            if first <= frame <= last:
                flag_read = True

        if flag_read:
            for key, s in keys.items():
                if s in l:
                    l = l.strip()
                    val = map(lambda x:float(x.strip()), l[l.index(s)+len(s):].split())
                    data[key] = val

    beam = data["beam direction"]
    rotaxis = data["rotation axis"]
    distance = data["distance"][0]
    orgx, orgy = data["beam center"]
    spacegroup = data["spacegroup"][0]
    a, b, c, alpha, beta, gamma = data["cell"]
    aaxis = data["a axis"]
    baxis = data["b axis"]
    caxis = data["c axis"]

    xp = XPARM("XPARM.XDS")

    xp.rotation_axis = numpy.array(rotaxis)
    xp.incident_beam = numpy.array(beam)
    xp.spacegroup = spacegroup
    xp.unit_cell = numpy.array((a, b, c, alpha, beta, gamma))
    xp.a_axis = numpy.array(aaxis)
    xp.b_axis = numpy.array(baxis)
    xp.c_axis = numpy.array(caxis)
    xp.origin = numpy.array((orgx, orgy))
    xp.distance = distance

    return xp.xparm_str()
# get_xparm_from_integrate_lp()
