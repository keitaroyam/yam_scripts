"""
mtzutil.py

(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import sys
import os
import struct
import math
import subprocess

##
# MTZ COLUMN TYPES
#H	 index h,k,l
#J	 intensity
#F	 structure amplitude, F
#D	 anomalous difference
#Q	 standard deviation of J,F,D or other (but see L and M below)
#G	 structure amplitude associated with one member of an hkl -h-k-l pair, F(+) or F(-)
#L	 standard deviation of a column of type G
#K	 intensity associated with one member of an hkl -h-k-l pair, I(+) or I(-)
#M	 standard deviation of a column of type K
#E	 structure amplitude divided by symmetry factor ("epsilon"). Normally scaled as well to give normalised structure factor
#P	 phase angle in degrees
#W	 weight (of some sort)
#A	 phase probability coefficients (Hendrickson/Lattman)
#B	 BATCH number
#Y	 M/ISYM, packed partial/reject flag and symmetry number
#I	 any other integer
#R	 any other real

class Error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# class Error

class MtzFile:

    def __init__(self, mtzfile):
        self._filename = mtzfile

        self.mtzheader = ""
        self.get_header()
        #print "\n".join(self.mtzheader)
        self.get_datasets()

        self.columns_to_be_used = {"Fini":None, "SIGFini":None,
                                   "Fobs":None, "SIGFobs":None,
                                   "FREE":None, "PHI":None, "FOM":None,
                                   "HLA":None, "HLB":None, "HLC":None, "HLD":None}

    # __init__()

    def get_filename(self):
        return self._filename
    # get_name()

    def get_header(self):
        
        def doSwap(m):
            s, n = 0, 4
            f = (m[0]>>4)& 0x0f
            if sys.byteorder != 'little': n = 1
            if f != 0 and f != n: s = 1
            return s == 1
        # doSwap()

        int_size = struct.calcsize("i")

        if not os.path.isfile(self._filename):
            raise Error, "The file does not exist: %s" % self._filename

        try:
            mtz = open(self._filename, 'rb')
            
            if mtz.read(4) != "MTZ ":
                raise Error, "This file is not MTZ format: " + self._filename

            mtz.seek(8, 0)
            to_swap = doSwap(struct.unpack('4B', mtz.read(4)))
            
            mtz.seek(4, 0)
    
            if to_swap:
                byteorder = "<>"[int(sys.byteorder == 'little')]
                hdrst, = struct.unpack(byteorder+'i', mtz.read(int_size))
            else:
                hdrst, = struct.unpack('i', mtz.read(int_size))

            mtz.seek(4*(hdrst-1), 0)
            header = mtz.read()

        finally:
            mtz.close()
        
        self.mtzheader = []

        while header:
            self.mtzheader.append(header[:80])
            header = header[80:]
         
    # get_header()

    def get_datasets(self):
        flag_read = False
        self.datasets = {}

        for line in self.mtzheader:
            if line.startswith("NDIF"):
                flag_read = True
                #self.datasets[int(line.split()[1])] = {}
                continue

            if flag_read:
                if line.startswith("END"):
                    flag_read = False
                else:
                    item = line.split()
                    self.datasets.setdefault(int(item[1]), {})[item[0]] = " ".join(item[2:])

    # get_datasets()

    def get_fullname(self, lab):
        assert self.datasets != []

        for line in self.mtzheader:
            item = line.split()
            if item[0] == "COLUMN":
                col, coltype, colmin, colmax, colid = item[1:]
                if lab == col:
                    dataset = self.datasets[int(colid)]
                    return "/".join(("",dataset["CRYSTAL"], dataset["DATASET"], col))

        return None

    # get_fullname()

    def get_column(self):
        column = {}

        if not self.mtzheader:
            return None

        for line in self.mtzheader:
            item = line.split()
            if item[0] == "COLUMN":
                column.setdefault(item[2],[]).append(item[1])

        return column

    # get_column()

    def get_column_minmax(self, col):
        for line in self.mtzheader:
            item = line.split()
            if item[0] == "COLUMN" and item[1] == col:
                minmax = map(lambda x:float(x), item[3:5])
                return min(minmax), max(minmax)

        return None, None
    # get_column_minmax()

    def get_ndif(self, lab):
        if not self.mtzheader:
            return None

        for line in self.mtzheader:
            item = line.split()
            if item[0] == "COLUMN" and item[1] == lab:
                return int(item[5])
    # get_ndif()

    def get_spacegroup(self):

        for line in self.mtzheader:
            if line.startswith("SYMINF"): # Number of Symmetry operations, number of primitive operations, lattice type, SG number, SG name, PG name
                sgno = line.split()[4]
                sgname = line[line.find("'")+1:line.rfind("'")].replace("'","").strip()
                return [ int(sgno), sgname ]

        raise Error, "No spacegroup info."

    # get_spacegroup()

    def get_z(self):
        syminf1 = filter(lambda l: l.startswith("SYMINF"), self.mtzheader)[0]
        return int(syminf1.split()[1])

    # get_z()

    def get_resolution(self):

        for line in self.mtzheader:
            if line.startswith("RESO"):
                return map(lambda x: math.sqrt(1./float(x)), line.split()[1:])

        raise Error, "No resolution info."
            
    # get_resolution()

    def get_cell(self):
        return map(lambda x: float(x), self.get_cell_str().split())
    # get_cell()

    def get_dcell(self, lab):
        ##
        # lab: column name. such as FP, ..
        #
        ndif = self.get_ndif(lab)
        dcell = self.datasets[ndif]["DCELL"]
        return map(lambda x: float(x), dcell.split())
    # get_dcell()

    def get_cell_str(self):
        for line in self.mtzheader:
            if line.startswith("CELL"):
                return " ".join(line.split()[1:])
    # get_cell_str()

# class mtzutil

def mtzdmp(filename):
    
    cmd = os.path.join(os.environ['CBIN'],"mtzdump") + ' hklin "%s"'%os.path.normpath(filename)
    #to8_3filename(os.path.normpath(filename))
    
    p = subprocess.Popen( cmd, shell=True, cwd=os.getcwd(), stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
                         # close_fds=True )
    
    #p.stdin.write("HEAD\nGO\n")
    p.stdin.write("\n")
    p.stdin.close()
    output = p.stdout.read()
    if p.returncode is None:
        p.wait()
        
    #p.wait()
    return output

# mtzdmp()


if __name__ == "__main__":
	mtzin = sys.argv[1]
        u = MtzFile(mtzin)
        print "Cell parameters: ", u.get_cell()
        print "Space group: ", u.get_spacegroup()
        print "COLUMN: ", u.get_column()
        print "Resolution range: ", u.get_resolution()
        print
        print "header="
        print "\n".join(u.mtzheader)
        
        for c in reduce(lambda x,y:x+y, u.get_column().values()):
            print u.get_fullname(c)

        print u.datasets
