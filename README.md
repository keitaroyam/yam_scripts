# yam_scripts
Shell or python (v2.7) scripts for crystallographic works (mostly for XDS)

Most of the scripts require [phenix](http://www.phenix-online.org/) installed (uses phenix.python)

## Contents
### Executables

command                  | description
------------------------ | ----------------------------
adxvx                    | adxv launcher for XDS-cbf files (read wavelength etc from XDS.INP in the same directory)
compare_XDS.INP          | to see difference between two XDS.INP files
copy_free_R_flag.py      | copy (and extend) R-free flag from other reflection file
create_free_R_flag.py    | create R-free flag (with ccp4 or phenix/cns style)
xds2mtz.py               | convert XDS file (XDS_ASCII.HKL or XSCALE output) to MTZ which contains all data (FP & IMEAN or F(+/-) & FP/DANO/ISYM & I(+/-))
xds_beamcenter_search.py | simple script to perform grid search of beam center coordinates in IDXREF
xds_merge_framecbf.py    | merge two FRAME.cbf of same image but with different crystal orientation to see two predictions with different colors (in adxv)
xds_plot_integrate.py    | make a plot from INTEGRATE.LP, which can be opened with loggraph in CCP4
xds_predict_mitai.py     | make a prediction file (FRAME.cbf) for arbitrary frame number
xds_profile_mitai.py     | make a 3D map file of 3D profiles from INTEGRATE.LP, which can be seen with PyMOL
xscale_simple.py         | to easily run XSCALE to merge XDS files
load_xparm.py            | PyMOL plugin to see crystal orientation visually from XPARM.XDS

### Libraries the executables depend on:

* make_adx.py
* mtzutil.py
* util.py
* xds.py
* xds_files.py
* xparm.py
