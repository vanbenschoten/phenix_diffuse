#Introduction

phenix.diffuse uses Guinier’s equation to calculate diffuse scattering intensity maps from Protein Data Bank (PDB)-formatted structural ensembles. This program is released as a part of the PHENIX X-ray crystallography software suite (phenix-online.org) distributions #1894 and later. 

For further documentation see Van Benschoten et.al, 2015 (Acta Crystallographica, Section D).

**NOTE:** This code is written within the framework of the Computational Crystallography Toolbox python library (http://cctbx.sourceforge.net/) and thus will not work with standard python distributions.



#Contents

diffuse.py:                  Final code checked in to Phenix distribution #1894

diffuse_fractional.py:       Code for unit cell expansion calculations discussed in Van Benschoten et.al.

diffuse_p1.py:               Positive control for diffuse_fractional.py



#Instructions

###Command line usage:###
phenix.diffuse pdb=test.pdb probabilities=0.3,0.3,0.3,0.1 resolution=4.0 prefix=my_map


###Program arguments###
*pdb*: input PDB file. At least one CRYST1 symmetry description is required.

*probabilities*: the weighted probability for each model in the corresponding PDB file. The values listed in this option will be assigned to the PDB models starting with MODEL 0. Alternatively, leaving out the probabilities option will lead to equal weights for each model in the PDB file.

*resolution*: the desired d_min for the output MTZ map

*prefix*: the desired MTZ file name


###Output###
MTZ file containing diffuse scattering intensity values (I) at each Bragg peak



#Authors
Andrew Van Benschoten (andrew.vanbenschoten@ucsf.edu)

Pavel Afonine (pafonine@lbl.gov)

Alexandre Urzhumtsev (sacha@igbmr.fr)
