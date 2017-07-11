qa-dti
------

Spherical phantom quality assurance proceedures for DTI sequences.

All code written by Sofia Chavez 2015-2017. Packaged and maintained by Joseph D Viviano.

**quickstart**

- requires the installation of MATLAB (core), python 2.7, and FSL.
- not tested with Octave (reach out if you want to help!)
- to install, add this folder to your path: `export PATH="${PATH}:/path/to/qa-dti"`.
- `qa-dti` performs the procedures described in (Chavez et al. 2017; Ref). Please reference this paper if this QC procedure proves helpful.

**details**

- the full pipeline is found in `matlab/analyze_dti_phantom.m`
- the subprocesses are found alongside this pipeline in `matlab/*.m`. They can be run individually if desired.

