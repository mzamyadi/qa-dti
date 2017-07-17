qa-dti
------

Spherical phantom quality assurance procedures for DTI sequences, described in

> A novel DTI-QA tool: automated metric extraction exploiting the sphericity of an agar filled phantom. 2017. Sofia Chavez, Joseph Viviano, Mojdeh Zamyadi, Peter B Kingsley, Peter Kochunov, Stephen Strothers, Aristotle Voineskos. Magnetic Resonance Imaging.

All code written by Sofia Chavez 2014-2017. Packaged and maintained by Joseph D Viviano.

**quickstart**

- requires the installation of MATLAB (core), python 2.7, and FSL.
- to install, add this folder to your path: `export PATH="${PATH}:/path/to/qa-dti"`.
- `qa-dti --help` for usage.
- `qa-dti --verbose` for debugging outputs.

**details**

- the full pipeline is found in `matlab/analyze_dti_phantom.m`
- the subprocesses are found alongside this pipeline in `matlab/*.m`. They can be run individually if desired.
- not tested with Octave (reach out if you want to help!)

