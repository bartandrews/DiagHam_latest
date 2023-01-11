DiagHam_latest
==============

DiagHam fork to accompany the paper "Self-similarity of spectral response functions for fractional quantum Hall states" `[Andrews22] <https://arxiv.org/abs/2201.04704>`__.

- DiagHam wiki: https://nick-ux.org/diagham/index.php/Main_Page
- DiagHam website: http://www.phys.ens.fr/~regnault/diagham/

Getting started
---------------

These installation instructions are in addition to those listed on the `DiagHam wiki <https://nick-ux.org/diagham/index.php/Install>`__.

1) Intel libraries need to be installed for optimal performance. At the time of writing, you need the Intel oneAPI Base Toolkit for the C/C++ compiler and MKL library and the Intel oneAPI HPC Toolkit for the Fortran compiler and MPI library. After the installation, you can remove the installer if desired, which is in e.g. ``/tmp/root/`` or ``~/Downloads/``. You can also add the following two lines to your ``~/.bashrc``: ``source /opt/intel/oneapi/compiler/latest/env/vars.sh``; ``source /opt/intel/oneapi/mkl/latest/env/vars.sh``. This saves some time when starting a shell, compared to sourcing the entire ``/opt/intel/oneapi/setvars.sh``.

2. The configure command to automatically use the latest instruction set (``-xHost``) is:

``../configure --enable-fqhe --enable-fti --enable-lapack --enable-gmp --enable-lapack-only --with-lapack-libs="" --with-blas-libs="-mkl" CC=icc CXX=icpc --enable-debug CFLAGS="-O3 -xHOST" CXXFLAGS="-O3 -xHOST"``

...or for using both AVX and AVX2 instruction sets...

``../configure --enable-fqhe --enable-fti --enable-lapack --enable-gmp --enable-lapack-only --with-lapack-libs="" --with-blas-libs="-mkl" CC=icc CXX=icpc --enable-debug CFLAGS="-O3 -xAVX -axCORE-AVX2" CXXFLAGS="-O3 -xAVX -axCORE-AVX2"``

NB: BLAS will always be called when LAPACK is called, and it contains LAPACK, so no need to duplicate mkl flags. Since MKL is provided by both BLAS and LAPACK, it’s sufficient to give one or the other – but both options are required so one can also cope with separate libraries.

A guide to Intel compiler flags can be found here: https://www.bu.edu/tech/support/research/software-and-programming/programming/compilers/intel-compiler-flags/

References
----------

`[Andrews22] <https://arxiv.org/abs/2201.04704>`__ "Self-similarity of spectral response functions for fractional quantum Hall states" by Bartholomew Andrews and Gunnar Möller, arXiv:2201.04704 [cond-mat.str-el].
