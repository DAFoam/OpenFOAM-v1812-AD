OpenFOAM-v1812-AD
=================

This repository contains the OpenFOAM-v1812 source codes differentiated by automatic differentiation (AD) in forward and reverse modes.

Installation
------------

The installation of OpenFOAM-v1812-AD is similar to that of OpenFOAM-v1812. One needs to first install all prerequisites, source the OpenFOAM-v1812-AD/etc/bashrc file, and then run `./Allwmake.sh`.

The default build will be for forward mode AD. To compile reverse mode AD, change `WM_CODI_AD_MODE` to `reverse` in OpenFOAM-v1812-AD/etc/bashrc, source it, and rebuild.

NOTE: OpenFOAM-v1812-AD only differentiates necessary libraries for computing partial derivatives and matrix-vector products for [DAFoam](https://dafoam.github.io), it has NOT differentiated the entire OpenFOAM code yet. In other words, some functionalities are still missing (e.g. combustion models).

Acknowledgement
---------------

OpenFOAM-v1812-AD uses the AD tools [CoDiPack](https://github.com/scicompkl/codipack) and [MeDiPack](https://github.com/scicompkl/medipack), developed by Dr. Nicolas Gauger's group at TU Kaiserslautern. In addition, the differentiation is inspired by [discreteAdjointOpenFOAM](https://www.stce.rwth-aachen.de/research/software/discreteadjointopenfoam), developed by Dr. Markus Towara et al. from RWTH Aachen University.

License
-------

Distributed using the GNU General Public License (GPL), version 3; see the COPYING file for details.