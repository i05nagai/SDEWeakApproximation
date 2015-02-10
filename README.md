SDE\_WA\_MLMC
====

Overview

This library is an implementation of a Multi-Level Monte-Carlo(MLMC) Algorithm proposed by Mike Giles in a [paper(pdf file)](https://people.maths.ox.ac.uk/gilesm/files/OPRE_2008.pdf).

For practical use, this library includes two other libraries:

* [sde\_wa.tar.gz](https://sites.google.com/site/marikoninomiya/system/app/pages/admin/revisions?wuid=wuid:gx:3d39bf39c9952da6) 
developed by [Mariko Ninomiya](https://sites.google.com/site/marikoninomiya/mariko-ninomiya) 
and [Syoiti Ninomiya](https://sites.google.com/site/craft-titech-ac-jp-ninomiya-craft-e-teams-edition-backup/).
    * sde\_wa is a library for Weak Approximations of Stochastic Differential Equations.
	* LGPL v2.1
* [mt19937ar.tgz](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.tgz) developed by Takuji Nishimura and Makoto Matsumoto.
    * mt19937ar is a Mersenne Twister pseudorandom number generator with period 2^19937-1 with improved initialization scheme.
	* Modified BSD Licence

## Description

The MLMC algorithm is heuristic implementation of the MLMC method, which is a technique to reduce computational complexity of weak approximations of Stochastic Differential Equations(SDEs).
For more detailed description of MLMC method, please refer to a [site](https://people.maths.ox.ac.uk/gilesm/mlmc_community.html).

## Usage

To use the library, please see a repository, samples of MLMC estimation by Kusuoka Schemes and Euler Maruyama Scheme.

## Install

To install the library, please follow three steps:

1. Git clone this repository.

    `git clone https://github.com/i05nagai/SDEWeakApproximation.git`

2. Create Makefiles.

    `make Makefiles`

3. Make the library.

    `make`

## Licence

This software is released under the LGPL v2.1, see LICENSE.

## Author

[i05nagai](https://github.com/i05nagai)

