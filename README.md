# A Coulomb attenuating method range separated functional implemented quantum esspresso

This is a modified verison of quantum espresso(base version is the devel branch[https://gitlab.com/QEF/q-e/-/tree/089dc4e7cfe0f2f6fdbec2fc07857fca4e31c76e]), in which the coulomb attenuating method range separated functionals can be used in pw calculation by obataining parameters from libxc[https://gitlab.com/libxc/libxc] or setting parameters in input file. Besides, rvv10 paramters can also be obtained from libxc.


## 1. Coulomb attenuating method range separated functional

A additional input item has been added, "exx_fraction_lr", which is used to set the fraction of long range part of exact HF exchange potential. Combined with input item of "exx_fraction", the fraction of exact HF exchange potential, range separated functional calculation can be conducted. The exact HF exchange potential can be expressed as 

$$V_{exx} = exx\\_fraction\cdot\frac{1}{|r-r'|} + exx\\_fraction\\_lr\cdot\frac{efr(\mu |r-r'|)}{|r-r'|}$$

By using the libxc, the parameters can be obatained automatically. There is no need to set those parameters in input files. 

Such as the example input in the test[https://github.com/dcccc/quantum-espresso/blob/master/test/N2_camb3lyp/pw.in]:
```
    ! it is the cam-b3lyp functional 
    input_dft='xc-000i-000i-433l-000i-000i-000i',
```

In the output file[https://github.com/dcccc/quantum-espresso/blob/master/test/N2_camb3lyp/pw.out] of pw calculaiton, you can check the parameters from libxc used 

```
......
     Current dimensions of program PWSCF are:
     Max number of different atomic species (ntypx) = 10
     Max number of k-points (npk) =  40000
     Max angular momentum in pseudopotentials (lmaxx) =  4
     screening_parameter from libxc    0.3300000
     exx_fraction from libxc   0.1900000
     exx_fraction_lr from libxc   0.4600000

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation= XC-000I-000I-433L-000I-000I-000I
                           (   0   0   0 433   0   0   0)
     EXX-fraction              =        0.19
     Any further DFT definition will be discarded
     Please, verify this is what you really want
......
```

Besides, the parametes("exx_fraction", "exx_fraction_lr", "screening_parameter") can also modified for a better result by tuning the fractions and screening length.

For example input in the test[https://github.com/dcccc/quantum-espresso/blob/master/test/N2_camb3lyp/pw_m.in]:
```
    input_dft='xc-000i-000i-434l-000i-000i-000i', exx_fraction=0.19,  exx_fraction_lr=0.46,screening_parameter=0.33,
```
The xc functional in "input_dft" is the TUNED_CAM_B3LYP[https://doi.org/10.1016/j.jphotochem.2012.03.003], not the CAM_B3LYP[https://doi.org/10.1016/j.cplett.2004.06.011]. However, by setting the "exx_fraction", "exx_fraction_lr", "screening_parameter" parameters in input files(those parameters are also updated in libxc at the begaining of calculaitons), the functional is actually bacomes the CAM_B3LYP. cheking the output file[https://github.com/dcccc/quantum-espresso/blob/master/test/N2_camb3lyp/pw_m.in]

```
......
     screening_parameter from libxc    0.1500000
     exx_fraction from libxc   0.0799000
     exx_fraction_lr from libxc   0.9201000
......
     exx_fraction to libxc    0.6500000
     exx_fraction_sr to libxc   -0.4600000
     screening_parameter to libxc    0.3300000
......
```

we can see that "exx_fraction", "exx_fraction_lr", "screening_parameter" parameters are updated, and the final energies is identical to the energy with CAM_B3LYP functional without paramters setting in input file. 

As to the name convention used here is different from libxc. So the values updated into libxc are not the same with values setted in input files.


## 2. wb97mrv, wb97xrv and b97mrv functionals

The wb97mv(id=531), b97mv(id=254) and wb97xv(id=466) funtionals is not possiable at present with  quantum espresso, because the vv10[https://doi.org/10.1063/1.3521275] nonlocal functional is not implemented. But, as the rvv10[https://doi.org/10.1103/PhysRevB.87.041108] can be used, a rvv10 substitution version is possiable. The rvv10 paramters in those functionals had been optimized by Narbe Mardirossian et al.[https://pubs.acs.org/doi/10.1021/acs.jpclett.6b02527]. And there are available in a modified version of libxc[https://github.com/dcccc/libxc] (wb97mrv(id=768), b97mrv(id=769), wb97xrv(id=767)).

Taking the case in the test as example[https://github.com/dcccc/quantum-espresso/blob/master/test/N2_wb97mv/pbe0_wb97mrv.in]

'''
    ! the wb97mrv
    input_dft='xc-000i-000i-000i-000i-768l-000i-vv10',
'''

In the output file, the parameters of can be checked

```
......
     rvv10 b value from libxc   6.2000000
     rvv10_c value from libxc   0.0093000
     screening_parameter from libxc    0.3000000
     exx_fraction from libxc   0.1500000
     exx_fraction_lr from libxc   0.8500000

     IMPORTANT: XC functional enforced from input :
     Exchange-correlation= XC-000I-000I-000I-000I-768L-000I-VV10
                           (   0   0   0   0  26   0 768)
     EXX-fraction              =        0.15
......
```


## NOTICE


1.  "vcut_spherical" method to treat the Coulomb potential divergencies at small q vectors for range separate functional is added. 'gygi-baldereschi' and 'vcut_ws' method are unmodified. You can refer to discription in vasp manual[https://doi.org/10.1103/PhysRevB.77.193110] or original article[https://doi.org/10.1103/PhysRevB.77.193110] for details.

2. In the calculations, "exx_fraction" can not be setted to a zero value, which will cause a numerial error in the calculations. If only long range part of exact exchange potential is needed, a small value of "exx_fraction", like 0.0001, maybe reasonable. In that case, a additional negligible exact exchange potential may have little influence on the final result.

3. Threr is **NO GUARANTY**. Before production calcuations you should always do test calcultions and make sure this is really what you want. A good calculaiton result comes from a reasonable model, and the calculation is always a rubbish in rubbish out style work.


Following is the original readme content

=====================================================================================================================================

![q-e-logo](logo.jpg)

This is the distribution of the Quantum ESPRESSO suite of codes (ESPRESSO:
opEn-Source Package for Research in Electronic Structure, Simulation, and
Optimization)

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## USAGE
Quick installation instructions for CPU-based machines. For GPU execution, see
file [README_GPU.md](README_GPU.md). Go to the directory where this file is. 

Using "make"
(`[]` means "optional"):
```
./configure [options]
make all
```
"make" alone prints a list of acceptable targets. Optionally,
`make -jN` runs parallel compilation on `N` processors.
Link to binaries are found in bin/.

Using "CMake" (v.3.14 or later):

```
mkdir ./build
cd ./build
cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc [-DCMAKE_INSTALL_PREFIX=/path/to/install] ..
make [-jN]
[make install]
```
Although CMake has the capability to guess compilers, it is strongly recommended to specify
the intended compilers or MPI compiler wrappers as `CMAKE_Fortran_COMPILER` and `CMAKE_C_COMPILER`.
"make" builds all targets. Link to binaries are found in build/bin.
If `make install` is invoked, directory `CMAKE_INSTALL_PREFIX`
is prepended onto all install directories.

For more information, see the general documentation in directory Doc/, 
package-specific documentation in \*/Doc/, and the web site 
http://www.quantum-espresso.org/. Technical documentation for users and
developers 
can be found on [Wiki page on gitlab](https://gitlab.com/QEF/q-e/-/wikis/home).

## PACKAGES

- PWscf: structural optimisation and molecular dynamics on the electronic ground state, with self-consistent solution of DFT equations;
- CP: Car-Parrinello molecular dynamics;
- PHonon: vibrational and dielectric properties from DFPT (Density-Functional Perturbation Theory);
- TD-DFPT: spectra from Time-dependent DFPT;
- HP: calculation of Hubbard parameters from DFPT;
- EPW: calculation of electron-phonon coefficients, carrier transport, phonon-limited superconductivity and phonon-assisted optical processes;
- PWCOND: ballistic transport;
- XSpectra: calculation of X-ray absorption spectra;
- PWneb: reaction pathways and transition states with the Nudged Elastic Band method;
- GWL: many-body perturbation theory in the GW approach using ultra-localised Wannier functions and Lanczos chains;
- QEHeat: energy current in insulators for thermal transport calculations in DFT.
- KCW: Koopmans-compliant functionals in a Wannier representation

## Modular libraries
The following libraries have been isolated and partially encapsulated in view of their release for usage in other codes as well:

- UtilXlib: performing basic MPI handling, error handling, timing handling.
- FFTXlib: parallel (MPI and OpenMP) distributed three-dimensional FFTs, performing also load-balanced distribution of data (plane waves, G-vectors and real-space grids) across processors.
- LAXlib: parallel distributed dense-matrix diagonalization, using ELPA, SCALapack, or a custom algorithm.
- KS Solvers: parallel iterative diagonalization for the Kohn-Sham Hamiltonian (represented as an operator),using block Davidson and band-by-band or block Conjugate-Gradient algorithms.
- LRlib: performs a variety of tasks connected with (time-dependent) DFPT, to be used also in connection with Many-Body Perturbation Theory.
- upflib: pseudopotential-related code.
- devXlib: low-level utilities for GPU execution

## Contributing
Quantum ESPRESSO is an open project: contributions are welcome.
Read the [Contribution Guidelines](CONTRIBUTING.md) to see how you
can contribute.

## LICENSE

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.
