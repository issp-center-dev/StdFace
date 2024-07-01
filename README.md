# StdFace

An input file generator for
[HPhi](https://github.com/issp-center-dev/HPhi),
[mVMC](https://github.com/issp-center-dev/mVMC),
[UHF](https://github.com/issp-center-dev/UHF-dev), and
[H-wave](https://github.com/issp-center-dev/H-wave).

## Requirement
C compiler (intel, Fujitsu, GNU, etc. )  

## Get StdFace

### With Git 

``` bash
$ git clone https://github.com/issp-center-dev/StdFace
```

### Without Git (supported from ver.1.0)

You can download StdFace from [release page](https://github.com/issp-center-dev/StdFace/releases).

## Install StdFace

``` bash
$ cmake -B build [options]
$ cmake --build build
$ cmake --install install
```

Options need to be specified which programs are to be built:

- ``-DHPHI=ON``

    HPhi mode is enabled. hphi_dry.out and libStdFace_hphi.a will be generated.

- ``-DMVMC=ON``

    mVMC mode is enabled. mvmc_dry.out and libStdFace_mvmc.a will be generated.

- ``-DUHF=ON``

    UHF mode is enabled. uhf_dry.out and libStdFace_uhf.a will be generated.

- ``-DHWAVE=ON``

    H-wave mode is enabled. hwave_dry.out and libStdFace_hwave.a will be generated.

One or more options may be set simultaneously. Defaults are OFF.


## Usage

xxx_dry.out (xxx = hphi, mvmc, uhf, hwave) can be done with the following command with input.in as input file.

``` bash
$ xxx_dry.out input.in
```

After execution, input files for executing each tool will be generated in the executed directory.
Descriptions of StdFace input files can be found in the manuals for each tool:

[HPhi](https://github.com/issp-center-dev/HPhi),
[mVMC](https://github.com/issp-center-dev/mVMC),
[UHF](https://github.com/issp-center-dev/UHF-dev), and
[H-wave](https://github.com/issp-center-dev/H-wave).

## Licence

The distribution of the program package and the source codes for StdFace follow GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)). 


## Author
Kazuyoshi Yoshimi, Mitsuaki Kawamura, Kota Ido, Yuichi Motoyama, and Tatsumi Aoyama.
