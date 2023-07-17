# StdFace

An input files generator for
[HPhi](https://github.com/issp-center-dev/HPhi),
[mVMC](https://github.com/issp-center-dev/mVMC),
[UHF](https://github.com/issp-center-dev/UHF-dev), and
[H-wave](https://github.com/issp-center-dev/H-wave).

## Requirement
C compiler (intel, Fujitsu, GNU, etc. )  

## Get StdFace

### With Git 

```
bash
$ git clone https://github.com/issp-center-dev/StdFace
```


### Without Git (supported from ver.1.0)

You can download StdFace from a [release note](https://github.com/issp-center-dev/StdFace/releases).

## Install StdFace

``` bash
$ cmake -B build [options]
$ cmake --build build
$ cmake --install install
```

If you define -DUHF=OFF as [options], uhf_dry.out will not be made.
If you define -DMVMC=OFF as [options], mvmc_dry.out will not be made.
If you define -DHPHI=OFF as [options], hphi_dry.out will not be made.
If you define -DHWAVE=OFF as [options], hawave_dry.out will not be made.

## Licence

The distribution of the program package and the source codes for StdFace follow GNU General Public License version 3 ([GPL v3](http://www.gnu.org/licenses/gpl-3.0.en.html)). 


## Author
Kazuyoshi Yoshimi, Mitsuaki Kawamura.

## Manual

under construction
