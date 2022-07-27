# supersolubility

If you wish to run the program locally, you can follow the directions below. 

## configuring, building, and testing FLINT

3) FLINT requires GMP and MPFR. download GMP (version 5.1.1 or later) and MPFR (3.0.0 or later). e.g. for ubuntu
```
apt-get install -y libgmp-dev libmpfr-dev
```
1) download [FLINT](http://www.flintlib.org/downloads.html)
2) untar it on your system
```
tar -xvzf flint-2.9.0.tar.gz
```
4) in the main directory of the FLINT directory tree, type
```
./configure
```
5) once FLINT is configured, type
```
make
make check
make install
```
6) FLINT comes with example programs to demonstrate current and future FLINT features. To build the example programs, type
```
make examples
```
7) you may need to set your LD_LIBRARY_PATH or equivalent for the FLINT, MPIR and MPFR libraries. e.g. for ubuntu
```
export LD_LIBRARY_PATH="/usr/local/lib"
```

For further information, [visit](http://www.flintlib.org/doc/building.html).
