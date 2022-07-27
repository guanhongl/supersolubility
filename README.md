# supersolubility

If you wish to run the program locally, you can follow the directions below. 

## configuring, building, and testing FLINT

1) FLINT requires GMP and MPFR. download GMP (version 5.1.1 or later) and MPFR (3.0.0 or later). e.g. for ubuntu
```
apt-get install -y libgmp-dev libmpfr-dev
```
2) download [FLINT](http://www.flintlib.org/downloads.html)
3) untar it on your system
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

## running the program

1) download Makefile and ss.c
2) in the main directory, type
```
make
```
3) to compute the # of supersoluble groups of order n for n = 10^k, type
```
./ss k 1
```
4) to compute the # of supersoluble groups of order n to m, type
```
./ss n m 1
```
5) to control the maximum number of threads, num_threads, type
```
./ss k num_threads
```

## gotchas

If it happens that GMP and MPFR are not in a standard location on your system (e.g. not in /usr/include/ and /usr/lib/), you need to tell the configure script where they are with the options --with-gmp=/path/to/gmp or --with-mpir=/path/to/mpir and --with-mpfr=/path/to/mpfr, e.g.
```
./configure --with-gmp=/path/to/gmp --with-mpfr=/path/to/mpfr
```
If it happens that FLINT is not in a standard location on your system, you need to tell the compiler and linker (e.g. gcc) where they are. Edit Makefile, e.g.
```
LFLAGS=-lflint -lmpfr -lgmp -lpthread -L/path/to/flint/lib
CFLAGS=-Wall -O3 -I/path/to/flint/include
```

## future work

Rewrite the program using shared-memory parallel programming models like OpenMPI to utilize multiple nodes.
