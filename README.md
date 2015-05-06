# badjust
Multiple testing correction with beta distribution

## How to install

You need to install GSL (http://www.gnu.org/software/gsl/gsl.html) in your environment first.  Then, please use the following example code to install the package:

    CFLAGS="-I/usr/include/gsl\ -I/usr/include"
    LDFLAGS="-L/usr/lib"
    LIBR="-lgsl\ -llapack\ -lm"
    MAKEFLAGS="CFLAGS=$CFLAGS LDFLAGS=$LDFLAGS LIBR=$LIBR" R CMD INSTALL badjust --no-multiarch

Note that, you may be required to install lapack (http://www.netlib.org/lapack/) as well.
