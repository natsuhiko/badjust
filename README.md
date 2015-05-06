# badjust
Multiple testing correction with beta distribution

## How to install

  CFLAGS="-I/usr/include/gsl\ -I/usr/include"
  LDFLAGS="-L/usr/lib"
  LIBR="-lgsl\ -llapack\ -lm"
  MAKEFLAGS="CFLAGS=$CFLAGS LDFLAGS=$LDFLAGS LIBR=$LIBR" R CMD INSTALL master --no-multiarch
