#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(wcd-express,0.6.3,Scott.Hazelhurst@wits.ac.za)
AM_INIT_AUTOMAKE()
AM_CONFIG_HEADER(src/config.h)
AC_CONFIG_SRCDIR(src/wcd.c)
AC_CONFIG_FILES([src/Makefile doc/Makefile])
# Checks for programs.


echo "Checking MPI"
AC_ARG_ENABLE(MPI,[   --enable-MPI     enable MPI support], [CC=mpicc],[])
AC_ARG_ENABLE(mpi,[   --enable-mpi     enable MPI support], [CC=mpicc],[])



echo "Checking if Pthreads asked for"
AC_ARG_ENABLE(PTHREADS,
             [   --enable-pthreads     enable pthreads support],\
             [PTHREADS=true],[])
AC_ARG_ENABLE(pthreads,
             [   --enable-pthreads     enable pthreads support],\
             [PTHREADS=true],[])


echo "Checking whether pthread library is there"
AC_CHECK_LIB(pthread,pthread_create,[HAVE_LIBPTHREAD=1])


echo "Checking whether auxiliary info wanted"
AC_ARG_ENABLE(NOAUXINFO,[   --enable-noauxinfo    don't store auxinfo], [NOAUXINFO=true],[])
AC_ARG_ENABLE(noauxinfo,[   --enable-noauxinfo    dont' store auxinfo], [NOAUXINFO==true],[])


echo "Checked whether pthread library is there"
case ${PTHREADS} in
true) echo "You have asked for Pthreads"
case ${HAVE_LIBPTHREAD} in
1) echo "have it so wcd-express should Vroom.............."
;;
*) echo "but I CAN'T find it. Please either install or get path right, or don't enable  it"
echo "Sorry!"
exit -1
;;
esac
;;
*) echo " You have not asked for pthreads"

;;
esac


echo "Checking compiler"

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_MAKE_SET


echo "Now doing other checks"


# Check whether big endian

AC_C_BIGENDIAN

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_HEADER_TIME
AC_HEADER_STDBOOL


# Checks for library functions.
AC_FUNC_MALLOC
AC_HEADER_STDC
AC_CHECK_FUNCS([memset strcasecmp getopt_long])
AM_CONDITIONAL(COMPILER, test ${CC}=gcc)
AM_CONDITIONAL(MPI, test x${CC} = xmpicc)
AM_CONDITIONAL(PTHREADS, test x${PTHREADS} = xtrue)
AM_CONDITIONAL(NOAUXINFO, test x${NOAUXINFO} = xtrue)
AM_CONDITIONAL(LONGOPT, test -n ${HAVE_GETOPT_LONG})
AM_CONDITIONAL(LONGOPT2, test -f /usr/include/getopt.h)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
