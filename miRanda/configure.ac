
AC_PREREQ(2.60)
AC_INIT(libpymiranda, 4.0, prashant {at} one.ducic.ac.in)
AC_CONFIG_HEADER([config.h])
AC_SUBST(VERSION)
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)
AM_INIT_AUTOMAKE
AC_CANONICAL_HOST

LT_INIT
AC_PROG_CXX

#AM_PATH_PYTHON(3.5.0)
PYTHON_LOC=`python3-config --prefix`
PYTHON_LIBS=`python3-config --ldflags`
PYTHON_CFLAGS=`python3-config --cflag`
AC_SUBST([PYTHON_LIBS])
AC_SUBST([PYTHON_CFLAGS])

# Checks for programs.
AC_PROG_CC
if test "$GCC" = yes; then
   AC_DEFINE(UNUSED, __attribute__ ((unused)), avoid warning about unused variables)
else
   AC_DEFINE(UNUSED,, only for gcc)
fi
AC_PROG_INSTALL
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memset sqrt strchr strdup strrchr strstr])

AC_CONFIG_FILES([Makefile
                 man/Makefile
                 src/Makefile])
AC_OUTPUT
