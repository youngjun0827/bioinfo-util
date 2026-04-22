AC_PREREQ([2.69])
AC_INIT([bioinfo-util], [1.0])
AM_INIT_AUTOMAKE([1.11 foreign -Wall])
AC_CONFIG_SRCDIR([src/fastq_stats.cpp])

AC_PROG_CXX
AC_LANG([C++])

AC_CHECK_HEADER([zlib.h], [],
    [AC_MSG_ERROR([zlib headers not found. Install zlib-dev or zlib-devel.])])
AC_CHECK_LIB([z], [gzopen],
    [LIBS="-lz $LIBS"],
    [AC_MSG_ERROR([zlib library not found. Install zlib-dev or zlib-devel.])])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
