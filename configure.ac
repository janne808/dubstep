#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([dubstep], [0.003], [janne808@radiofreerobotron.net])
AC_CONFIG_SRCDIR([src/dubstep.c])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [main])
AC_CHECK_LIB([rt], [main])
AC_CHECK_LIB([pthread], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([clock_gettime floor pow sqrt])

# Package options
AC_ARG_ENABLE([sdl], AS_HELP_STRING([--enable-sdl],[Enable SDL]))

AS_IF([test "x$enable_sdl" = "xyes"],[
	    AC_DEFINE([ENABLE_SDL], [1], [Define to enable SDL])
	    
	    # Check for SDL and GL libs
	    SDL_VERSION=1.2.0
	    AM_PATH_SDL($SDL_VERSION,
		:,
   	    	AC_MSG_ERROR([*** SDL version $SDL_VERSION not found!])
	    )
	    CFLAGS="$CFLAGS $SDL_CFLAGS"
	    LIBS="$LIBS $SDL_LIBS"
	    
	    AC_CHECK_LIB([GL], [main])
])

AC_ARG_ENABLE([tiffio], AS_HELP_STRING([--enable-tiffio],[Enable TIFFIO]))

AS_IF([test "x$enable_tiffio" = "xyes"],[
	    AC_DEFINE([ENABLE_TIFFIO], [1], [Define to enable TIFFIO])	    
	    AC_CHECK_LIB([tiff], [main])
])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
