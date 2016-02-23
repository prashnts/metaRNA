
libpymirandadir = $(pkgpythondir)

pyexec_LTLIBRARIES = libpymiranda.la
libpymiranda_la_LDFLAGS = -module -avoid-version -export-symbols-regex PyInit_libpymiranda
libpymiranda_la_CFLAGS = $(PYTHON_CFLAGS) -I/usr/local/include/ViennaRNA
libpymiranda_la_LIBADD = $(PYTHON_LIBS) -L/usr/local/lib -lRNA -lm
libpymiranda_la_SOURCES = pymiranda.c output.c scan.c swat.c thermo.c utils.c ExpString.c
libpymiranda_la_HEADERS = miranda.h
libpymiranda_ladir = .dir


EXTRA_DIST = libpymiranda.h
