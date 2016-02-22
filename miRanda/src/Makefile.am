#--------------------------------------------------------------------------------
# miRanda- An miRNA target scanner, aims to predict mRNA targets for microRNAs,
# using dynamic-programming alignment and thermodynamics
#
# Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center, New York
#
# Distributed under the GNU Public License (GPL)
# See the files 'COPYING' and 'LICENSE' for details
#
# Authors: Anton Enright, Bino John, Chris Sander and Debora Marks
# Email: mirnatargets (at) cbio.mskcc.org - reaches all authors
#
# Written By: Anton Enright
#
# Please send bug reports to: miranda (at) cbio.mskcc.org
#
# If you use miRanda in your research please cite:
# Enright AJ, John B, Gaul U, Tuschl T, Sander C and Marks DS;
# (2003) Genome Biology; 5(1):R1.
#
# This software will be further developed under the open source model,
# coordinated by Anton Enright and Chris Sander:
# miranda (at) cbio.mskcc.org (reaches both).
#--------------------------------------------------------------------------------
#
# Copyright (C) (2003) Memorial Sloan-Kettering Cancer Center
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#

bin_PROGRAMS = miranda
miranda_SOURCES = miranda.c output.c scan.c swat.c thermo.c utils.c ExpString.c
noinst_HEADERS = miranda.h

LDADD = -L/usr/local/lib -lRNA -lm

INCLUDES = -I/usr/local/include/ViennaRNA

libpymirandadir = $(pkgpythondir)

pyexec_LTLIBRARIES = libpymiranda.la
libpymiranda_la_LDFLAGS = -module -avoid-version -export-symbols-regex PyInit_libpymiranda
libpymiranda_la_CFLAGS = $(PYTHON_CFLAGS)
libpymiranda_la_LIBADD = $(PYTHON_LIBS)
libpymiranda_la_SOURCES = pymiranda.c

EXTRA_DIST = libpymiranda.h
