#
# Copyright 2013-2021 INRIA
#

##
## definition of compilation constants
##

## ---------------------------------------------
## compiler to use
## ---------------------------------------------

#CPP              = g++-4.8
#CPP              = clang++

## ---------------------------------------------
## compilation flags
## ---------------------------------------------
#
#CFLAGS      = -g
CFLAGS      = -O3 -pedantic -Wall -Wsign-compare -msse2 -std=c++98
#CFLAGS      = -O0
#EIGEN_FLAGS = -DEIGEN_DONT_VECTORIZE
EIGEN_FLAGS =

## ---------------------------------------------
## includes and libraries
## ---------------------------------------------

LEXLS_INCLUDE = ../../include
EIGEN_INCLUDES = -I/usr/include/eigen3

###EOF
