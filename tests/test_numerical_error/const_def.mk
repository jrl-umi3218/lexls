##
## definition of compilation constants
##

## ---------------------------------------------
## compiler to use
## ---------------------------------------------

#CPP              = g++-4.8
CPP              = clang++

## ---------------------------------------------
## compilation flags
## ---------------------------------------------
#
#CFLAGS      = -g
#CFLAGS      = -O3 -pedantic -Wall -Wsign-compare -msse2 -std=c++98
CFLAGS      = -O0
EIGEN_FLAGS = -DEIGEN_INITIALIZE_MATRICES_BY_ZERO

## ---------------------------------------------
## define object and library directories in the
## home directory (referred to as MAIN_DIR)
## ---------------------------------------------

LEXLS_INCLUDE = ../../include

## ---------------------------------------------
## includes and libraries
## ---------------------------------------------

EIGEN_INCLUDES = -I/usr/local/include/eigen3

###EOF
