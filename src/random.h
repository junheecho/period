/*
the header file defines programs for (psuedo-)random real numbers following uniform distribution in [0,1] and gaussian normal:
*/
#ifndef RANDOM_H
#define RANDOM_H

#include "iRRAM/lib.h"
#include "iRRAM/core.h"
#include <utility>
using namespace iRRAM;


REAL gaussian_real();
REAL uniform_real();

#endif
