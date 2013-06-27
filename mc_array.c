/*
   MODULE MC_Array
   Implements dynamic array variables in one two and three Dimensions
   : Implementation
*/

#ifndef INCLUDE_MC_ARRAY
#include "mc_array.h"
#include "mc_array.t"
#endif

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

MemoryHolder::MemoryHolder(void *buff)
{
  base = buff;
  ncopy = new int;
  *ncopy = 1;
}

MemoryHolder::MemoryHolder(MemoryHolder *h)
{
  base = h->base;
  ncopy = h->ncopy;
  (*ncopy)++;
}

MemoryHolder::~MemoryHolder(void)
{
  if (!--(*ncopy))
    delete base, ncopy;
}

template class Array1D<double>;
template class Array2D<double>;
template class Array3D<double>;
