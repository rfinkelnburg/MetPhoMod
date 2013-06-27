/*

   This file is a part of

   MetPhoMod v2.0, the comprehensive three dimensional, prognostic
		   mesoscale atmospheric summer smog model.

   Copyright (C) 1996-1999, Silvan Perego


   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, refer to the WEB page

   http://www.gnu.org/copyleft/gpl.html

*/

#include <stddef.h>
#include <stdlib.h>
#include <math.h>

const int nspecialformula = 5;

char *specialformula[] = {
   "K = M * 6.e-34 * pow(T/300., -2.3)",
   "K = 2.3e-13*exp(600./T) + 1.7e-33*M*exp(1000./T)",
   "K = 3.22e-34*exp(2800./T) + 2.38e-54*M*exp(3200./T)",
   "k1=7.2e-15*exp(785./T); k2=4.1e-16*exp(1440./T); k3=1.9e-33*exp(725./T)*M; K = k1+k3/(1.+k3/k2)",
   "K = 1.5e-13*(1.+2.439e-20*M)"};

double SpecialConst(int idx, double T, double M)
{
  double K, k1, k2, k3;
  switch (idx)  {
#line 39 "racm.inp"
    case 0 : K = M * 6.e-34 * pow(T/300., -2.3);
             break;
#line 48 "racm.inp"
    case 1 : K = 2.3e-13*exp(600./T) + 1.7e-33*M*exp(1000./T);
             break;
#line 49 "racm.inp"
    case 2 : K = 3.22e-34*exp(2800./T) + 2.38e-54*M*exp(3200./T);
             break;
#line 62 "racm.inp"
    case 3 : k1=7.2e-15*exp(785./T); k2=4.1e-16*exp(1440./T); k3=1.9e-33*exp(725./T)*M; K = k1+k3/(1.+k3/k2);
             break;
#line 74 "racm.inp"
    case 4 : K = 1.5e-13*(1.+2.439e-20*M);
             break;
  }
  return (K);
}

