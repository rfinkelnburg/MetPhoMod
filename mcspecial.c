#include <stddef.h>
#include <stdlib.h>
#include <math.h>

int nspecialformula = 5;

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
#line 39 "/home/perego/study/plateau/racm.chem"
    case 0 : K = M * 6.e-34 * pow(T/300., -2.3);
             break;
#line 48 "/home/perego/study/plateau/racm.chem"
    case 1 : K = 2.3e-13*exp(600./T) + 1.7e-33*M*exp(1000./T);
             break;
#line 49 "/home/perego/study/plateau/racm.chem"
    case 2 : K = 3.22e-34*exp(2800./T) + 2.38e-54*M*exp(3200./T);
             break;
#line 64 "/home/perego/study/plateau/racm.chem"
    case 3 : k1=7.2e-15*exp(785./T); k2=4.1e-16*exp(1440./T); k3=1.9e-33*exp(725./T)*M; K = k1+k3/(1.+k3/k2);
             break;
#line 76 "/home/perego/study/plateau/racm.chem"
    case 4 : K = 1.5e-13*(1.+2.439e-20*M);
             break;
  }
  return (K);
}

