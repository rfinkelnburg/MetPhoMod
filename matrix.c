/*
  IMPLEMENTATION MODULE Matrix
  Dieses Modul fasst einige wichtige Matrixoperationen zusammen.
*/

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"

Matrix AllocateMatrix(int ni, int nj)
{
  double *a, **ap;
  if (!(a = (double *)malloc(ni * nj * sizeof(double))))  return (NULL);
  if (!(ap = (double **)malloc(ni * sizeof(double *))))  return (NULL);
  while (ni--)  ap[ni] = a + ni * nj;
  return (ap);
}


int FreeMatrix(Matrix a)
{
  free(*a); free(a);
  return (0);
}

int ZeroMatrix(int ni, int nj, Matrix m)
{
  while (ni--)  memset((char *)m[ni], 0, nj * sizeof(double));
  return (0);
}

int UnityMatrix(int ni, Matrix m)
{
  double *p;
  int i, j;
  for (i = ni; i--; )  {
    p = *m++;
    for (j = ni; j--; p++)
      *p = (double)(i == j);
  }
  return (0);
}

int CopyMatrix(int ni, int nj, Matrix a, Matrix b)
{
  nj *= sizeof(double);
  while (ni--)  memcpy((char *)b++, (char *)a++, nj);
  return (0);
}

int TransposeMatrix(int ni, int nj, Matrix a, Matrix b)
{
  int i, j, n;
  double swap;
  n = (ni < nj ? ni : nj);
  for (i = n; i--; )
    for (j = i + 1; j < n; j++)  {
      swap = a[i][j];
      a[i][j] = b[i][j];
      b[i][j] = swap;
    }
  for (i = n; i < ni; )
    for (j = nj; j--; )
      b[j][i] = a[i][j];
  for (j = n; j < nj; )
    for (i = ni; i--; )
      b[j][i] = a[i][j];
  return (0);
}

int AddMatrix(int ni, int nj, Matrix a, Matrix b, Matrix c,
	      double af, double bf)
{
  int i, j;
  for (i = ni; i--; )
    for (j = nj; j--; )
      c[i][j] = af * a[i][j] + bf * b[i][j];
  return (0);
}

int MultMatrix(int ni, int nj, int nk, Matrix a, Matrix b, Matrix c)
{
  int i, j, k;
  double sum;
  for (i = ni; i--; )
    for (j = nj; j--; )  {
      sum = 0.;
      for (k = nk; k--; )
        sum += a[i][k] * b[k][j];
      c[i][j] = sum;
    }
  return (0);
}

#define TINY 1.E-100

int LUdecomp(int n, Matrix a, int *indx, double *d, double *vv)
{
  int i, imax, j, k;
  double big, dum, sum, temp;
  *d = 1.0;
/*  for (i = 0; i < n; i++)  {
    big = 0.0;
    for (j = 0; j < n; j++)
      if ((temp = fabs(a[i][j])) > big)  big = temp;
    if (big == 0.)  return (1);  Die Matrix ist singulaer.
    vv[i] = 1./big;
  }  */
  for (j = 0; j < n; j++)  {
    for (i = 0; i < j; i++)  {
      sum = a[i][j];
      for (k = 0; k < i; k++)  sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++)  {
      sum = a[i][j];
      for (k = 0; k < j; k++)  sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum = fabs(sum)) > big)  {
        big = dum;
        imax = i;
      }
    }
    if (j != imax)  {
      for (k = 0; k < n; k++)  {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -*d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0)  return (1);  /* a[j][j] = TINY; */
    if (j < n-1)  {
      dum = 1.0/(a[j][j]);
      for (i = j+1; i < n; i++)  a[i][j] *= dum;
    }
  }
  return (0);
}

void LUbackSub(int n, Matrix a, int *indx, double *b)
{
  int i, ii, ip, j;
  double sum;
  ii = -1;
  for (i = 0; i < n; i++)  {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii >= 0)  for (j=ii; j <= i-1; j++)  sum -= a[i][j] * b[j];
    else if (sum)  ii = i;
    b[i] = sum;
  }
  for (i = n-1; i >= 0; i--)  {
    sum = b[i];
    for (j = i+1; j < n; j++)  sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}

void InvertMatrix(int n, Matrix a, Matrix inv, int *indx)
{
  int i;
  UnityMatrix(n, inv);
  for (i = n; i--; )
    LUbackSub(n, a, indx, inv[i]);
  TransposeMatrix(n, n, inv, inv);
}
