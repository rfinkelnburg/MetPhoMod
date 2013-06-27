/*
  DEFINITION MODULE Matrix
  Dieses Modul fasst einige wichtige Matrixoperationen zusammen.
*/

#ifndef INCLUDE_MATRIX
#define INCLUDE_MATRIX

typedef double **Matrix;

Matrix AllocateMatrix(int ni, int nj);
/* Alloziert den Speicherplatz fuer eine Matrix mit ni Zeilen und nj
   Spalten. Eine Matrix wird als ein ARRAY OF POINTER TO ARRAY OF LONGREAL
   definiert. (Gemaess dem Vorschlag aus den numerical recipes) */

int FreeMatrix(Matrix a);
/* Gibt den Speicherplatz der Matrix a wieder frei. */
   
int ZeroMatrix(int ni, int nj, Matrix m);
/* Loescht den Inhalt der Matrix m. */

int UnityMatrix(int ni, Matrix m);
/* Erzeugt eine quadratische Einheitsmatrix der Dimension ni*ni. */

int CopyMatrix(int ni, int nj, Matrix a, Matrix b);
/* Kopiert die Matrix a nach b */

int TransposeMatrix(int ni, int nj, Matrix a, Matrix b);
/* Erstellt in b die Transponierte Matrix zu a. Die Dimension von
   a ist ni*nj, diejenige von b nj*ni. a und b duerfen identisch sein. */

int AddMatrix(int ni, int nj, Matrix a, Matrix b, Matrix c, double af, double bf);
/* Addiert Matrix af * a und bf * b.  c ist das Resultat diese Matrix darf
   mit a oder b identisch sein. a oder b wird dann allerdings zerstoert. */

int MultMatrix(int ni, int nj, int nk, Matrix a, Matrix b, Matrix c);
/* Multipliziert die Matrix a mit der Matrix b. Das Resultat wird in c
   gespeichert. Das Format von a ist ni*nk, dasjenige von b, nk*nj und das
   von c somit ni*nj. c darf nicht mit a oder b identisch sein! */

int LUdecomp(int n, Matrix a, int *indx, double *d, double *vv);
/* Fuehrt eine LU-Decomposition nach Crout durch. Dieser Algorythmus stammt
   aus den numerical recipes.
   n ist die Dimension der Matrix, a die Matrix, indx gibt die Vertauschung
   der Reihen an. d ist 1. falls eine grade Anzahl Reihen vertauscht wurde,
   sonst -1. vv ist ein Vektor der Groesse n. Er muss von der aufrufenden
   Prozedur zur Verfuegung gestellt werden. */

void LUbackSub(int n, Matrix a, int *indx, double *b);
/* Fuehrt eine Ruecksubstitution aus. Zusammen mit LUdecomp loest diese
   Prozedur das Gleichungssystem Ax=b. */
   
void InvertMatrix(int n, Matrix a, Matrix inv, int *indx);
/* Invertiert die Matrix a. Die invertierte wird in i abgelegt. */

#endif
