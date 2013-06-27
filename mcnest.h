/*
   MODULE mcnest.h
   Dieses Modul unterstuetzt das Herstellen von Nest-Files. Berphomod in
   seiner aktuellen Version unterstuetzt zwar das Nesten nicht wirklich.
   Mit den Routinen in diesem Modul ist es aber einfach moeglich, Subdomains
   zu definieren, fuer die die entsprechenden Randwerte dann in entsprechende
   Dateien rausgeschrieben werden.
   
   This modules supports writing of Nest-Files. While Berphomod in its actual
   version is not able to do a real time nesting, the routines defined here,
   support the creation of files, which describe the borders of a subdomain,
   and which can easily be used for a new run on this new domain.
   
   mcnest is a part of Berphomod, the comprehensive eulerian atmospheric
   air pollution simulation program.
*/

#define NNESTVAR  8

extern char nestvarname[NNESTVAR][20];
extern double dummynestvar;


char *OnNestDomain(RelyCmmd cmd, VarDesc *v);
/* Tests if a subdomain was properly established. Designed to work with
   Module "mcparse" */

int DefineNestDomain(char *basename, int section);
/* Establishes a Subdomain. */

int EndOfNesting(void);

void CreateDomainFiles(BOOL reopen);
/* Creates the CDF-Files needed for nesting */

void WriteToDomainFiles(long actime);
/* Writes actual values into the nest-files */

void CloseDomainFiles(void);
/* Closes all domain files that were opened previousely */
