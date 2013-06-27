/*
   DEFINITION MODULE mcfloat
   Exception Handling von mcfloat
*/

extern BOOL dumpnow, mailtime;

#ifndef NOPROTO

void InstallFpeHandler(void);
/* Installiert den normalen FPE-Handler */

void InstallSignalHandler(void);
/* Installiert die verschiedenen Signal-Handler */

void MailTheTime(void);

#else

void InstallFpeHandler();
/* Installiert den normalen FPE-Handler */

void InstallSignalHandler();
/* Installiert die verschiedenen Signal-Handler */

void MailTheTime();

#endif
