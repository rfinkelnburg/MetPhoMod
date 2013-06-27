/*
   Implementation Module mcfloat
*/

#define NOPROTO

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <pwd.h>
#ifdef SVR4
 #include <floatingpoint.h>
#endif
#include <signal.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#ifdef AIX
  #include <fptrap.h>
#endif
#include "mcglobal.h"
#include "mcfloat.h"
#include "mcprint.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

#define _NO_PROTO

void matherr(int sig)
{
  printf("Floatingpoint exception occured!\n");
  if (timeddump)  {
    printf("Making a dumpfile\n");
    MakeAFullDump();
  }
  CloseOutputFiles();
  exit (1);
}

#ifdef SVR5

void InstallFpeHandler()
{
   ieee_handler("set", "common", SIGFPE_ABORT);
}

#else

void InstallFpeHandler()
{
  struct sigaction oldaction, action;
  action.sa_handler = matherr;
  memset(&action.sa_mask, 0, sizeof(action.sa_mask));
  action.sa_flags = 0;
  //  sigaction(SIGFPE, &action, &oldaction); not really useful
}

#endif

BOOL mailtime = FALSE, dumpnow = FALSE;

void SignalUSR1Handler(int)
{
  int i;
  printf("Signal USR1 received - Updating output-files!\n");
  FlushOutputFiles();
  mailtime = 1;
}

void SignalUSR2Handler(int)
{
  int i;
  printf("Signal USR2 received - Creating a dumpfile!\n");
  dumpnow = 1;
}

void SignalTERMHandler(int)
{
  printf("TERM Signal received.\n");
#ifdef PARALLEL
  if (parallel && master)  {
    printf("sending TERM to all workers...\n");
    SendSignal(SIGTERM);
  }
#endif
  if (timeddump)  {
    printf("Creating dumpfile...\n");
    MakeAFullDump();
  }
  printf("Closing output files...\n");
  CloseOutputFiles();
  printf("and exit with error-code!\n");
  exit (1);
}

void SignalFPEHandler(int)
{
  abort();
  while (1);
}


void InstallSignalHandler()
{
  struct sigaction oldaction, action;
  action.sa_handler = SignalUSR1Handler;
  action.sa_flags = 0;
  memset(&action.sa_mask, 0, sizeof(action.sa_mask));
  sigaction(SIGUSR1, &action, &oldaction);
  action.sa_handler = SignalUSR2Handler;
  sigaction(SIGUSR2, &action, &oldaction);
  action.sa_handler = SignalTERMHandler;
  sigaction(SIGTERM, &action, &oldaction);
  sigaction(SIGQUIT, &action, &oldaction);
  sigaction(SIGINT, &action, &oldaction);
  action.sa_handler = SignalFPEHandler;
  sigaction(SIGFPE, &action, &oldaction);
}

void MailTheTime()
{
  int p[2];
  char cmmt[256], *uid;
  static char uname[30];
  pipe(p);
  if (fork())  {
    close(p[0]);
    sprintf(cmmt, "The program meteochem is successfully working at\n%li seconds.\n", actime);
    write(p[1], cmmt, strlen(cmmt)+1);
    close(p[1]);
    wait(NULL);
  }
  else  {
    close(0); dup(p[0]);
    close(p[0]); close(p[1]);
    uid = getpwuid(geteuid())->pw_name;
    printf("Sending mail to %s\n", uid);
    execlp("mail", "mail", uid, (char *)0);
    printf("Start of mail did not succeed!\n");
    exit (0);
  }
}
