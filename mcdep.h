/*
   DEFINITION MODULE mcdep
   In diesem Modul wird die Deposition von Schadstoffen gerechnet.
*/

#ifdef PARALLEL

void SendDepositionData(void);

void GetDepositionData(void);

#endif

VarDesc *DepositionVariable(VarDesc *v, VarDesc *vh);

void InitDeposition(void);

void Deposit(long tinc);

