/*
   MODULE mcparallel.h
   Die meteochem-Schnittstelle zur PVM3-Bibliothek.
*/

typedef enum {NODEINIT=1, TIMEDATA, GRIDDATA, OPTIONDATA, ENVDATA,
	      INITDATA, BORDERDATA, GROUNDDATA, COMMAND, STATUS, AVERAGES,
	      BORDERVALUES, VALUES, TOTHELEFT, TOTHERIGHT, ADVECTIONDATA,
              CHEMISTRYDATA, EMISSIONDATA, DEPOSITIONDATA, NOTIFYDATA,
              RESTARTDATA}
	   CommunicationClass;

typedef enum {DO_TIME_STEP, SENDGRIDVAL, SENDMESHVAL, SENDGRID,
              SENDMESH, SENDGROUND, EXIT}  Command;

#define OVERLAP 2

extern int master, workers, parallel, myparent, leftbrother, rightbrother,
    	   mfirstx, mnx, mlastx;
extern BOOL leftest, rightest;
extern char nodename[80];
extern double *meshcache;

void IsMaster(void);

void InitializeWorkers(char *progname, BOOL debug);

void SendDataToWorkers(void);

void ReadDataFromMaster(void);

void SendCommand(void);

void SendStatus(int i);

int IsPlausible(double *energy);

Command GetCommand(void);

void GetAveragesFromWorkers(double *sum);

void SendAveragesToWorkers(double *bar);

void ChangeAveragesWithMaster(double *bar);

void ChangeTai(int *tai);

void ChangeVectMult(double *sum);

#ifdef MCPARSE
double GetValFromWorker(VarDesc *v, int i, int j, int k, BOOL cache);
#endif

void ReadMeshFromWorkers(double *m, int z0, int dims);

void SendMeshToWorkers(double *m, int z0);

void SendMeshToMaster(double *m, int z0, int dims);

void GetMeshFromMaster(double *m, int z0);

void SendWallToWorker(int net, int x, int tid, CommunicationClass c);

void SendPressureWall(double *x, int xl, int tid, CommunicationClass c);

double *GetGridVarFromWorker(int et);

void GetWallFromWorker(int net, int x, int tid, CommunicationClass c);

void GetPressureWall(double *x, int xl, int tid, CommunicationClass c);

void ReadGroundFromWorker(void);

void GetGroundFromWorker(void);

void GetGroundFromMaster(void);

void SendGroundToMaster(void);

void SendGroundToWorkers(void);
/*
void GetLayerFromWorkers(double *l);

void SendLayerToWorkers(double *l);

void ChangeLayerWithMaster(double *l);
*/

#ifdef MCCDFIN

void PackVariable(double *d, BorderVarType t, InputSection s);

void PackCoordVariable(VarDesc *v);

void UnpackVariable(double *d);

#endif

void SendRestartName(char *name);

BOOL GetRestartName(char *name);

void SendToAll(int id);

void SendSignal(int sig);
