/* program fxy.c */
#include  "vector_n.h"
#include  <math.h>
#include  <pthread.h>
#include  <time.h>
#include  <mpi.h>

#define   AmpLong 200.0f
#define   MAX_THREADS  2

int  socketFile;
char  scratchName[25][96];			/* name of scratch disks*/

typedef  struct procStruct
{
	int threadID;                        /* thread number for debugging */
	int il;
	int nc;
	int xc;
	int free;
	int win;
    float *radon;
    float *xij;                       /* save X matrix numDesign**2*2 */
	float *trace2;
} procStruct;

static const int  success = 1;			/* constants for threads to return */
static const int  failure = -1;

pthread_t  workThreads[MAX_THREADS];	/* worker threads */
procStruct  proc[MAX_THREADS];          /* pointers to processing structures */

pthread_mutex_t  printfMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  fprintfMutex = PTHREAD_MUTEX_INITIALIZER;

float *LLsave;
float *traceIn;
float  *trace,*trace1,*trace2;            /* input trace */
int  *ithdr,*ithdr1,*ithdr2;                      /* integer trace header */
float  *rthdr,*rthdr1,*rthdr2;                    /* float trace header */
float  sampleRate;                /* sample rate of data */
float timee;                       /* trace length (ms) */
int  samples;                     /* number of samples in data */
int  headers;                     /* number of entries in trace header */

int  indexCdp;                    /* cdp header */
int  indexTrc_type;                    /* trc_type header */
int  minCDP;                      /* min CDP in data */
int  minIlin;                     /* first inline */
int  minXlin;                     /* first xline */

int  ncrbi;                       /* number of cdps inline */
int  icrb;                        /* number of cdps inline to fxy */
int  nlines;                      /* number of cdps xline to fxy */
int idmap;                        /* icrb*nlines */

int   lmin;                       /* minimum inline no. to fxy */
int   xmin;                       /* minimum xline no. to fxy */
int   tleng;                      /* length of operator time window default = 1000 ms */
int   taper;                      /* time window taper, default = 200 ms */
int  numDesiF;                    /* number of lines in filter design default=25 */
int  numFilt;                     /* number of lines in filter default=5 */
int  fxyF;                        /* fxy decon option */
int  radonF;                      /* high dip filter option */
int  numThreads;				  /* number of threads to run */

int  numDesiI,numDesiX;                  /* number of lines in filter design default=25 */
float hiradon;                   /* higher radon transform accurace */
int  dipI,dipX;                         /* dip ms/trace */
int  method;
int  overlapF,overlapI,overlapX;  /* overlap between filter calculating */
int  overlap1;                     /* overlap between filter calculating */
int OuPut,nbin;

int numWin;                       /* number of time windows */
int sampleWin;                    /* samples in a time window */
int sampleTap;                    /* samples in a time taper */
int ac,half;                      /* ac=numFilter*numFilter-1, half=ac/2 */
int lenft;                        /* t direction no. for fft */
int *lhi,*low,*lowd,*lhid,*lhw;               /* lhw=lhi-low+1 */
int *diplimI,*diplimX;
float *stDipI,*stDipX;
float *resoluI,*resoluX;
int low1,lhi1;                    /* get maximum lhi and minimum low for all windows */ 
int lhwmax,lhwmax2;               /* lhwmax=max */

int freTrace;                     /* numWin*lhw*2 */
float *buff,*buff1;               /* save input traces after FFT of each window
                                     freTrace*icrb*(numDesign+8) */
int   *ihead;                     /* save input headers headers*icrb*(numDesign+8) */
float *rhead;                     

int  cdp,il,xl,oldcdp;                   /* global for cdp, no of inline, no of xline */
int outLine,ouStart;                      /* number of output lines */
int numGroup,numGroupR,numGroupR_s;          /* total group number, current group number */
int xGroup,xGroupR,xcGroup;               /* total group number for xline, current group number */

int *idx;                         /* used in solving equation */
float vv[512];                    /* used in solving equation */
float *rrij;                      /* save correlation matrix numFilter**4 */
float *RR,*G;                     /* matrix and right of equation, dimension=ac */ 
int failSolve;

   int fre;                       /* frequency */
   int ndipI,ndipX;                /* radon dip number */
float dpI,dpX;

int  numFiles,numHeads,numAmpls,numBins;                    /* number of files */
int  tracesPerFile;               /* number of traces in each file */
int  tracesPerHead,tracesPerAmpl,tracesPerBin;               /* number of traces in each file */
char  tracN[200][128];            /* name of tracF */
char  trarN[200][128];            /* name of tracF */
char  headN[200][128];                 /* name of headF */
char  amplN[200][128];                 /* name of headF */
char  binN[100][200][128];                 /* name of headF */

FILE *tracF[200],*trarF[200],*headF[200],*amplF[200],*binF[100][200];
int fxyDone,radonDone;

FILE *fp,*fp1; 
int *mute,sign;
float *absR;
int numDesi,numDesi2,numDesiR2;

/* 02,06,11 */
int  indexpstmBin;                /* location of pstm offset bin number header */
int  indexEnd_ens;                /* end os ensamples flag */
int  indexSeqno;                  
int idtyp;                        /* input data type */
int bin;
int shtotal,shtotal1,shtotal4;
int xlS,ileng,xlSo;
int saveLine, xsave, saveLine_old,xsave_old;
int sx1;
int memo;
int outS[30],outE[30],saveS[30],saveE[30],saveS0[30],saveE0[30];
int *offDisk,*offDisk0;
int nbuff0,nbuff10;
  typedef struct cdpbin            
  {
    int    seq;                    
    int    cdp;                    
    int bin;
  } cdpbin;
  cdpbin   *Cdpbin;                 

int AmpBa; // amplitude balance 
float AmpParm;
int sampleAs; // balancing samples 
int gap,ntolA;
float *ampl,tolAmin,tolAmax;

int FreBa; /* 2004,12 */
float FreParm;
int lenft0,naveFreq,lhi0;
float *aveFreq,*aveFout;

int timeseq[7];
float sDipI[7],sDipX[7];
float rsoluI[7],rsoluX[7];
float qmax[7],qmaxm[7],qminp[7],qmin[7],qmax1,qmin1;  
int numgroup;

/* after MPI */
int numprocs,myid;
MPI_Status status;
float *fconv;
char *cconv,*cint;
int *iint;

void sentParm();
void  GetSocketData( char *theData, int length );
void  ProgramErr( char *errString );
void  PrintTime( void );
void  GetGlobals( void );
void  lubksb( float *a, int n, int *indx, float *b );
int  ludcmp( float *a, int n, int *indx, float *vv );
void  PreFxy();
void binStart();
void getAmpl();
void fixAmpl(int xl);
void FreAmp();
void  Fxy();
void  traceDisk(int cdp,int kk);
   void radov(procStruct *ps);
   void iradon(procStruct *ps);
   void iradov(procStruct *ps);
   void radov3d(procStruct *ps);
void  procThread2();		/* routine threads run */
   void iradon3d(procStruct *ps);
   void iradov3d(procStruct *ps);
void Xij3d(procStruct *ps);
void Radon3d(int i,int j);
void FxyPerform(int numG,int xG,int i,int j);
void Xij(int j,int xcG);
void Rij();
void Get_rij(int i1,int i2,int j1, int j2,int numrj2,int numxi2,int numxj2);
void Rsolve(int jj);
void strongD(procStruct  *ps,int np,int ni);
void  procThread();		/* routine threads run */
void XijLine(procStruct *ps);
void reXijLine(procStruct *ps);
void lineRadon(int i,int j);
void ouTrace(int numG,int xG,int i,int j);
void Freq(int ii,int xx,int ww,int ff,float *a);
void caddr(float *trace,float *trace1,float *a,int j);
void caddc(float *trace,float *trace1,float *a,int j);
void caddi(float *trace,float *trace1,float *a,int j);
void caddir(float *trace,float *trace1,float *a,int j);
void vclrf(float *trace,int i);
void vclri(int *trace,int i);
void vmovf(float *trace,float *trace1,int i);
void vmovi(int *trace,int *trace1,int i);
void vabsf(float *trace,int j,int k);
void  get_trace( float *trace, int seq,int num,int kk);
void  put_trace( float *trace, int seq,int num,int kk);
void  get_ampl( float *trace, int seq, int num );
void  put_ampl( float *trace, int seq, int num );
void  get_head( float *trace, int seq,int num );
void  put_head( float *trace, int seq,int num );
void  get_bin( float *trace, int seq,int num,int bin );
void  put_bin( float *trace, int seq,int num,int bin );
void  OpenTracFile();
void  OpenHeadFile();
void  OpenAmplFile();
void  OpenBinFile();
void stoepf (int n, float *r, float *g, float *f, float *a);
void get_LLH1(procStruct *ps);
void get_LLH(procStruct *ps);
void  procThreadLL();		/* routine threads run */
void get_LLsave();
void Get_aveFreq();
void Get_aveFout();
void Get_tol();
void Perform();
void ouPuts(int numG,int i);
void get_buff(int numG,int xG,int numcG,int xcG,int kk);
void put_buff(int numG,int numcG,int kk);
void freparm(int *fre);
void Abalanc(float *trace,int xl);
void get_out_one();
void getMaster(float *trace,int nmaxfre,int lh,int lh01);
void get_out();
void get_freq();
void radonParm();
void diskOpen();
void  masterOut();
void Master0(float *buff,cdpbin *Cdpbin,int *kkb,int *kk,int il);
void  slaveOut();
int compare( const void *a1, const void *a2 );
void  masterOne();
void  slaveOne();

int  main( int argc, char *argv[])
{
	int i,j,result, sockPort;
	char temp[9];
    struct sockaddr_in  socketAddr;
    struct hostent  *host;
	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Get_processor_name(processor_name,&namelen);
	if(myid==0) printf("number of process used %d \n",numprocs);
	fprintf(stderr,"Process %d on %s \n",myid,processor_name);

	if(myid==0)
	{
	host = gethostbyname( argv[1] );
	sockPort = atoi( argv[2] );
/*	printf( "sockPort %d argv[0] %s\n",sockPort,argv[0] );
	printf( "sockPort %d argv[1] %s\n",sockPort,argv[1] );
	printf( "sockPort %d argv[2] %s\n",sockPort,argv[2] );*/
	strcpy( scratchName[0], argv[2] );

	i = getpid();
	sprintf( temp, "%d", i );
	strcat( scratchName[0], temp );

	if( host == ( struct hostent * )NULL )
		{
		perror( "get host by name error" );
		exit( 2 );
		}

	memset( &socketAddr, 0, sizeof ( socketAddr ) );
	socketAddr.sin_family = AF_INET;
	memcpy( &socketAddr.sin_addr, host->h_addr, host->h_length );
	socketAddr.sin_port = htons( sockPort );

	if( ( socketFile = socket( AF_INET, SOCK_STREAM, 0 ) ) < 0 )
		{
		perror( "generate error" );
		exit( 3 );
		}
	printf( "socket created\n" );

	if( connect( socketFile, ( struct sockaddr * )&socketAddr, sizeof( socketAddr ) ) < 0 )
		{
		perror( "connect error" );
		exit( 4 );
		}
	printf( "socket connected\n" );

/* AIXTHREAD_SCOPE=S;*/
   GetGlobals();
	}
	sentParm();
    MPI_Bcast(&scratchName[0],100,MPI_CHAR,0,MPI_COMM_WORLD);

     PreFxy();
     Fxy();

   if(myid>0)
   {
   for( i = 0; i < numFiles; i++ )
   {
     fclose( tracF[i] );
     remove( tracN[i] );
   }
   if((fxyF && radonF) || radonF==3)
   for( i = 0; i < numFiles; i++ )
   {
     fclose( trarF[i] );
     remove( trarN[i] );
   }

   for( i = 0; i < numHeads; i++ )
   {
     fclose( headF[i] );
     remove( headN[i] );
   }
   if(AmpBa)
   for( i = 0; i < numAmpls; i++ )
   {
     fclose( amplF[i] );
     remove( amplN[i] );
   }
   for( j = 0; j < nbin; j++ )
   {
	 for(i=0;i<numBins;i++)
	 {
     fclose( binF[j][i] );
     remove( binN[j][i] );
	 }
   }

   }

	if(myid==0)
	{
/*	MPI_Recv(&i,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);*/
    result = 1;
    write( socketFile, ( char* )&result, 4 );

/*	sleep( 10 );*/
	}

	MPI_Finalize();
    close( socketFile );

    if(myid==0) 
	{
      PrintTime();
      printf("client quitting \n");
	}
    return(0);
}             

void  GetGlobals( void )
{
  int i;
//   if((fp=fopen("argp","w"))==NULL)
//   printf("cannot open file \n");
//   if((fp1=fopen("argpp","w"))==NULL)
//   printf("cannot open file \n");
   fconv = FloatVector( 0, 1);
   cconv = (char *)fconv;
   iint = IntVector( 0, 1);
   cint = (char *)iint;
   GetSocketData( ( char* )&sampleRate, 4 );
   GetSocketData( ( char* )&samples, 4 );
   GetSocketData( ( char* )&headers, 4 );
   GetSocketData( ( char* )&shtotal, 4 );
   GetSocketData( ( char* )&minCDP, 4 );
   GetSocketData( ( char* )&minIlin, 4 );
   GetSocketData( ( char* )&minXlin, 4 );
   GetSocketData( ( char* )&lmin, 4 );
   GetSocketData( ( char* )&xmin, 4 );

   GetSocketData( ( char* )&idtyp, 4 );
   GetSocketData( ( char* )&icrb, 4 );
   GetSocketData( ( char* )&nlines, 4 );
   GetSocketData( ( char* )&idmap, 4 );
   GetSocketData( ( char* )&ncrbi, 4 );

   GetSocketData( ( char* )&tleng, 4 );
   GetSocketData( ( char* )&taper, 4 );
   GetSocketData( ( char* )&fxyF, 4 );
   GetSocketData( ( char* )&numDesiF, 4 );
   GetSocketData( ( char* )&numFilt, 4 );
   GetSocketData( ( char* )&overlapF, 4 );

   overlap1=3;
   GetSocketData( ( char* )&radonF, 4 );
   numThreads=1;
   GetSocketData( ( char* )&hiradon, 4 );
   GetSocketData( ( char* )&method, 4 );
   GetSocketData( ( char* )&numDesiI, 4 );
   GetSocketData( ( char* )&dipI, 4 );
   GetSocketData( ( char* )&overlapI, 4 );
   GetSocketData( ( char* )&numDesiX, 4 );
   GetSocketData( ( char* )&dipX, 4 );
   GetSocketData( ( char* )&overlapX, 4 );
   GetSocketData( ( char* )&numgroup, 4 );
   GetSocketData( ( char* )&timeseq[0], 4*numgroup );
   GetSocketData( ( char* )&qmin[0], 4*numgroup );
   GetSocketData( ( char* )&qminp[0], 4*numgroup );
   GetSocketData( ( char* )&qmax[0], 4*numgroup );
   GetSocketData( ( char* )&qmaxm[0], 4*numgroup );
   GetSocketData( ( char* )&sDipI[0], 4*numgroup );
   GetSocketData( ( char* )&rsoluI[0], 4*numgroup );
   GetSocketData( ( char* )&sDipX[0], 4*numgroup );
   GetSocketData( ( char* )&rsoluX[0], 4*numgroup );
/*   numDesiX=numDesiI;
   stDipX=stDipI;
   resoluX=resoluI;
   dipX=dipI;
   overlapX=overlapI;*/
   if(radonF==4) overlapI=3;
   GetSocketData( ( char* )&AmpBa, 4 );
   GetSocketData( ( char* )&AmpParm, 4 );
   GetSocketData( ( char* )&FreBa, 4 );
   GetSocketData( ( char* )&FreParm, 4 );
   GetSocketData( ( char* )&OuPut, 4 );
   GetSocketData( ( char* )&nbin, 4 );

   GetSocketData( ( char* )&indexCdp, 4 );
   GetSocketData( ( char* )&indexTrc_type, 4 );
   if(!idtyp)
   {
     GetSocketData( ( char* )&indexSeqno, 4 );
     GetSocketData( ( char* )&indexEnd_ens, 4 );
     GetSocketData( ( char* )&indexpstmBin, 4 );
   }

   printf( "\n Linux Socket VHF 6.0 job parameters\n" );
   printf( "--------------\n\n" );
   printf( "minimum inline no. to fxy  %7d\n", lmin );
   printf( "minimum xline no. to fxy  %7d\n", xmin );
   printf( "length of operator time window %7d\n", tleng );
   printf( "time window taper %7d\n", taper );
   printf( "total put values = %7d in each sequence \n",numgroup);
   printf( "time sequence \n");
   for(i=0;i<numgroup;i++) printf( " %7d ", timeseq[i] );
   printf( " \n");
   printf( "qmin minimum frequency sequence \n");
   for(i=0;i<numgroup;i++) printf( " %6.1f ", qmin[i] );
   printf( " \n");
   printf( "larger than minimum frequency sequence for taper \n");
   for(i=0;i<numgroup;i++) printf( " %6.1f ", qminp[i] );
   printf( " \n");
   printf( "maximum frequency sequence \n");
   for(i=0;i<numgroup;i++) printf( " %6.1f ", qmax[i] );
   printf( " \n");
   printf( "less than maximum frequency sequence for taper \n");
   for(i=0;i<numgroup;i++) printf( " %6.1f ", qmaxm[i] );
   printf( " \n");
   if(fxyF==0) printf( "no fxy decon filter used\n" );
   else
   {
     printf( "number of lines for fxy filter design %7d\n", numDesiF );
     printf( "number of lines for fxy decon filter used %7d\n", numFilt );
     printf( "overlap between FXY group = %5d \n",overlapF);
   }
   if(radonF==0) printf( "no radon filter used \n" );
   if(radonF>0) 
   {
     printf( "Multiplier for Radon transform accuracy %7.3f\n", hiradon );     
	 if(method==0) printf( "you use more noisy Radon filter method \n");
	 else printf( "you use better noise suppression Radon filter method \n");
	 if(radonF!=2)
	 {
       printf( "the following parms using inline direction \n");
       printf( "number of lines for radon filter design %7d\n", numDesiI );
	   printf( "dip(ms/trace) for line radon filter used %7d\n", dipI );
       printf( "overlap between Radon group in inline direction = %5d \n",overlapI);
       printf( "increasing resolutoin coefficient sequence for inline direction \n");
       for(i=0;i<numgroup;i++) printf( " %6.2f ", rsoluI[i] );
       printf( " \n");
       printf( "enhancing stronger dips sequence for inline direction \n");
       for(i=0;i<numgroup;i++) printf( " %6.2f ", sDipI[i] );
       printf( " \n");
	 }
	 if(radonF==2 || radonF==3)
	 {
       printf( "the following parms using xline direction \n");
       printf( "number of lines for radon filter design %7d\n", numDesiX );
	   printf( "dip(ms/trace) for line radon filter used %7d\n", dipX );
       printf( "overlap between Radon group = %5d \n",overlapX);
       printf( "increasing resolutoin coefficient sequence for xline direction \n");
       for(i=0;i<numgroup;i++) printf( " %6.2f ", rsoluX[i] );
       printf( " \n");
       printf( "enhancing stronger dips sequence for xline direction \n");
       for(i=0;i<numgroup;i++) printf( " %6.2f ", sDipX[i] );
       printf( " \n");
	 }
   }
   if(AmpBa) printf( "amplitude balacing performing uisng parm %8.2f \n",AmpParm);
   else printf( "amplitude balacing not performing \n");
   if(FreBa==1) printf( "low frequency balacing performing, the parm = %8.2f \n",FreParm);
   else printf( "low frequency balacing not performing \n");
   if(!idtyp)
   {
   if(OuPut)
	   printf( "output: first sort is pstmBin, second sort is CDP \n");
   else
   {
	   printf( "output: first sort is CDP, second sort is pstmBin  \n");
	   printf( "bin number = %7d \n",nbin);
   }
   }

   printf( "\n Linux Socket VHF 6.0 geometry parameters\n" );
   printf( "--------------\n\n" );
   printf( "first CDP               %7ld\n", minCDP );
   printf( "first inline            %7ld\n", minIlin );
   printf( "first xline             %7ld\n", minXlin );
   printf( "sampleRate              %8.2f\n",sampleRate );
   printf( "samples                 %7d\n",  samples );
   printf( "headers                 %7d\n",  headers );

   printf("you process total inline number %7d \n",nlines);
   printf("you process total crossline number %7d \n",icrb);
   printf("you process total CRB number %7d \n",idmap); 
   if(overlapI>overlapX) overlapI=overlapX;
   if(overlapX>overlapI) overlapX=overlapI;

}

void sentParm()
{

/*   if(myid==8) 
   {
   if((fp=fopen("argp8","w"))==NULL)
   printf("cannot open file \n");
   }
   if(myid==7) 
   {
   if((fp=fopen("argp7","w"))==NULL)
   printf("cannot open file \n");
   }
   if(myid==6) 
   {
   if((fp=fopen("argp6","w"))==NULL)
   printf("cannot open file \n");
   }
   if(myid==5) 
   {
   if((fp=fopen("argp5","w"))==NULL)
   printf("cannot open file \n");
   }*/
/*   if(myid==4) 
   {
   if((fp=fopen("argp4","w"))==NULL)
   printf("cannot open file \n");
   }
   if(myid==3) 
   {
   if((fp=fopen("argp3","w"))==NULL)
   printf("cannot open file \n");
   }
   if(myid==2) 
   {
   if((fp=fopen("argp2","w"))==NULL)
   printf("cannot open file \n");
   }
   if(myid==1) 
   {
   if((fp=fopen("argp1","w"))==NULL)
   printf("cannot open file \n");
   }*/
   if(myid>0)
   {
   fconv = FloatVector( 0, 1);
   cconv = (char *)fconv;
   iint = IntVector( 0, 1);
   cint = (char *)iint;
   }
  MPI_Bcast(&sampleRate,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&samples,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&headers,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&shtotal,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&minCDP,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&minIlin,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&minXlin,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&lmin,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&xmin,1,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Bcast(&idtyp,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&icrb,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nlines,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&idmap,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&ncrbi,1,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Bcast(&tleng,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&taper,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&fxyF,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&numDesiF,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&numFilt,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&overlapF,1,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Bcast(&radonF,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&hiradon,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&method,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&numDesiI,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&dipI,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&overlapI,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&numDesiX,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&dipX,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&overlapX,1,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Bcast(&numgroup,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&timeseq[0],numgroup,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&qmin[0],numgroup,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&qminp[0],numgroup,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&qmax[0],numgroup,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&qmaxm[0],numgroup,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&sDipI[0],numgroup,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&rsoluI[0],numgroup,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&sDipX[0],numgroup,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&rsoluX[0],numgroup,MPI_FLOAT,0,MPI_COMM_WORLD);

  MPI_Bcast(&AmpBa,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&AmpParm,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&FreBa,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&FreParm,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&OuPut,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&nbin,1,MPI_INT,0,MPI_COMM_WORLD);

  shtotal1=shtotal*10;
  shtotal4=shtotal*40;
  timee=sampleRate*samples;
  overlap1=3;
  if(numDesiI-overlapI*2<3) overlap1=numDesiI-overlapI*2; 
  if(numDesiX-overlapX*2<3) overlap1=numDesiX-overlapX*2; 
  if(overlap1==1 && fxyF && numFilt>3) numFilt=3;
//  numThreads=2;
  numThreads=1;
  if(!radonF) numThreads=1;
  MPI_Bcast(&indexCdp,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&indexTrc_type,1,MPI_INT,0,MPI_COMM_WORLD);
  if(!idtyp)
  {
    MPI_Bcast(&indexSeqno,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&indexEnd_ens,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&indexpstmBin,1,MPI_INT,0,MPI_COMM_WORLD);
  }
  else
  {
    indexSeqno=0;
    indexEnd_ens=0;
    indexpstmBin=0;
  }
/*   printf("finish sending file \n");*/
}

void get_freq()
{
  int i,j,k;
  int ts,te,ti;
  float a1,a2,a3;

/* frequency calculating for each window */
  low = IntVector(0, numWin );
  lowd = IntVector(0, numWin);
  lhi = IntVector(0, numWin );
  lhid = IntVector(0, numWin );
  lhw = IntVector(0, numWin );
  stDipI = FloatVector(0, numWin );
  stDipX = FloatVector(0, numWin );
  resoluI = FloatVector(0, numWin );
  resoluX = FloatVector(0, numWin );
  lenft=RealPrimeFactorN( sampleWin+100, sampleWin+300 );
  a1=lenft/1000.*sampleRate;
//   printf( "maximum frequency sequence \n");
//   for(i=0;i<numgroup;i++) printf( " %6.1f ", qmax[i] );
//   printf( " \n");
  qmax1=qmax[0];
  for(i=1;i<numgroup;i++) 
  {
	  if(qmax1<qmax[i]) qmax1=qmax[i];
//      fprintf(fp, "i,qmax1 %5d %6.1f \n", i,qmax1 );
  }
  qmin1=qmin[0];
  for(i=1;i<numgroup;i++) if(qmin1>qmin[i]) qmin1=qmin[i];
  low1=(int)(qmin1*a1);
  if(low1<2) low1=2;
  lhi1=(int)(qmax1*a1)+1;
//  fprintf(fp,"lenft,low1,lhi1,qmax1,qmin1 %7d %7d %7d %7.1f %7.1f \n",lenft,low1,lhi1,qmax1,qmin1);

  lhwmax=lhi1-low1+1;
  lhwmax2=lhwmax*2;
  freTrace=numWin*lhwmax2;
  for(i=0;i<numWin;i++)
  {
	ts=i*(sampleWin-sampleTap);
	te=ts+sampleWin;
	ti=(ts+te)/2;
	k=(int)(ti*sampleRate);
	for(j=0;j<numgroup;j++) if(k<=timeseq[j]) break;
	if(j==numgroup) 
	{
	  low[i]=qmin[numgroup-1]*a1;
	  if(low[i]<low1) low[i]=low1;
	  lhi[i]=qmax[numgroup-1]*a1;
	  if(lhi[i]>lhi1) lhi[i]=lhi1;
	  lowd[i]=qminp[numgroup-1]*a1;
	  if(lowd[i]<low[i]) lowd[i]=low[i];
	  lhid[i]=qmaxm[numgroup-1]*a1;
	  if(lhid[i]>lhi[i]) lhid[i]=lhi[i]; 
	  stDipI[i]=sDipI[numgroup-1];
	  stDipX[i]=sDipX[numgroup-1];
	  resoluI[i]=rsoluI[numgroup-1];
	  resoluX[i]=rsoluX[numgroup-1];
	}
	if(j==0) 
	{
	  low[i]=qmin[0]*a1;
	  if(low[i]<low1) low[i]=low1;
	  lhi[i]=qmax[0]*a1;
	  if(lhi[i]>lhi1) lhi[i]=lhi1;
	  lowd[i]=qminp[0]*a1;
	  if(lowd[i]<low[i]) lowd[i]=low[i];
	  lhid[i]=qmaxm[0]*a1;
	  if(lhid[i]>lhi[i]) lhid[i]=lhi[i]; 
	  stDipI[i]=sDipI[0];
	  stDipX[i]=sDipX[0];
	  resoluI[i]=rsoluI[0];
	  resoluX[i]=rsoluX[0];
	}
	if(j>0 && j<numgroup)
	{
	  a3=(float)(timeseq[j]-timeseq[j-1]);
	  if(a3==0.) a2=0.;
	  else a2=(float)(k-timeseq[j-1])/(float)(timeseq[j]-timeseq[j-1]);
	  low[i]=(qmin[j-1]+(qmin[j]-qmin[j-1])*a2)*a1;
	  if(low[i]<low1) low[i]=low1;
	  lhi[i]=(qmax[j-1]+(qmax[j]-qmax[j-1])*a2)*a1;
	  if(lhi[i]>lhi1) lhi[i]=lhi1;
	  lowd[i]=(qminp[j-1]+(qminp[j]-qminp[j-1])*a2)*a1;
	  if(lowd[i]<low[i]) lowd[i]=low[i];
	  lhid[i]=(qmaxm[j-1]+(qmaxm[j]-qmaxm[j-1])*a2)*a1;
	  if(lhid[i]>lhi[i]) lhid[i]=lhi[i]; 
	  stDipI[i]=sDipI[j-1]+(sDipI[j]-sDipI[j-1])*a2;
	  stDipX[i]=sDipX[j-1]+(sDipX[j]-sDipX[j-1])*a2;
	  resoluI[i]=rsoluI[j-1]+(rsoluI[j]-rsoluI[j-1])*a2;
	  resoluX[i]=rsoluX[j-1]+(rsoluX[j]-rsoluX[j-1])*a2;
	}
	lhw[i]=lhi[i]-low[i]+1;
//  fprintf(fp,"i,ts,te,ti,j,k,low[i],lhi[i],lhw[i] %7d %7d %7d %7d %7d %7d %7d %7d %7d \n",i,ts,te,ti,j,k,low[i],lhi[i],lhw[i]);
//  fprintf(fp,"i,stDipI[i],stDipX[i],resoluI[i],resoluX[i] %7d %7.3f %7.3f %7.3f %7.3f \n",i,stDipI[i],stDipX[i],resoluI[i],resoluX[i]);
  }
  if(FreBa>0)
  {
    lenft0=RealPrimeFactorN( samples+300, samples+1300 );
    a1=lenft0/1000.*sampleRate;
    lhi0=(int)(qmax1*a1)+1; 
    aveFreq = FloatVector(0, lhi0 );
    aveFout = FloatVector(0, lhi0 );
//    fprintf(fp,"lenft0,lhi0 %7d %7d \n",lenft0,lhi0);
  }
//  fprintf(fp,"finish freq lhwmax,low1,lhi1,freTrace %7d %7d %7d %7d \n",lhwmax,low1,lhi1,freTrace);
//    fflush(fp);
}

void radonParm()
{
  int i,ndips;
  float a1,a2,dipm;

  diplimI = IntVector(0, lhwmax );
  diplimX = IntVector(0, lhwmax );
//  fprintf(fp,"lhwmax,low1,lhi1 %7d %7d %7d \n",lhwmax,low1,lhi1);

  dipm=500/(qmax1+5.);
  a1=numDesiI*1.4142136;
  dpI=dipm*2.0/a1;
  ndips=a1*0.5+2;
  if(ndips<numDesiI) ndips=numDesiI;

  a2=log(dipI)-log(dipm);
/*    fprintf(fp,"a2,dip %10.4f %5d\n",a2,dip);
    fflush(fp);*/
  if(a2<0.) a2=0.;
  ndipI=0.5*a1*(1.+a2*1.05)+2;
  if((float)ndipI>ndips*1.3) ndipI=(int)(ndips*1.3);
  ndips=(ndips+ndipI)/2;
  ndips=(int)((float)ndips*hiradon);
  ndipI=(int)((float)ndipI*hiradon);
  dpI=dipI/(ndips*sampleRate);

  a1=numDesiX*1.4142136;
  dpX=dipm*2.0/a1;
  ndips=a1*0.5+2;
  if(ndips<numDesiX) ndips=numDesiX;
  a2=log(dipX)-log(dipm);
/*    fprintf(fp,"a2,dip %10.4f %5d\n",a2,dip);
    fflush(fp);*/
  if(a2<0.) a2=0.;
  ndipX=0.5*a1*(1.+a2*1.05)+2;
  if((float)ndipX>ndips*1.3) ndipX=ndips*1.3;
  ndips=(ndips+ndipX)/2;
  ndips=(int)((float)ndips*hiradon);
  ndipX=(int)((float)ndipX*hiradon);
  dpX=dipX/(ndips*sampleRate);

/*    fprintf(fp,"dpI,dpX,ndipI,ndipX %10.6f %10.6f %7d %7d \n",dpI,dpX,ndipI,ndipX);
    fflush(fp);
    fprintf(fp,"dipm,dpX,ndipX,ndips %10.2f %10.4f %7d %7d   \n",dipm,dpX,ndipX,ndips);
    fflush(fp);*/
//    dp=dip/(ndips*sampleRate);
/*    fprintf(fp,"dipI,dipX,numDesiI,numDesiX,overlapI,overlapX %5d %5d %5d %5d %5d %5d \n",dipI,dipX,numDesiI,numDesiX,overlapI,overlapX);
    fprintf(fp,"resoluI,resoluX,stDipI,stDipX %10.2f %10.2f %10.2f %10.2f \n",resoluI,resoluX,stDipI,stDipX);
    fflush(fp);*/

  a1=0.5*lenft/dpI;
  diplimI[0]=ndipI;
  for(i=1;i<lhwmax;i++)
  {
	diplimI[i]=(int)(a1/(i+low1))+1;     /* each freq ndip should < diplim */
	if(diplimI[i]>ndipI) diplimI[i]=ndipI;
//      fprintf(fp,"i,diplimI %5d %5d \n",i,diplimI[i]);
  }
  a1=0.5*lenft/dpX;
  diplimX[0]=ndipX;
  for(i=1;i<lhwmax;i++)
  {
	diplimX[i]=(int)(a1/(i+low1))+1;     /* each freq ndip should < diplim */
	if(diplimX[i]>ndipX) diplimX[i]=ndipX;
//      fprintf(fp,"i,diplimX %5d %5d \n",i,diplimX[i]);
  }
}

void diskOpen()
{
  int i,j,k,m;
  int overI,overX,num;

  i=freTrace+2;
  if(i<samples+2) i=samples+2;

  num=numDesiF;
  if(num<numDesiX) num=numDesiX;
  j=150000000/(i*num);
//  j=1500000/(i*num);
  k=j/(icrb-5);
  if(k>1) ileng=icrb+2*overlap1;
  else 
  {
	overI=overlapF;
	if(overI<overlapI) overI=overlapI;
	if(icrb/j<2) ileng=icrb/2+overlap1+overI;
	else ileng=icrb/4+overlap1+overI;
  }

  if(radonF==4) 
  {
	overI=overX=3;
	if(fxyF>0 && overlapF>3) overI=overlapF;
  }
  else
  {
    if(fxyF>0) overI=overX=overlapF;
	else overI=overX=3;
	if((radonF==1 || radonF==3) && overI<overlapI) overI=overlapI;
	if(radonF>0 && overX>overlapX) overX=overlapX;
  }
//    fprintf(fp,"i,ileng,num,overX,icrb %5d %5d %5d %5d %5d\n",i,ileng,num,overX,icrb);
//    fflush(fp);
  nbuff0=i*num*ileng;
  buff = FloatVector(0,nbuff0);
  vclrf(buff,i*num*ileng );
  nbuff10=i*(num-overX-overlap1)*icrb;
  buff1 = FloatVector(0,nbuff10);
  vclrf(buff1,i*(num-overX-overlap1)*icrb);
  ihead = IntVector(0, (headers+2)*(num-overX-overlap1)*icrb);
  rhead = ( float * )ihead;
  vclrf(rhead,(headers+2)*(num-overX-overlap1)*icrb);

  if(fxyF>0)
  {
/* array for fxy used */
    half=numFilt*numFilt;
    ac=2*half;
    G = FloatVector(0,(half-1)*lhwmax*numWin); /*8000*/
    rrij = FloatVector(0,ac*half/2);
    RR=FloatVector(0,half*half);
  }

    if(radonF==0)  idx=IntVector(0,ac);
    else idx=IntVector(0, numDesi2);
  if(radonF>0 && method==1)
  {
    if(radonF<4) 
    {
   	  LLsave = FloatVector(0,lhwmax*numDesiR2);
	  vclrf(LLsave,lhwmax*numDesiR2);
    }
    else 
    {
	  LLsave = FloatVector(0,lhwmax*numDesiR2*numDesiR2);
	  vclrf(LLsave,lhwmax*numDesiR2*numDesiR2);
    }
  }
  if(radonF==0) proc[0].xij=FloatVector(0,numDesiF*numDesiF*2);
  else
  {
	i=ndipI; if(ndipI<ndipX) i=ndipX;
    for( j = 0; j < numThreads; j++ )
//    for( j = 0; j < 2; j++ )
    {
      proc[j].threadID=j;
	  proc[j].trace2=FloatVector(0,lenft*2);
      if(radonF<4)
      {
	    proc[j].xij=FloatVector(0,numDesi2*2);
	    proc[j].radon=FloatVector(0,i*2*lhwmax2);
      }
      else
      {
	    proc[j].xij=FloatVector(0,numDesiR2*numDesiR2);
	    proc[j].radon=FloatVector(0,ndipI*ndipI*4*lhwmax2);
      }
    }
  }

  if(AmpBa)
  {
	if(sampleRate>3.) { sampleAs=samples/2+2; gap=2;}
	else { sampleAs=samples/4+2; gap=4;}
    tracesPerAmpl = 400000000 / (sampleAs*4);
    numAmpls = 0;
    OpenAmplFile();
    i=((outE[myid]-outS[myid])*icrb-1)/tracesPerAmpl+1;
    if(i>1) for(j=1;j<i;j++) OpenAmplFile();
    if(myid==1)
    printf( "ampl per ampl file %d\n", tracesPerAmpl );
    ampl = FloatVector(0,icrb*sampleAs );
    vclrf(ampl,icrb*sampleAs );
//    fprintf(fp,"sampleAs,gap %7d %5d \n",sampleAs,gap);
//    fflush(fp);
  }

  tracesPerHead = 100000000 / (headers);
  numHeads = 0;
  OpenHeadFile();
  i=((outE[myid]-outS[myid])*icrb-1)/tracesPerHead+1;
  if(i>1)
  {
	for(j=1;j<i;j++) OpenHeadFile();
  }
  if(myid==1)
  printf( "head number of per header file %d\n", tracesPerHead );
  tracesPerFile = 100000000 / (freTrace);
  numFiles = 0;
  OpenTracFile();
  k=((saveE[myid]-saveS[myid])*icrb-1)/tracesPerFile+1;
  if(k>1)
  {
	for(j=1;j<k;j++) OpenTracFile();
  }
  if(myid==1)
  printf( "trace number of per trace file %d\n", tracesPerFile );
  tracesPerBin = 100000000 / (shtotal);
  numBins = 0;
  OpenBinFile();
  m=((outE[myid]-outS[myid])*icrb-1)/tracesPerBin+1;
  if(m>1)
  {
	for(j=1;j<m;j++) OpenBinFile();
  }
  if(myid==1)
  {
  printf( "trace number of per bin file %d\n", tracesPerBin );
  printf( "%5d header files, %5d trace files, bin number =%5d, %5d files opened for each bin\n",i,k,nbin,m );
  }

}

void  PreFxy( void )
{
  int i,j,k;

  sampleWin = (int)((float)tleng/sampleRate);
  sampleTap = (int)((float)taper/sampleRate);
  numWin=samples/(sampleWin-sampleTap);
  if(samples%(sampleWin-sampleTap)>sampleTap) numWin++;
  while(samples-(numWin-1)*(sampleWin-sampleTap)<sampleWin/2)
  {
	tleng+=5;
    sampleWin = (int)((float)tleng/sampleRate);
    numWin=samples/(sampleWin-sampleTap);
    if(samples%(sampleWin-sampleTap)>sampleTap) numWin++;
  }
  numDesi=numDesiF;
  if(numDesiF<numDesiI) numDesi=numDesiI;
  if(numDesi<numDesiX) numDesi=numDesiX;
  numDesi2=numDesi*numDesi;
  numDesiR2=numDesiI*numDesiI;

/*  fprintf(fp,"overlapF,numDesiF,%6d %6d \n",overlapF,numDesiF);
  fprintf(fp,"overlapI,numDesiI,dipI,resoluI,stDipI %7d %7d %7d %7.2f %7.2f \n",overlapI,numDesiI,dipI,resoluI,stDipI);
  fprintf(fp,"overlapX,numDesiX,dipX,resoluX,stDipX %7d %7d %7d %7.2f %7.2f \n",overlapX,numDesiX,dipX,resoluX,stDipX);
  fprintf(fp,"sampleWin,sampleTap,numWin %7d %7d %7d   \n",sampleWin,sampleTap,numWin);
  fprintf(fp,"numDesi,numDesiF,numDesiR,numDesi2,numDesiR2 %7d %7d %7d %7d %7d  \n",numDesi,numDesiF,numDesiI,numDesi2,numDesiR2  );
  fflush(fp);*/

  vclri(outS,30);
  vclri(outE,30);
  vclri(saveS,30);
  vclri(saveE,30);
  {
	i=numprocs-1;
	j=nlines/i;
	outE[1]=j;
	k=nlines%i;
//	fprintf(fp,"i,j,k %5d %5d %5d \n",i,j,k); 
//  fflush(fp);
	if(k>0) outE[1]+=k/2+k%2;
	for(k=2;k<i;k++) outE[k]=outE[k-1]+j;
	outE[i]=nlines;
	for(k=2;k<=i;k++) outS[k]=outE[k-1];
  }
  if(fxyF) 
  {
	for(k=0;k<numprocs;k++)
	{
	  saveS[k]=outS[k]-overlapF;
	  if(saveS[k]<0) saveS[k]=0;
	  saveE[k]=outE[k]+overlapF;
	  if(saveE[k]>nlines) saveE[k]=nlines;
	}
  }
  else 
  {
	for(k=0;k<numprocs;k++) saveS[k]=outS[k];
	for(k=0;k<numprocs;k++) saveE[k]=outE[k];
  }
 /* for(i=0;i<numprocs;i++)
  {
	fprintf(fp,"(1)i,saveS[i],saveE[i],outS[i],outE[i] %5d %5d %5d %5d %5d \n",i,saveS[i],saveE[i],outS[i],outE[i]); 
  fflush(fp);
  }
  fprintf(fp,"radonF,overlapX %7d %7d \n",radonF,overlapX);
  fflush(fp);*/
  if(radonF>0 && fxyF==0) 
  {
	for(k=0;k<numprocs;k++)
	{
	  saveS[k]=outS[k]-overlapX;
	  if(saveS[k]<0) saveS[k]=0;
	  saveE[k]=outE[k]+overlapX;
	  if(saveE[k]>nlines) saveE[k]=nlines;
	}
	if(radonF==3)
	{
	for(k=0;k<numprocs;k++)
	{
	  saveS[k]=outS[k]-overlapX*2;
	  if(saveS[k]<0) saveS[k]=0;
	  saveE[k]=outE[k]+overlapX*2;
	  if(saveE[k]>nlines) saveE[k]=nlines;
	}
	}
  }
  if(radonF>0 && fxyF>0) 
  {
	for(k=0;k<numprocs;k++)
	{
	  saveS[k]-=overlapX;
	  if(saveS[k]<0) saveS[k]=0;
	  saveE[k]+=overlapX;
	  if(saveE[k]>nlines) saveE[k]=nlines;
	}
	if(radonF==3)
	{
	for(k=0;k<numprocs;k++)
	{
	  saveS[k]-=overlapX;
	  if(saveS[k]<0) saveS[k]=0;
	  saveE[k]+=overlapX;
	  if(saveE[k]>nlines) saveE[k]=nlines;
	}
	}
  }
  saveE[0]=0;
  vmovi(saveS,saveS0,numprocs);
  vmovi(saveE,saveE0,numprocs);
  for(i=0;i<numprocs;i++)
//	fprintf(fp,"(2)i,saveS[i],saveE[i],outS[i],outE[i] %5d %5d %5d %5d %5d \n",i,saveS[i],saveE[i],outS[i],outE[i]); 
//  fflush(fp);

    numGroup=(nlines+(-overlapF+overlap1)*2)/(numDesiF-2*overlapF);
    if((nlines+(-overlapF+overlap1)*2)%(numDesiF-2*overlapF)>0) numGroup++;
    xGroup=(icrb+(-overlapF+overlap1)*2)/(numDesiF-2*overlapF);
    if((icrb+(-overlapF+overlap1)*2)%(numDesiF-2*overlapF)>0) xGroup++;

  if(fxyF>0 && myid>0)
  {
    i=saveE[myid]-saveS[myid]-overlapF*2;
	if(saveS[myid]==0 || myid==numprocs-1) i+=overlap1;
	if(saveS[myid]==0 && myid==numprocs-1) i+=overlap1;
	numGroup=i/(numDesiF-2*overlapF);
	if(i%(numDesiF-2*overlapF)>0) numGroup++;
    xGroup=(icrb+(-overlapF+overlap1)*2)/(numDesiF-2*overlapF);
    if((icrb+(-overlapF+overlap1)*2)%(numDesiF-2*overlapF)>0) xGroup++;
//    fprintf(fp,"numGroup,xGroup %7d %7d \n",numGroup,xGroup);
//    fflush(fp);
  }

  if(radonF>0 && myid>0)
  {
	if(fxyF==0)
	{
      i=saveE[myid]-saveS[myid]-overlapX*2;
	  if(saveS[myid]==0 || myid==numprocs-1) i+=overlap1;
	  if(saveS[myid]==0 && myid==numprocs-1) i+=overlap1;
	  numGroupR=i/(numDesiX-2*overlapX);
	  if(i%(numDesiX-2*overlapX)>0) numGroupR++;
	  numGroupR_s=numGroupR;
	}
//    numGroupR=(nlines+(-overlapX+overlap1)*2)/(numDesiX-2*overlapX);
//    if((nlines+(-overlapX+overlap1)*2)%(numDesiX-2*overlapX)>0) numGroupR++;
    xGroupR=(icrb+(-overlapI+overlap1)*2)/(numDesiI-2*overlapI);
    if((icrb+(-overlapI+overlap1)*2)%(numDesiI-2*overlapI)>0) xGroupR++;
//    fprintf(fp,"numGroupR,xGroupR %7d %7d \n",numGroupR,xGroupR);
//    fflush(fp);
  }

  get_freq();
  if(radonF) radonParm();
//  fprintf(fp,"lhwmax,low1,lhi1 %7d %7d %7d \n",lhwmax,low1,lhi1);

/* alloc arry */
  traceIn=FloatVector(0, shtotal1);
  trace = FloatVector(0, samples*2 );
  trace1 = FloatVector(0, samples*2+1000 );
  trace2 = FloatVector(0, samples*2+1000 );
  vclrf(trace,samples*2);
  vclrf(trace1,samples*2);
  vclrf(trace2,samples*2);
  ithdr = IntVector(0, headers+3 );
  rthdr = ( float * )ithdr;
  ithdr1 = IntVector(0, headers+3 );
  rthdr1 = ( float * )ithdr1;
  ithdr2 = IntVector(0, headers+3 );
  rthdr2 = ( float * )ithdr2;
  vclrf(rthdr1,headers+3);
  ithdr1[indexTrc_type]=2;
  mute = IntVector(0, idmap*2 );
  absR = FloatVector(0, idmap );
  vclrf(absR,idmap);

//	  fprintf(fp,"nlines,nbin %5d %5d %5d %5d \n",nlines,nbin,outS[myid],outE[myid]);
  if(myid==0 && nbin==1)
  {
    offDisk0=IntVector(0,nlines*nbin);
    vclri(offDisk0,nlines*nbin);
  }
  if(myid>0)
  {
    offDisk=IntVector(0,(saveE[myid]-saveS[myid])*nbin);
    vclri(offDisk,(saveE[myid]-saveS[myid])*nbin);
  }

  if(myid>0) 
  {
	diskOpen(); 
    if(radonF>0 && method==1) get_LLsave();
  }
// fprintf(fp,"done preFxy \n");
//       fflush(fp);
}

/* binStart */
void binStart()
{

// fprintf(fp,"binStart numline,lmin,lmax %6d %6d %6d \n",numline,lmin,lmax);
//       fflush(fp);
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0)
  {
    if(!idtyp) 
	{
	  bin=ithdr[indexpstmBin];
	}
	else bin=0;
    PrintTime();
    printf("start input trace for bin = %5d \n",bin);
  }
  MPI_Bcast(&bin,1,MPI_INT,0,MPI_COMM_WORLD);
  fxyDone=0;
  radonDone=0;
  oldcdp=-1;
  vclrf(absR,idmap);
  vclri(mute,idmap*2);

  vmovi(saveS0,saveS,numprocs);
  vmovi(saveE0,saveE,numprocs);
  if(myid==0 && nbin==1) vclri(offDisk0,nlines*nbin);
  if(myid>0 && nbin==1) vclri(offDisk,(saveE[myid]-saveS[myid])*nbin);
  if(FreBa)
  {
    vclrf(aveFreq,lhi0);
    vclrf(aveFout,lhi0);
    naveFreq=0;
  }
  numGroupR=numGroupR_s;
// fprintf(fp,"binStart done,bin,oldcdp  %7d %7d \n",bin,oldcdp);
//       fflush(fp);
}

void trace0Save()
{
  
/* fprintf(fp,"trace0Save,oldcdp,cdp %7d %7d \n",oldcdp,cdp);
       fflush(fp);*/
  cdp=idmap-1;
  while(cdp>=oldcdp+1)
  {
	oldcdp++;
	il=oldcdp/icrb;
	xl=oldcdp%icrb;
	traceDisk(oldcdp,0);        /* traces input to disk */
  }
}

void trace1Save()
{
  int i;
  
/* fprintf(fp,"trace1Save,oldcdp,cdp %7d %7d \n",oldcdp,cdp);
       fflush(fp);*/
  i=oldcdp;
  while(cdp>oldcdp+1)
  {
	oldcdp++;
	il=oldcdp/icrb;
	xl=oldcdp%icrb;
	traceDisk(oldcdp,0);        /* traces input to disk */
  }
  if(i<oldcdp)
  {
	il=cdp/icrb;
	xl=cdp%icrb;
  }
}

void freparm(int *fre)
{
  int i,j,k,nmaxfre,lh;
  float a1,maxfre;
  
  nmaxfre=*fre;
 //   fprintf(fp,"naveFreq aveFout %6d \n",naveFreq);
 /*   for(i=0;i<lhi0;i++)
    {
	  a1=(aveFout[i]);
      fprintf(fp,"%10.2f",a1);
      if(i%10==9) fprintf(fp,"i %8d \n",i);
    }
    fprintf(fp,"i %8d \n",i);
    fflush(fp);*/
	trace2[0]=(aveFreq[0]+aveFreq[1]+aveFreq[2])/3.;
	trace2[1]=(aveFreq[0]+aveFreq[1]+aveFreq[2]+aveFreq[3])*0.25;
	trace2[lhi0-1]=(aveFreq[lhi0-1]+aveFreq[lhi0-2]+aveFreq[lhi0-3])/3.;
	trace2[lhi0-2]=(aveFreq[lhi0-1]+aveFreq[lhi0-2]+aveFreq[lhi0-3]+aveFreq[lhi0-4])*0.25;
	for(i=2;i<lhi0-2;i++)
	  trace2[i]=(aveFreq[i-2]+aveFreq[i-1]+aveFreq[i]+aveFreq[i+1]+aveFreq[i+2])*0.2;
	vmovf(trace2,aveFreq,lhi0);
/*    fprintf(fp,"naveFreq aveFreq %6d \n",naveFreq);
    for(i=0;i<lhi0;i++)
    {
	  a1=(aveFreq[i]);
      fprintf(fp,"%10.2f",a1);
      if(i%10==9) fprintf(fp,"i %8d \n",i);
    }
    fprintf(fp,"i %8d \n",i);
    fflush(fp);*/
 /*   a1=lenft0/1000.*sampleRate;
    lh=(int)(qminp[0]*a1)+1;
	if(qminp[0]<7.) lh=(int)(7.*a1)+1; 
    else lh=(int)(qminp[0]*a1)+1;*/
	maxfre=0.;
	for(i=0;i<lhi0;i++) if(aveFreq[i]>maxfre) { maxfre=aveFreq[i]; nmaxfre=i; }
	a1=nmaxfre*FreParm;
      printf("a1 %10.2f \n",a1);
	  if(a1/(float)lhi0>0.25) a1=0.25*(float)lhi0;
	nmaxfre=(int)a1;
      printf("changed a1,nmaxfre %10.2f %7d \n",a1,nmaxfre);
	nmaxfre=(int)a1;
	trace2[0]=(aveFout[0]+aveFout[1]+aveFout[2])/3.;
	trace2[1]=(aveFout[0]+aveFout[1]+aveFout[2]+aveFout[3])*0.25;
	trace2[lhi0-1]=(aveFout[lhi0-1]+aveFout[lhi0-2]+aveFout[lhi0-3])/3.;
	trace2[lhi0-2]=(aveFout[lhi0-1]+aveFout[lhi0-2]+aveFout[lhi0-3]+aveFout[lhi0-4])*0.25;
	for(i=2;i<lhi0-2;i++)
	  trace2[i]=(aveFout[i-2]+aveFout[i-1]+aveFout[i]+aveFout[i+1]+aveFout[i+2])*0.2;
	vmovf(trace2,aveFout,lhi0);
/*    fprintf(fp,"naveFreq aveFout %6d \n",naveFreq);
    for(i=0;i<lhi0;i++)
    {
	  a1=(aveFout[i]);
      fprintf(fp,"%10.2f",a1);
      if(i%10==9) fprintf(fp,"i %8d \n",i);
    }
    fprintf(fp,"i %8d \n",i);
    fflush(fp);*/
	for(i=0;i<nmaxfre;i++) 
	{
	  if(aveFout[i]>0.)	
	  {
	    aveFreq[i]=(aveFreq[i]*aveFout[nmaxfre]/(aveFreq[nmaxfre]*aveFout[i])); 
//		if(FreBa==2) aveFreq[i]=sqrt(aveFreq[i]);
	  }
	  else aveFreq[i]=0.;
	}

	*fre=nmaxfre;
/*    fprintf(fp,"nmaxfre aveFreq %6d \n",nmaxfre);
    for(i=0;i<lhi0;i++)
    {
	  a1=(aveFreq[i]);
      fprintf(fp,"%10.2f",a1);
      if(i%10==9) fprintf(fp,"i %8d \n",i);
    }
    fprintf(fp,"i %8d \n",i);
    fflush(fp);*/
}

void  get_out()
{
  int i,j,k,m,seq,i1,k1,i2,il;
  float a1;

  MPI_Barrier(MPI_COMM_WORLD);
 //   fprintf(fp,"traceIn[0] %8.1f \n",traceIn[0]);
 //   fflush(fp);
  if(myid==0)
  {
	PrintTime();
    printf("start output traces by CDP order \n");
    masterOut();
  }
  else slaveOut();
}

void  masterOut()
{
  int i,j,k,k1,seq,il,sender,i2,kkb,kk[100],num,num0;
  float a1,a2,*buf,*buf1;

	buf = FloatVector(0,shtotal*icrb*2 );
	buf1 = FloatVector(0,shtotal*20*2 );
//  num0=0;
  a2=0.;
  PrintTime();
  printf("finished output %6.1f percent of inlines \n",a2);
  a2=10.0;
  j=0;
  for(i=1;i<numprocs;i++)
  {
    {
	  MPI_Send(&i,1,MPI_INT,i,1,MPI_COMM_WORLD);
//  fprintf(fp,"send i %5d \n",i);
//  fflush(fp);
      for(il=outS[i];il<outE[i];il++)
      {
	    if((float)(il*100.)/(float)nlines>a2) 
	    {
          PrintTime();
		  printf("finished output %6.1f percent of inlines \n",a2);
		  a2+=10.;
		}
	    MPI_Recv(&num,1,MPI_INT,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
//  fprintf(fp,"il,get num %5d %5d \n",il,num);
//  fflush(fp);
	    if(num==0) continue;
		while(num>0)
		{
		  if(num>icrb) { k=icrb; num-=icrb;}
		  else { k=num; num=0;}
	      MPI_Recv(buf,icrb*shtotal,MPI_FLOAT,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
//  fprintf(fp,"recv (1) k,num,j %6d %6d %6d \n",k,num,j );
//  fflush(fp);
			for(k1=0;k1<k;k1++)
			{
			  vmovf(&buf[k1*shtotal],&buf1[j*shtotal],shtotal);
			  j++;
		      if(j>10)
		      {
//  fprintf(fp,"out,j,num0 %6d %6d \n",j,num0 );
//  fflush(fp);
                write( socketFile, ( char* )buf1, shtotal4 );
		        vmovf(&buf1[shtotal1],buf1,(j-10)*shtotal);
		        j-=10;
//				num0++;
		      }
			}
		}
	  }
	}
  }
  buf1[(j-1)*shtotal]=-1.;
//  fprintf(fp,"j,k,buf1[(j-1)*shtotal] %5d %5d %8.1f \n",j,k,buf1[(j-1)*shtotal]);
//  fflush(fp);
  while(j>0) 
  {
//	fprintf(fp,"last j,num0 %5d %5d \n",j,num0);
//	fflush(fp);
//	num0++;
    write( socketFile, ( char* )buf1, shtotal4 );
	j-=10;
//	fprintf(fp,"last j,num0 %5d %5d \n",j,num0);
//	fflush(fp);
	if(j>0) vmovf(&buf1[shtotal1],buf1,(j-10)*shtotal);
  }
}

void Master0(float *buff,cdpbin *Cdpbin,int *kkb,int *kk,int il)
{
  int i,i1,j,k,bin,cdp,oldcdp,seqno;

// fprintf(fp,"Master0 il %5d \n",il);
//       fflush(fp);
// fprintf(fp,"sender,outE[myid],outS[myid] %7d %7d %7d \n",sender,outE[myid],outS[myid]);
//  for(i=1;i<outE[myid]-outS[myid];i++) offDisk[i]+=offDisk[i-1];
  *kkb=0;
  for(bin=0;bin<nbin;bin++)
  {
	i1=bin*(outE[0]-outS[0])+il;
	if(offDisk[i1]==0) continue;
    get_bin(&buff[(*kkb)*shtotal],kk[bin],offDisk[i1],bin);
	for(j=0;j<offDisk[i1];j++)
	{
      k=j*shtotal;
      vmovf(&buff[(*kkb+j)*shtotal+1],rthdr,headers);
	  Cdpbin[(*kkb+j)].seq=*kkb+j;
	  Cdpbin[(*kkb+j)].bin=bin;
	  cdp=ithdr[indexCdp];
	  Cdpbin[(*kkb+j)].cdp=cdp;
	}
    *kkb+=offDisk[i1];
	kk[bin]+=offDisk[i1];
  }
  if(*kkb==0) return; 
  qsort(Cdpbin,*kkb,sizeof(cdpbin),compare);  /* cdp,bin order */
  oldcdp=0;
  for(j=0;j<*kkb;j++)
  {
	i1=Cdpbin[j].seq*shtotal;
	vmovf(&buff[i1+1],rthdr,headers);
	cdp=ithdr[indexCdp];
	if(oldcdp<cdp)
	{
	  ithdr[indexSeqno]=1;
	  ithdr[indexEnd_ens]=0;
	  ithdr2[indexEnd_ens]=1;
	  if(oldcdp>0)
	  {
		  vmovf(rthdr2,&buff[Cdpbin[j-1].seq*shtotal+1],headers);
	  }
	  seqno=2;
	  oldcdp=cdp;
	}
	else ithdr[indexSeqno]=seqno++;
	if(j==*kkb-1) 
	{
	  ithdr[indexEnd_ens]=1;
	}
	vmovf(rthdr,&buff[i1+1],headers);
	vmovf(rthdr,rthdr2,headers);
  }
}

void  slaveOut()
{
  int i,i1,j,il,k,sender,bin,kk[100],kkb,cdp,oldcdp,seqno;
  float *buf,*buf1;

// fprintf(fp,"slaveOut \n");
//       fflush(fp);
  free( buff);
  free( buff1);
  buf=FloatVector(0,shtotal*icrb*nbin);
  buf1=FloatVector(0,shtotal*icrb*2);
  Cdpbin =(cdpbin *)malloc(icrb*nbin*sizeof( cdpbin ) );
/*  for(bin=0;bin<nbin;bin++)
  {
	for(i=0;i<outE[myid]-outS[myid];i++)
	{
	  fprintf(fp,"bin,i,offDisk[bin*(outE[myid]-outS[myid])+i] %5d %5d %5d \n",bin,i,offDisk[bin*(outE[myid]-outS[myid])+i]);
	}
	fflush(fp);
  }*/
  MPI_Recv(&sender,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
// fprintf(fp,"sender,outE[myid],outS[myid] %7d %7d %7d \n",sender,outE[myid],outS[myid]);
//  for(i=1;i<outE[myid]-outS[myid];i++) offDisk[i]+=offDisk[i-1];
  vclri(kk,nbin);
  for(il=0;il<outE[myid]-outS[myid];il++)
  {
	kkb=0;
	for(bin=0;bin<nbin;bin++)
	{
	  i1=bin*(outE[myid]-outS[myid])+il;
	  if(offDisk[i1]==0) continue;
// fprintf(fp,"il,bin,kkb,i1,offDisk[i1] %7d %7d %7d %7d %7d \n",il,bin,kkb,i1,offDisk[i1]);
// fflush(fp);
      get_bin(&buf[kkb*shtotal],kk[bin],offDisk[i1],bin);
	  for(j=0;j<offDisk[i1];j++)
	  {
		vmovf(&buf[(kkb+j)*shtotal+1],rthdr,headers);
		Cdpbin[(kkb+j)].seq=kkb+j;
		Cdpbin[(kkb+j)].bin=bin;
		cdp=ithdr[indexCdp];
		Cdpbin[(kkb+j)].cdp=cdp;
 // fprintf(fp,"j,Cdpbin[kkb+j].seq,Cdpbin[kkb+j].cdp,Cdpbin[kkb+j].bin %7d %7d %7d %7d \n",j,Cdpbin[kkb+j].seq,Cdpbin[kkb+j].cdp,Cdpbin[kkb+j].bin);
// fflush(fp);
	  }
      kkb+=offDisk[i1];
	  kk[bin]+=offDisk[i1];
// fprintf(fp,"kkb,il,bin %7d %7d %7d \n",kkb,il,bin);
// fflush(fp);
	}
    if(kkb>0) qsort(Cdpbin,kkb,sizeof(cdpbin),compare);  /* cdp,bin order */
    MPI_Ssend(&kkb,1,MPI_INT,0,1,MPI_COMM_WORLD);
	if(kkb==0) continue;
    oldcdp=0;
	for(j=0;j<kkb;j++)
	{
	  i1=Cdpbin[j].seq*shtotal;
	  vmovf(&buf[i1+1],rthdr,headers);
	  cdp=ithdr[indexCdp];
//  fprintf(fp,"j,cdp %7d %7d  \n",j,cdp);
// fflush(fp);
	  if(oldcdp<cdp)
	  {
	    ithdr[indexSeqno]=1;
	    ithdr[indexEnd_ens]=0;
	    ithdr2[indexEnd_ens]=1;
		if(oldcdp>0)
		{
		  vmovf(rthdr2,&buf[Cdpbin[j-1].seq*shtotal+1],headers);
		}
		seqno=2;
		oldcdp=cdp;
	  }
	  else ithdr[indexSeqno]=seqno++;
	  if(j==kkb-1) 
	  {
	    ithdr[indexEnd_ens]=1;
	  }
	  vmovf(rthdr,&buf[i1+1],headers);
	  vmovf(rthdr,rthdr2,headers);
	}
	k=0;
//  fprintf(fp,"start send kkb %7d \n",kkb);
// fflush(fp);
	for(j=0;j<kkb;j++)
	{
	  vmovf(&buf[Cdpbin[j].seq*shtotal],&buf1[k*shtotal],shtotal);
	  k++;
	  if(k==icrb)
	  {
        MPI_Ssend(buf1,shtotal*icrb,MPI_FLOAT,0,1,MPI_COMM_WORLD);
		k=0;
	  }
	}
	if(k>0) MPI_Ssend(buf1,shtotal*icrb,MPI_FLOAT,0,1,MPI_COMM_WORLD);
  }
}

int compare( const void *a1, const void *a2 )
{
   cdpbin *a3,*a4;
   a3=(cdpbin *)a1; a4=(cdpbin *)a2;           
   if(a3->cdp == a4->cdp) return(a3->bin - a4->bin);
   else return(a3->cdp - a4->cdp );
}   

void Abalanc(float *trace,int xl)
{
  int i,j,k,ii,am,ai;
  float amax,amin,a1,ama,ami;

  ii=xl*sampleAs;
  amin=buff1[ii];
  amax=buff1[ii+1];
  if(amin>0.) a1=amax/amin;
  else a1=0.;
  vmovf(&buff1[ii],trace2,sampleAs);

  for(i=2;i<sampleAs;i++) 
  {
	if(trace2[i]<tolAmin && trace2[i]>0.) trace2[i]=tolAmin;
	if(trace2[i]>tolAmax && trace2[i]>0.) trace2[i]=tolAmax;
  }
  trace2[0]=trace2[1]=trace2[2];

  for(i=0;i<sampleAs;i++) if(trace2[i]/tolAmin>AmpParm) trace2[i]=AmpParm*tolAmin;

  for(i=0;i<sampleAs;i++) if(trace2[i]>0.) {j=i; break;}
  if(i<sampleAs)
  {
    for(i=0;i<j;i++) trace2[i]=trace2[j];
    for(i=sampleAs-1;i>=0;i--) if(trace2[i]>0.) {j=i; break;}
    if(j<sampleAs-1) for(i=j+1;i<sampleAs;i++) trace2[i]=trace2[j];
  }

  trace[samples-1]*=trace2[(samples-1)/gap];
  for(k=0;k<samples-1;k++) 
  {
    if(k/gap<=2) trace[k]*=trace2[2];
	else
	{
	  if(k%gap==gap/2) trace[k]*=(trace2[k/gap]+trace2[k/gap+1])*0.5;
	  if(k%gap==0) trace[k]*=trace2[k/gap];
	  if(gap>2)
	  {
        if(k%gap==1) trace[k]*=trace2[k/gap]*0.75+trace2[k/gap+1]*0.25;
        if(k%gap==3) trace[k]*=trace2[k/gap]*0.25+trace2[k/gap+1]*0.75;
	  }
    }
  }
//  fprintf(fp,"il,xl,amin,amax,amax/amin %6d %6d %10.4f %10.4f %10.2f \n",il,xl,amin,amax,a1);
//  fflush(fp);
}

void  get_out_one()
{
  int i,j,k,m,seq,nmaxfre,lh,lh01,mm;
  float a1,a2,a3;
  int i2,cdp2;
  float b1;

  MPI_Barrier(MPI_COMM_WORLD);
// fprintf(fp,"get_out_one \n");
//       fflush(fp);
  if(myid==0) { masterOne();}
  else slaveOne();
}

void  masterOne()
{
  int i,j,k,m,seq,il,k1,sender,kk,lh,lh01,nmaxfre,mm;
  float a1;

    PrintTime();
    printf("start output traces for bin = %5d \n",bin);
	mm=shtotal*icrb*2;
	  buff = FloatVector(0,shtotal*icrb*4 );

    for(i=1;i<numprocs;i++)
	{
//	  fprintf(fp,"i,outS[i],outE[i] %5d %5d %5d \n",i,outS[i],outE[i]);
//	fflush(fp);
	  MPI_Recv(&offDisk0[outS[i]],outE[i]-outS[i],MPI_INT,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
//	  fprintf(fp,"i,outS[i],offDisk0[outS[i]] %5d %5d %5d \n",i,outS[i],offDisk0[outS[i]]);
//	fflush(fp);
	}
/*	for(i=0;i<nlines;i++)
	{
	  fprintf(fp,"i,offDisk0[i] %5d %5d \n",i,offDisk0[i]);
	}
	fflush(fp);*/

	j=0;
	seq=1;
	for(i=1;i<numprocs;i++)
	{
	  MPI_Send(&i,1,MPI_INT,i,1,MPI_COMM_WORLD);
      for(il=outS[i];il<outE[i];il++)
      {
	    if(offDisk0[il]==0) continue;
		  if(offDisk0[il]*shtotal<1000000) 
	      MPI_Recv(buff,offDisk0[il]*shtotal,MPI_FLOAT,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		  else
		  {
			k=0; m=offDisk0[il]/5;
			while(k<5)
			{
			  if(k<4) 
	          MPI_Recv(&buff[k*m*shtotal],m*shtotal,MPI_FLOAT,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			  else
	          MPI_Recv(&buff[k*m*shtotal],(offDisk0[il]-4*m)*shtotal,MPI_FLOAT,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			  k++;
			}
		  }
//	      MPI_Recv(buff,offDisk0[il]*shtotal,MPI_FLOAT,i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	    for(k=0;k<offDisk0[il];k++)
	    {
		  k1=k*shtotal;
          if(i==0 || !idtyp) vmovf(&buff[k1+1],rthdr2,headers);
		  if(!idtyp) 
		  {
		    ithdr2[indexSeqno]=seq++;
            vmovf(rthdr2,&buff[k1+1],headers);
		  }
          vmovf(&buff[k1],&buff[mm+j*shtotal],shtotal);
		  j++;
		  if(j>10)
		  {
            write( socketFile, ( char* )&buff[mm], shtotal4 );
		    vmovf(&buff[mm+shtotal1],&buff[mm],(j-10)*shtotal);
		    j-=10;
		  }
		}
	  }
	}
	buff[mm+(j-1)*shtotal]=-1.;
	if(!idtyp)
	{
       vmovf(&buff[mm+(j-1)*shtotal+1],rthdr2,headers);
	   ithdr2[indexEnd_ens]=1;
       vmovf(rthdr2,&buff[mm+(j-1)*shtotal+1],headers);
	}
    while(j>0) 
	{
//	  fprintf(fp,"last j %5d \n",j);
//	fflush(fp);
      write( socketFile, ( char* )&buff[mm], shtotal4 );
	  j-=10;
	  if(j>0) vmovf(&buff[mm+shtotal1],&buff[mm],(j-10)*shtotal);
	}
}

void  slaveOne()
{
  int i,k,il,k1,sender,kk,nmaxfre,lh,lh01,m,low1;
  float a1,a2,a3;
  int i2,cdp2,cdp,il1,xl1;
  float b1;

// fprintf(fp,"slaveOne \n");
//       fflush(fp);
  if(myid>0)
  {
	MPI_Ssend(offDisk,outE[myid]-outS[myid],MPI_INT,0,1,MPI_COMM_WORLD);
/*	for(i=0;i<outE[myid]-outS[myid];i++)
	{
	  fprintf(fp,"i,offDisk[i] %5d %5d \n",i,offDisk[i]);
	}
	fflush(fp);*/
    MPI_Recv(&sender,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
  }
// fprintf(fp,"sender,outE[myid],outS[myid] %7d %7d %7d \n",sender,outE[myid],outS[myid]);
//  for(i=1;i<outE[myid]-outS[myid];i++) offDisk[i]+=offDisk[i-1];
  else vmovi(offDisk,offDisk0,outE[0]-outS[0]);
  if(FreBa>0)
  {
	i2=(int)(100./sampleRate);
	b1=1./(float)i2;
    a1=lenft0/1000.*sampleRate;
    lh=(int)(qminp[0]*a1)+1;
    low1=(int)(qmin[0]*a1);
    lh01=(int)(qmaxm[0]*a1)+1;
    a2=1./((float)(lhi0-lh01));
	freparm(&nmaxfre);
// fprintf(fp,"i2,lh,lh01,nmaxfre,a1,a2,b1 %5d %5d %5d %5d %7.4f %7.4f %7.4f \n",i2,lh,lh01,nmaxfre,a1,a2,b1);
// fflush(fp);
  }
  kk=0;
  for(il=0;il<outE[myid]-outS[myid];il++)
  {
	if(offDisk[il]==0) continue;
    get_bin(buff,kk,offDisk[il],0);
/*		if(il==10)
		{
	    for(m=400;m<500;m++)
	    {
	      fprintf(fp,"%9.1f",buff[shtotal*10+1+headers+m]);
	      if(m%10==9) fprintf(fp," %8d \n",m);
	    }
	    fprintf(fp,"%6d %6d \n",m,il);
		}*/

    if(AmpBa) get_ampl(buff1,il*icrb,icrb);
	for(k=0;k<offDisk[il];k++)
	{
	  k1=k*shtotal+1+headers;
	  m=(int)buff[k1-1-headers]%icrb;
	  if(AmpBa) Abalanc(&buff[k1],m);
	  if(!FreBa)
	  {
	    continue;
	  }
	  else
	  {
        vmovf(&buff[k1-headers],rthdr2,headers);
		vclrf(trace1,lenft0+10);
	    vmovf(&buff[k1],trace1,samples);
        ForwardPrimeFFT( trace1, lenft0 );
		vclrf(&trace1[lhi0*2],lenft0+10-lhi0*2);

        for(m=lh01;m<lhi0;m++)
		{
			  a3=a2*(lhi0-m);
		      trace1[m+m]*=a3;
		      trace1[m+m+1]*=a3;
		}
		vclrf(trace1,low1+low1);
       for(m=low1;m<lh;m++)
		{
//			  a3=(m/(float)lh+1.0)*0.5;
			  a3=(m-low1+1)/(float)(lh-low1+1);
		      trace1[m+m]*=a3;
		      trace1[m+m+1]*=a3;
		}
        for(m=low1;m<nmaxfre;m++)
		{
			if(aveFreq[m]>1.)
			{
		      trace1[m+m]*=(aveFreq[m]-1.)*1.2+1.;
		      trace1[m+m+1]*=(aveFreq[m]-1.)*1.2+1.;
			}
		}
        ReversePrimeFFT( trace1, lenft0 );
//	    if(muteop)
	    {
		  cdp=ithdr2[indexCdp];
		  cdp-=minCDP;
          il1=cdp/ncrbi+minIlin-lmin;
          xl1=cdp%ncrbi+minXlin-xmin;
          cdp=il1*icrb+xl1;
	      cdp2=cdp+cdp;
	      if(mute[cdp2]>0) vclrf(trace1,mute[cdp2]);
		  if(mute[cdp2+1]<samples-1) 
			vclrf(&trace1[mute[cdp2+1]+1],samples-mute[cdp2+1]-1);
		  for(m=mute[cdp2];m<mute[cdp2]+i2;m++)
            trace1[m]*=b1*(m-mute[cdp2]);
	    }
	    vmovf(trace1,&buff[k1],samples);
	  }
	}
	if(myid>0) 
	{
	  if(offDisk[il]*shtotal<1000000) 
	  MPI_Ssend(buff,shtotal*offDisk[il],MPI_FLOAT,0,1,MPI_COMM_WORLD);
	  else
	  {
		k=0; m=offDisk[il]/5;
		while(k<5)
		{
		  if(k<4) 
	      MPI_Ssend(&buff[k*m*shtotal],m*shtotal,MPI_FLOAT,0,1,MPI_COMM_WORLD);
		  else
	      MPI_Ssend(&buff[k*m*shtotal],(offDisk[il]-4*m)*shtotal,MPI_FLOAT,0,1,MPI_COMM_WORLD);
		  k++;
		}
	  }
    }
//	if(myid>0) MPI_Ssend(buff,shtotal*offDisk[il],MPI_FLOAT,0,1,MPI_COMM_WORLD);
    else put_bin(buff,kk,offDisk[il],0);
	kk+=offDisk[il];
  }
}

void Get_aveFreq()
{
  int i,nn;
  float a1;

  nn=0;
  vclrf(trace2,lhi0);
  if(myid>0)
	for(i=0;i<lhi0;i++) aveFreq[i]*=naveFreq;
  MPI_Reduce(&naveFreq,&nn,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(aveFreq,trace2,lhi0,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Bcast(&nn,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(trace2,lhi0,MPI_FLOAT,0,MPI_COMM_WORLD);

//  fprintf(fp,"nn,naveFreq %5d %5d \n",nn,naveFreq);
  a1=1./(float)nn;
  for(i=0;i<lhi0;i++) aveFreq[i]=trace2[i]*a1;
/*  for(i=0;i<lhi0;i++)
  {
	fprintf(fp,"%9.1f ",aveFreq[i]);
	if(i%10==9) fprintf(fp,"  %6d \n",i);
  }
  fprintf(fp,"  %6d \n",i);
  fflush(fp);*/
}

void Get_aveFout()
{
  int i,nn;
  float a1;

  nn=0;
  vclrf(trace2,lhi0);
  if(myid>0)
	for(i=0;i<lhi0;i++) aveFout[i]*=naveFreq;
  MPI_Reduce(&naveFreq,&nn,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(aveFout,trace2,lhi0,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Bcast(&nn,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(trace2,lhi0,MPI_FLOAT,0,MPI_COMM_WORLD);
  a1=1./(float)nn;
  for(i=0;i<lhi0;i++) aveFout[i]=trace2[i]*a1;

//  fprintf(fp,"\n nn,naveFreq %5d %5d \n",nn,naveFreq);
  a1=1./(float)nn;
  for(i=0;i<lhi0;i++) aveFout[i]=trace2[i]*a1;
/*  for(i=0;i<lhi0;i++)
  {
	fprintf(fp,"%9.1f ",aveFout[i]);
	if(i%10==9) fprintf(fp,"  %6d \n",i);
  }
  fprintf(fp,"  %6d \n",i);
  fflush(fp);*/
}

void Get_tol()
{
  int i,nn;
  float a1,a2;

  nn=0;
  a1=a2=0.;
//	fprintf(fp,"ntolA,nn,a1,a2,tolAmin,tolAmax %5d %5d %11.1f %11.1f %9.1f %9.1f \n",ntolA,nn,a1,a2,tolAmin,tolAmax);
//	fflush(fp);
  MPI_Reduce(&ntolA,&nn,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&tolAmin,&a1,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&tolAmax,&a2,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Bcast(&nn,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&a1,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  MPI_Bcast(&a2,1,MPI_FLOAT,0,MPI_COMM_WORLD);
  tolAmin=a1/(float)nn;
  tolAmax=a2/(float)nn;
//	fprintf(fp,"ntolA,nn,a1,a2,tolAmin,tolAmax %5d %5d %11.1f %11.1f %9.1f %9.1f \n",ntolA,nn,a1,a2,tolAmin,tolAmax);
//	fflush(fp);
}

void FreAmp()
{
  int i,ii,k,il,k1,kk,nmaxfre,lh,lh01,m,low1;
  float a1,a2,a3;
  int i2,cdp2,cdp,il1,xl1;
  float b1;

//	fprintf(fp,"FreAmp ntolA,tolAmin,tolAmax,bin %5d %9.1f %9.1f %6d \n",ntolA,tolAmin,tolAmax,bin);
//	fflush(fp);
  if(FreBa)
  {
	i2=(int)(100./sampleRate);
	b1=1./(float)i2;
    a1=lenft0/1000.*sampleRate;
    lh=(int)(qminp[0]*a1)+1;
    low1=(int)(qmin[0]*a1);
    lh01=(int)(qmaxm[0]*a1)+1;
    a2=1./((float)(lhi0-lh01));
	freparm(&nmaxfre);
  }

//	fprintf(fp,"FreAmp ntolA,tolAmin,tolAmax,bin,nmaxfre %5d %9.1f %9.1f %6d %6d \n",ntolA,tolAmin,tolAmax,bin,nmaxfre);
//	fflush(fp);
  kk=0;
  for(il=0;il<outE[myid]-outS[myid];il++)
  {
	ii=(bin-1)*(outE[myid]-outS[myid])+il;
//	fprintf(fp,"il,kk,bin,offDisk[ii] %5d %6d %6d %6d \n",il,kk,bin,offDisk[ii]);
//	fflush(fp);
	if(offDisk[ii]==0) continue;
    get_bin(buff,kk,offDisk[ii],bin-1);
/*		if(il==10)
		{
	    for(m=400;m<500;m++)
	    {
	      fprintf(fp,"%9.1f",buff[shtotal*10+1+headers+m]);
	      if(m%10==9) fprintf(fp," %8d \n",m);
	    }
	    fprintf(fp,"%6d %6d \n",m,il);
		}*/

    if(AmpBa) get_ampl(buff1,il*icrb,icrb);
	for(k=0;k<offDisk[ii];k++)
	{
	  k1=k*shtotal+1+headers;
	  m=(int)buff[k1-1-headers]%icrb;
	  if(AmpBa) Abalanc(&buff[k1],m);
	  if(!FreBa)
	  {
	    continue;
	  }
	  else
	  {
        vmovf(&buff[k1-headers],rthdr2,headers);
		vclrf(trace1,lenft0+10);
	    vmovf(&buff[k1],trace1,samples);
        ForwardPrimeFFT( trace1, lenft0 );
		vclrf(&trace1[lhi0*2],lenft0+10-lhi0*2);

        for(m=lh01;m<lhi0;m++)
		{
			  a3=a2*(lhi0-m);
		      trace1[m+m]*=a3;
		      trace1[m+m+1]*=a3;
		}
		vclrf(trace1,low1+low1);
       for(m=low1;m<lh;m++)
		{
//			  a3=(m/(float)lh+1.0)*0.5;
			  a3=(m-low1+1)/(float)(lh-low1+1);
		      trace1[m+m]*=a3;
		      trace1[m+m+1]*=a3;
		}
        for(m=low1;m<nmaxfre;m++)
		{
			if(aveFreq[m]>1.)
			{
		      trace1[m+m]*=(aveFreq[m]-1.)*1.2+1.;
		      trace1[m+m+1]*=(aveFreq[m]-1.)*1.2+1.;
			}
		}
        ReversePrimeFFT( trace1, lenft0 );
//	    if(muteop)
	    {
		  cdp=ithdr2[indexCdp];
		  cdp-=minCDP;
          il1=cdp/ncrbi+minIlin-lmin;
          xl1=cdp%ncrbi+minXlin-xmin;
          cdp=il1*icrb+xl1;
	      cdp2=cdp+cdp;
	      if(mute[cdp2]>0) vclrf(trace1,mute[cdp2]);
		  if(mute[cdp2+1]<samples-1) 
			vclrf(&trace1[mute[cdp2+1]+1],samples-mute[cdp2+1]-1);
		  for(m=mute[cdp2];m<mute[cdp2]+i2;m++)
            trace1[m]*=b1*(m-mute[cdp2]);
	    }
	    vmovf(trace1,&buff[k1],samples);
	  }
	}
    put_bin(buff,kk,offDisk[ii],bin-1);
	kk+=offDisk[ii];
  }
//	fprintf(fp,"FreAmp ntolA,tolAmin,tolAmax,bin %5d %9.1f %9.1f %6d \n",ntolA,tolAmin,tolAmax,bin);
//	fflush(fp);
}

void Fxy()
{
  int i,k,seq;

  if(idtyp) binStart();
  while( 1 )                                 
  {
	if(myid==0) GetSocketData( ( char* )traceIn, shtotal4);
    MPI_Bcast(traceIn,shtotal1,MPI_FLOAT,0,MPI_COMM_WORLD);
// if(traceIn[0]!=1.) fprintf(fp,"traceIn[0] %7.1f \n",traceIn[0]);
  //     fflush(fp);
	if(traceIn[0]==-3.)  //all done
	{
	  if(nbin>1) 
	  {
		get_out();
		break;
	  }
	  else break;
	}
	if(traceIn[0]==-2.) 
	{
	  if(oldcdp<idmap-1) trace0Save();
 //fprintf(fp,"traceIn[0]=-2 \n");
  //     fflush(fp);
	  if(FreBa) Get_aveFreq();
	  Perform();
	  if(FreBa) Get_aveFout();
	  if(AmpBa) Get_tol();
      MPI_Barrier(MPI_COMM_WORLD);
	  if(nbin==1) 
	  {
		get_out_one();
		if(idtyp) break;
		else continue;
	  }
	  else 
		if((FreBa || AmpBa) && myid>0) 
	    {
		  FreAmp();
//		  if(bin==nbin) { get_out(); break; }
	    }
	}
	for(i=0;i<10;i++)
	{
	  k=i*shtotal;
	  vmovf(&traceIn[k+1],rthdr,headers);
	  vclrf(&trace[samples],200);
	  vmovf(&traceIn[k+headers+1],trace,samples);
      cdp=ithdr[indexCdp]-minCDP;                   /* no of cdp */
      il=cdp/ncrbi+minIlin;
      xl=cdp%ncrbi+minXlin;
	  if(k==0 && !idtyp)
	  {
		seq=ithdr[indexSeqno];
        if(seq==1) binStart();
	  }
      il=il-lmin;
      xl=xl-xmin;  
      cdp=il*icrb+xl;
	  trace1Save();
 //  fprintf(fp,"fxy cdp,il,xl %7d %7d %7d  \n",cdp,il,xl);
 //  fflush(fp);
	  traceDisk(cdp,1);               /* traces input to disk */
	  oldcdp=cdp;
// fprintf(fp,"i,traceIn[k] %5d %7.1f \n",i,traceIn[k]);
//       fflush(fp);
      if(traceIn[k]<0.)
	  {
		trace0Save();
	    if(FreBa) Get_aveFreq();
		Perform();
	    if(FreBa) Get_aveFout();
	    if(AmpBa) Get_tol();
        MPI_Barrier(MPI_COMM_WORLD);
		if(nbin==1) get_out_one();
		else 
		if((FreBa || AmpBa) && myid>0) 
	    {
		  FreAmp();
//		  if(bin==nbin) { get_out(); break; }
	    }
	    break;
      }
    }
//if(i!=10 || traceIn[k]!=1.) fprintf(fp,"i,traceIn[k] %5d %7.1f \n",i,traceIn[k]);
  //     fflush(fp);
	if(i<10 && idtyp) break;
//	if(i<10 && !idtyp && bin==nbin) { get_out(); break; }
  }
}

//trace input trace2 amplitude
void getAmpl()   
{
  int i,j,j1,j2,k,nOper,nOp2;
  float sum,a1;

  nOper=(int)(AmpLong/sampleRate);
  nOp2=nOper/2;
  nOper=nOp2*2+1;
  sum=0.; k=0;
  for(i=0;i<=nOp2;i++) 
  {
	a1=fabs(trace[i]);
	if(a1>0.)
	{
	  sum+=a1;
	  k++;
	}
  }
  if(k>0) trace2[0]=sum/(float)k;
  else trace2[0]=0.;
  for(i=gap;i<samples;i+=gap)
  {
	j1=i-nOp2-gap;
	j2=i+nOp2+1-gap;
	for(j=0;j<gap;j++)
	{
	  if(j1+j>0)
	  {
	    a1=fabs(trace[j1+j]);
	    if(a1>0.) { sum-=a1; k--; }
	  }
	  if(j2+j<samples)
	  {
	    a1=fabs(trace[j2+j]);
	    if(a1>0.) { sum+=a1; k++; }
	  }
	}
    if(k>0) trace2[i/gap]=sum/(float)k;
    else trace2[i/gap]=0.;
  }
/*  if(xl==10) 
  {
	fprintf(fp,"(1)il,nOper,nOp2,k,sum %5d %5d %5d %5d %10.1f \n",il,nOper,nOp2,k,sum);
	fprintf(fp,"i,j1,j2,k,sum,trace2[i/gap] %5d %5d %5d %5d %10.1f %10.1f \n",i,j1,j2,k,sum,trace2[i/gap]);
    fflush(fp);
  }
  if(xl==2) exit(0);*/
}

void traceDisk(int cdp1,int kk)
{
  int i,j,k,il1,xl1,cdpp;
  float a1,a2;
  static saveTrace=0;

  il1=cdp1/icrb;
  if(il1<saveS[myid] || il1>=saveE[myid]) return;
//   fprintf(fp,"traceDisk cdp1,il1,kk,saveTrace,saveHeader %7d %7d %7d %7d %7d  \n",cdp1,il1,kk,saveTrace,saveHeader);
//   fflush(fp);
  if(kk==0)
  {
	if(il1>=outS[myid] && il1<outE[myid]) 
	{
	  il1+=lmin;
      xl1=cdp1%icrb+xmin;
      cdpp=(il1-minIlin)*ncrbi+xl1-minXlin+minCDP;
	  ithdr1[indexCdp]=cdpp;
      vmovi(ithdr1,&ihead[saveTrace*headers],headers);
	  if(AmpBa) vclrf(&ampl[sampleAs*saveTrace],sampleAs);
	}
	vclrf(&buff[freTrace*saveTrace],freTrace);
  }
  else
  {
    for(i=0;i<samples;i++) if(trace[i]!=0.) break;
    mute[cdp1*2]=i;  
    if(i!=samples)        
    {
      for(j=samples-1;j>=i;j--) if(trace[j]!=0.) break;
      mute[cdp1*2+1]=j; 
	  for(k=i;k<=j;k++) 
	  {
		a1=fabs(trace[k]);
		if(a1>10e+8) 
		{
		  trace[k]=0.;
		  continue;
		}
		absR[cdp1]+=fabs(trace[k]);
	  }
	  if(k<j) 
	  {
        vclrf(trace,samples);
		ithdr[indexTrc_type]=0;
        absR[cdp1]=0.;
        mute[cdp1*2]=samples;  
        mute[cdp1*2+1]=0;  
        if(il1>=outS[myid] && il1<outE[myid] && AmpBa)
	       vclrf(&ampl[sampleAs*saveTrace],sampleAs);
        printf("trace killed cdp = %7d \n",ithdr[indexCdp]);
	  }
	  else
	  {
		if(il1>=outS[myid] && il1<outE[myid] && AmpBa)
	    {
		  getAmpl();
          vmovf(trace2,&ampl[sampleAs*saveTrace],sampleAs);
	    }
	  }
    }

    if(i<=(int)(timee/sampleRate)) 
    {
	  for(i=0;i<numWin;i++) 
	  {
	    vmovf(&trace[i*(sampleWin-sampleTap)],trace1,sampleWin);
	    vclrf(&trace1[sampleWin],lenft+10-sampleWin);
        ForwardPrimeFFT( trace1, lenft );
	    vmovf(&trace1[2*low[i]],&buff[saveTrace*freTrace+i*lhwmax2],lhw[i]*2);
	  }
    }
	else vclrf(&buff[freTrace*saveTrace],freTrace);
    if(il1>=outS[myid] && il1<outE[myid]) 
	{
	  vmovi(ithdr,&ihead[saveTrace*headers],headers);
	  if(FreBa>0 && kk>0)
	  {
	    vclrf(trace1,lenft0+10);
	    vmovf(trace,trace1,samples);
        ForwardPrimeFFT( trace1, lenft0 );
	    for(i=0;i<lhi0;i++) 
	    {
		  j=i+i;
		  trace1[i]=trace1[j]*trace1[j]+trace1[j+1]*trace1[j+1];
		  trace1[i]=sqrt(trace1[i]); 
	    }
	    if(naveFreq==0) { vmovf(trace1,aveFreq,lhi0); naveFreq=1;}
	    else
	    {
		  a2=1./(float)(naveFreq+1);
		  a1=(float)naveFreq*a2;
		  naveFreq++;
	      for(i=0;i<lhi0;i++) 
		    aveFreq[i]=aveFreq[i]*a1+trace1[i]*a2;
	    }
	  }
	}
  }
  saveTrace++;
  if(saveTrace==icrb || cdp1==idmap-1) 
  {
    if(cdp1>=outS[myid]*icrb+icrb-1 && cdp1<outE[myid]*icrb) 
	{
	  put_head(rhead,cdp1-outS[myid]*icrb-saveTrace+1,saveTrace);
      if(AmpBa) put_ampl(ampl,cdp1-outS[myid]*icrb-saveTrace+1,saveTrace);
	}
    put_trace(buff,cdp1-saveS[myid]*icrb-saveTrace+1,saveTrace,0);
    if((fxyF && radonF) || radonF==3)  
      put_trace(buff,cdp1-saveS[myid]*icrb-saveTrace+1,saveTrace,1);
/*	  if(cdp1==saveS[myid]*icrb+icrb*5-1)
	  {
		for(k=500;k<550;k++)
		{
		  fprintf(fp,"%8.2f",buff[k]);
		  if(k%10==9) fprintf(fp,"    cdp1,k %7d %7d \n",cdp1,k);
		}
		fprintf(fp,"    cdp1,k %7d %7d \n",cdp1,k);
	  }*/
/*   fprintf(fp,"after put_trace saveTrace,cdp,cdp1,cdpp,il,xl %7d %7d %7d %7d %7d %7d \n",saveTrace,cdp,cdp1,cdpp,il,xl);
   fflush(fp);*/
	saveTrace=0;
  }
//   fprintf(fp,"done traceDisk \n");
//   fflush(fp);
}

void get_buff(int numG,int xG,int numcG,int xcG,int kk)
{
  int i,j,k,numDI,numDX,overI,overX,over;
  int i1,i2,j1,si1,si2,si3,sx2,sx3;

  if(fxyF && fxyDone==0) 
  {
	overI=overX=overlapF;
	numDI=numDX=numDesiF;
  }
  else
  {
	overI=overlapI;
	overX=overlapX;
	numDI=numDesiI;
	numDX=numDesiX;
	if(radonF==4) 
	{
	  numDX=numDI; 
	  overX=overI;
	}
  }
  if(myid==numprocs-1) over=overlap1;
  else over=overX;

//   fprintf(fp,"numDI,numDX,overI,overX,over %7d %7d %7d %7d %7d \n",numDI,numDX,overI,overX,over);
//   fflush(fp);
  if(outS[myid]==0) 
  {
	saveLine=(numDX-2*overX)*numcG-overlap1;
    if(numcG==0) outLine=overlap1;
	else outLine=overX;
  }
  else 
  {
    saveLine=(numDX-2*overX)*numcG;
	outLine=overX;
  }
  saveLine_old=saveLine-numDX+2*overX;
  if(xcG==0) xsave_old=overlap1-numDI;
  else xsave_old=xsave;

//  if(outE[myid]==nlines || ((radonF==1 || radonF==3) && sign==0 && (fxyDone || !fxyF)))
  if(outE[myid]==nlines)
  {
    if(saveLine>saveE[myid]-saveS[myid]+overlap1-numDX) 
    {
	  saveLine=saveE[myid]-saveS[myid]+overlap1-numDX;
	  outLine=numDX-overX-saveLine+saveLine_old;
    }
  }
  else
  {
    if(saveLine>saveE[myid]-saveS[myid]-numDX) 
    {
	  saveLine=saveE[myid]-saveS[myid]-numDX;
	  outLine=numDX-overX-saveLine+saveLine_old;
    }
  }
//   fprintf(fp,"outS[myid],outLine,saveLine %7d %7d %7d \n",outS[myid],outLine,saveLine);
//   fflush(fp);
  //si1:start line in buff si2+si1: start in disk si3:last line in buff i2:header start
  if(outS[myid]==0 && saveLine<0) { si1=overlap1; si2=-overlap1; }
  else { si1=0; si2=saveLine; }
  if(outS[myid]==0) i2=si2; 
  else i2=si2-overX;
  if(outE[myid]==nlines && saveLine>saveE[myid]-saveS[myid]-numDX) 
    si3=saveE[myid]-saveS[myid]-saveLine;
  else si3=numDX;
  xsave=(numDI-2*overI)*xcG-overlap1;
  if(xcG==0) ouStart=overlap1;
  else ouStart=overI;
  if(xsave>icrb+overlap1-numDI) 
  {
	xsave=icrb+overlap1-numDI;
	ouStart=numDI-overI-xsave+xsave_old;
  }
  if(xsave<0) { sx1=0; sx2=overlap1; sx3=overlap1; }
  else 
  {  
	if(xsave>icrb-numDI) { sx1=xsave; sx2=0; sx3=numDI-icrb+xsave; }
	else { sx1=xsave; sx2=0; sx3=0; }
  }
//   fprintf(fp,"get_buff numcG,saveLine_old,saveLine,outLine,xlS,sx1 %7d %7d %7d %7d %7d %7d \n",numcG,saveLine_old,saveLine,outLine,xlS,sx1);
//   fflush(fp);
//  vclri(ihead,headers*numDesi2);
/*   fprintf(fp,"xcG,xsave_old,xsave,ouStart %7d %7d %7d %7d  \n",xcG,xsave_old,xsave,ouStart);
   fprintf(fp,"si1,si2,si3,sx1,sx2,sx3 %7d %7d %7d %7d %7d %7d  \n",si1,si2,si3,sx1,sx2,sx3);
   fflush(fp);*/
  if(xlS<sx1+numDI)
  {
	j1=ileng;
	if(sx1+ileng>icrb) j1=icrb-sx1;
//    fprintf(fp,"get_trace j1,sx1 %7d %7d  \n",j1,sx1);
//    fflush(fp);
    for(i=si1;i<si3;i++)
    {
//    fprintf(fp,"get_buff saveLine,si1,si2,si3,i,i2 %7d %7d %7d %7d %7d %7d \n",saveLine,si1,si2,si3,i,i2);
//    fflush(fp);
	  if(i>=nlines+overlap1) break;
	  get_trace(&buff[(i*ileng+sx2)*freTrace],(si2+i)*icrb+sx1,j1,kk);
	  if((kk==0 && (radonF==0 || (fxyF==0 && radonF!=3) || (fxyF>0 && radonF==3 && fxyDone))) ||
	  (kk>0 && radonF>0 && ((fxyF>0 && radonF!=3) || (fxyF==0 && radonF==3))))
	  if(xcG==0 && i>=outLine && 
	      ((numcG<numG-1 && i<numDX-overX) || (numcG==numG-1 && i<numDX-over)))
	  {
 //   fprintf(fp,"(1)i,outLine,i2 %7d %7d %7d \n",i,outLine,i2);
 //   fflush(fp);
	    get_head(&rhead[(i-outLine)*icrb*headers],(i2+i)*icrb,icrb);
	  }
	  if(xsave<0) 
	  {
	    for(j=0;j<overlap1;j++)
         vmovf(&buff[(i*ileng+2*overlap1-j)*freTrace],&buff[(i*ileng+j)*freTrace],freTrace);
	  }
	  if(sx1+j1>=icrb) 
	  {
	    if(xlS==0)
	    {
	      for(j=0;j<overlap1;j++)
            vmovf(&buff[(i*ileng+j1+sx2-j-2)*freTrace],
		          &buff[(i*ileng+j1+sx2+j)*freTrace],freTrace);
	    }
	    else
	    {
	      for(j=0;j<overlap1;j++)
		    if(j1+j<ileng)
              vmovf(&buff[(i*ileng+j1-j-2)*freTrace],
		          &buff[(i*ileng+j1+j)*freTrace],freTrace);
	    }
	  }
    }
    xlSo=sx1;
	if(j1<ileng) j1=ileng;
    xlS=sx1+j1;
	if(xlS>icrb+overlap1) xlS=icrb+overlap1; 
 //   fprintf(fp,"xlS,xlSo,j1,sx1 %7d %7d %7d %7d  \n",xlS,xlSo,j1,sx1);
 //   fflush(fp);
    if(saveLine<0 && outS[myid]==0) 
    {
	  i=ileng*freTrace;
	  for(j=0;j<overlap1;j++) vmovf(&buff[(2*overlap1-j)*i],&buff[j*i],i);
    }
    if(saveLine+numDX>saveE[myid] && myid==numprocs-1) 
    {
	  i=ileng*freTrace;
	  for(j=0;j<saveLine+numDX-saveE[myid]+saveS[myid];j++)
        vmovf(&buff[(numDX-2*(saveLine+numDX-saveE[myid]+saveS[myid])+j-1)*i],&buff[(numDX-j-1)*i],i);
    }
  }
  if(xcG==xG-1) saveLine_old=saveLine;
}

void put_buff(int numG,int numcG,int kk)
{
  int i,j,ii,numDX,overX;

  if(fxyF && fxyDone==0) 
  {
	numDX=numDesiF;
	overX=overlapF;
  }
  else
  {
	overX=overlapX;
	numDX=numDesiX;
  }

  if(numG==numcG+1 && myid==numprocs-1) ii=numDX-overlap1;
  else ii=numDX-overX;
//   fprintf(fp,"put_buff saveLine,xsave %7d %7d  \n",saveLine,xsave);
//   fflush(fp);
  for(i=outLine;i<ii;i++)
  {
/*    if(!fxyF || kk==0) j=saveLine+i;*/
    if(kk==0 && radonF!=3) j=saveLine+i;
	else
	{
      if(outS[myid]==0) j=saveLine+i;     
	  else j=saveLine+i-overX;
	}
	put_trace(&buff1[(i-outLine)*icrb*freTrace],j*icrb,icrb,kk);
  }
}

void ouPuts(int numG,int numcG)
{
  int i,i1,j,j1,k,m,n,n1,k1,ii,numDX,overX,cdp0;
  float a1,a2,a3;
  static int kk=0; 

  if(fxyF && fxyDone==0) 
  {
	numDX=numDesiF;
	overX=overlapF;
  }
  else
  {
	overX=overlapX;
	numDX=numDesiX;
  }

  if(numcG==0) kk=0;

 //  fprintf(fp,"ouPuts numG,numcG,outLine,saveLine %5d %5d %5d %5d \n",numG,numcG,outLine,saveLine );
 //  fflush(fp);
  a1=1./(float)sampleTap;
/*   fprintf(fp,"enter ouRadon saveLine,outLine,countGroup %7d %7d %7d  \n",saveLine,outLine,countGroup);
   fprintf(fp,"fxyDone,radonDone,velo %7d %7d %7d  \n",fxyDone,radonDone,velo);
   fflush(fp);*/
/*	  if(numcG==0)
	  {
        for(j=100;j<200;j++)
		{
		  fprintf(fp,"%9.1f",buff1[10*freTrace+j]);
		  if(j%10==9) fprintf(fp,"j %8d \n",j);
		}
        fprintf(fp,"buff1 j %8d \n",j);
	  }*/
  if(numG==numcG+1 && outE[myid]==nlines) ii=numDX-overlap1;
  else ii=numDX-overX;
//   fprintf(fp,"numG,numcG,outE[myid],ii %5d %5d %5d %5d \n",numG,numcG,outE[myid],ii);
//   fflush(fp);
  for(i=outLine;i<ii;i++)
  {
	i1=(i-outLine)*icrb;
	n1=m=0;
/*   if(numcG==2)
   {
	  fprintf(fp,"i,ii,i1,n1 %5d %5d %5d %5d \n",i,ii,i1,n1 );
   fflush(fp);
   }*/
	for(j=0;j<icrb;j++)
    {
	  vmovi(&ihead[(i1+j)*headers],ithdr2,headers);
	  cdp=ithdr2[indexCdp];
      il=(cdp-minCDP)/ncrbi+minIlin-lmin;
	  xl=j;
      cdp=il*icrb+j;
	  if(j==0) cdp0=cdp;
//	  if(numcG==2 && i==17) {fprintf(fp,"j,m,n1,cdp,il,indexcdp %5d %5d %5d %5d %5d %5d \n",j,m,n1,cdp,il,ithdr2[indexCdp]); fflush(fp);}
	  if(j==0 && AmpBa) get_ampl(ampl,cdp-outS[myid]*icrb,icrb);
	  if(ithdr2[indexTrc_type]!=1 || ithdr2[indexCdp]==0 ) continue;
	  j1=(i1+j)*freTrace;
	  for(k=0;k<numWin;k++)
	  {
	      vclrf(trace1,lenft+10);
		  k1=k*(sampleWin-sampleTap);
		  vmovf(&buff1[j1+k*lhwmax2],&trace1[low[k]*2],lhw[k]*2);
          ReversePrimeFFT( trace1, lenft );
	      if(k==0) vmovf(trace1,trace2,sampleWin);
	      if(numWin>1 && k>0) 
	      {
		    for(n=0;n<sampleTap;n++) 
		      trace2[k1+n]=(trace2[k1+n]*(sampleTap-n)+trace1[n]*n)*a1;
            vmovf(&trace1[sampleTap],&trace2[k1+sampleTap],sampleWin-sampleTap);
		  }
	  }
/*	  if(i==outLine && j==10) 
	  {
		  fprintf(fp,"mute[cdp*2],mute[cdp*2+1] %6d %6d \n",mute[cdp*2],mute[cdp*2+1]);
		for(k=500;k<550;k++)
		{
		  fprintf(fp,"%8.2f",trace2[k]);
		  if(k%10==9) fprintf(fp,"    i,j,k %7d %7d %7d \n",i,j,k);
		}
		fprintf(fp,"absR[cdp],i,j,k %10.1f %7d %7d %7d \n",absR[cdp],i,j,k);
	  }*/
//	  if(muteop)
	  {
	      if(mute[cdp*2]>0) vclrf(trace2,mute[cdp*2]);
		  if(mute[cdp*2+1]<samples-1) 
			  vclrf(&trace2[mute[cdp*2+1]+1],samples-mute[cdp*2+1]-1);
	  }
      vabsf(trace2,mute[cdp*2],mute[cdp*2+1]);
/*	  if(i==outLine && j==10) 
	  {
		  fprintf(fp,"mute[cdp*2],mute[cdp*2+1] %6d %6d \n",mute[cdp*2],mute[cdp*2+1]);
		for(k=500;k<550;k++)
		{
		  fprintf(fp,"%8.2f",trace2[k]);
		  if(k%10==9) fprintf(fp,"    i,j,k %7d %7d %7d \n",i,j,k);
		}
		fprintf(fp,"    i,j,k %7d %7d %7d \n",i,j,k);
	  }*/
	  buff[m]=(float)cdp;
//	  if(numcG==2 && i==17) fprintf(fp,"ithdr2[indexCdp] %7d  \n",ithdr2[indexCdp]);
	  vmovf(rthdr2,&buff[m+1],headers);
	  vmovf(trace2,&buff[m+headers+1],samples);
	  if(AmpBa)
	  {
        vmovf(trace2,trace,samples);
//	  if(numcG==4 && i==17) {fprintf(fp,"before getampl j %5d \n",j);    fflush(fp);}
		getAmpl();
//	  if(numcG==4 && i==17) {fprintf(fp,"before fixampl \n");    fflush(fp);}
		fixAmpl(j);
	  }
//	  if(numcG==4 && i==17) {fprintf(fp,"after fixampl \n");    fflush(fp);}
	  if(FreBa>0)
	  {
	    vclrf(trace1,lenft0+10);
	    vmovf(&buff[m+headers+1],trace1,samples);
        ForwardPrimeFFT( trace1, lenft0 );
	    for(k=0;k<lhi0;k++) 
	    {
		  k1=k+k;
		  trace1[k]=trace1[k1]*trace1[k1]+trace1[k1+1]*trace1[k1+1];
		  trace1[k]=sqrt(trace1[k]); 
	    }
	    if(naveFreq==0) { vmovf(trace1,aveFout,lhi0); naveFreq=1;}
	    else
	    {
		  a2=1./(float)(naveFreq+1);
		  a3=(float)naveFreq*a2;
		  naveFreq++;
	      for(k=0;k<lhi0;k++) 
		    aveFout[k]=aveFout[k]*a3+trace1[k]*a2;
	    }
	  }
	  m+=shtotal;
	  n1++;
//	  if(numcG==2 && i==17) {fprintf(fp,"j,m,n1 %5d %5d %5d \n",j,m,n1);    fflush(fp);}
	}
	if(n1>0)
	{
/*   if(numcG==2)
   {
	  fprintf(fp,"cdp0 %5d \n",cdp0 );
   fflush(fp);
   }*/
	  if(AmpBa) put_ampl(ampl,cdp0-outS[myid]*icrb,icrb);
/*   if(numcG==2)
   {
	  fprintf(fp,"cdp0(1) %5d \n",cdp0 );
   fflush(fp);
   }*/
	  if(nbin==1)	
	  {
		offDisk[il-outS[myid]]=n1;
		put_bin(buff,kk,n1,0);
/*		if(i==outLine+5)
		{
	    for(k=400;k<500;k++)
	    {
	      fprintf(fp,"%9.1f",buff[shtotal*10+1+headers+k]);
	      if(k%10==9) fprintf(fp," %8d \n",k);
	    }
	    fprintf(fp,"%6d %6d \n",k,i);
		}*/
	  }
	  else 
	  {
		offDisk[(outE[myid]-outS[myid])*(bin-1)+il-outS[myid]]=n1;
		put_bin(buff,kk,n1,bin-1);
 //  fprintf(fp,"il,n1,bin,kk %5d %5d %5d %5d \n",il,n1,bin,kk);
 //  fflush(fp);
	  }
	  kk+=n1;
	}
  }
//   fprintf(fp,"ouPuts(1) numG,numcG,outLine,saveLine %5d %5d %5d %5d \n",numG,numcG,outLine,saveLine );
//   fflush(fp);
}
  
void fixAmpl(int xl)
{
  int i,j,k,n,nOper,ami,ama,nOp2;
  float amin,amax,a1;

  nOper=(int)(AmpLong/sampleRate);
  nOp2=nOper/2;
  k=sampleAs*xl;
  amax=0.;
  amin=1000000000.;
/*  if(il==18 && xl==38)
  {
	    fprintf(fp,"trace2 %6d %6d %6d \n",cdp,il,xl);
	    for(i=0;i<sampleAs;i++)
	    {
	      fprintf(fp,"%9.1f",trace2[i]);
	      if(i%10==9) fprintf(fp," %8d \n",i);
	    }
	    fprintf(fp,"%6d %6d \n",i,il);

	    fprintf(fp,"ampl %6d %6d %6d \n",cdp,il,xl);
	    for(i=0;i<sampleAs;i++)
	    {
	      fprintf(fp,"%9.1f",ampl[k+i]);
	      if(i%10==9) fprintf(fp," %8d \n",i);
	    }
	    fprintf(fp,"%6d %6d \n",i,il);
		fflush(fp);
  }*/

  for(i=nOp2/gap;i<=sampleAs-nOp2/gap;i++)
  {
    if(trace2[i]<=0. || ampl[k+i]<=0.) 
	{
	  trace[i]=0.;
	  trace2[i]=0.;
	  continue;
	}
    trace[i]=ampl[k+i]/trace2[i];       // int/out
  }
  for(i=nOper/gap;i<=sampleAs-nOper/gap;i++)
  {
    if(trace2[i]<=0.) continue; 
	if(amin>trace[i]) { amin=trace[i]; ami=i; }
	if(amax<trace[i]) { amax=trace[i];  ama=i; }
  }
  if(amax==0.) { vclrf(trace,sampleAs); amin=0.;}
  for(i=0;i<nOp2/gap;i++)
  {
	trace[i]=trace[nOp2/gap];
	trace[sampleAs-i]=trace[sampleAs-nOp2/gap];
  }
  trace[0]=amin;
  trace[1]=amax;
  vmovf(trace,&ampl[k],sampleAs);
  if(amin>0.) { tolAmin+=amin; tolAmax+=amax; ntolA++; }

if(amin>0.) a1=amax/amin;
else a1=0.;
//	if(a1>AmpParm) { a1=AmpParm; amin=amax/AmpParm; }
//	    if(il>45 && il<=75 && xl>45 && xl<=75) 
//fprintf(fp,"il,xl,amin,ami,amax,ama,amax/amin %5d %5d %7.4f %5d %7.2f %5d %9.3f \n",il,xl,amin,ami,amax,ama,a1);
//fflush(fp);
/*		if(xl==10)
		{
	    fprintf(fp,"ampl %6d %6d %6d \n",cdp,il,xl);
	    for(n=0;n<sampleAs;n++)
	    {
	      fprintf(fp,"%9.1f",ampl[k+n]);
	      if(n%10==9) fprintf(fp," %8d \n",n);
	    }
	    fprintf(fp,"%6d %6d \n",n,il);
		}*/
/*		if(xl==10)
		{
	    fprintf(fp,"trace %6d %6d %6d \n",cdp,il,xl);
	    for(n=0;n<sampleAs;n++)
	    {
	      fprintf(fp,"%9.3f",trace[n]);
	      if(n%10==9) fprintf(fp," %8d \n",n);
	    }
	    fprintf(fp,"%6d %6d \n",n,il);
		}*/
/*  trace2[0]=trace[0];
  trace2[1]=trace[1];
  for(i=2;i<=sampleAs-3;i++)
	trace2[i]=(trace[i-2]+trace[i-1]+trace[i]+trace[i+1]+trace[i+2])*0.2;
  trace2[sampleAs-2]=trace[sampleAs-2];
  trace2[sampleAs-1]=trace[sampleAs-1];*/

//  vmovf(trace,trace2,sampleAs);      // trace==trace2==in/out
/*		if(xl==10)
		{
	    fprintf(fp,"trace2 %6d %6d %6d \n",cdp,il,xl);
	    for(n=0;n<sampleAs;n++)
	    {
	      fprintf(fp,"%9.3f",trace2[n]);
	      if(n%10==9) fprintf(fp," %8d \n",n);
	    }
	    fprintf(fp,"%6d %6d \n",n,il);
		}*/
/*	if(a1>=AmpParm) a1=AmpParm-1.;
	else a1-=1.;
  for(i=0;i<sampleAs;i++) trace2[i]=amin*(1.+(AmpParm-1.)*(trace2[i]-amin)/(amax-amin));*/
/*		if(xl==10)
		{
	    fprintf(fp,"trace2 %6d %6d %6d \n",cdp,il,xl);
	    for(n=0;n<sampleAs;n++)
	    {
	      fprintf(fp,"%9.3f",trace2[n]);
	      if(n%10==9) fprintf(fp," %8d \n",n);
	    }
	    fprintf(fp,"%6d %6d \n",n,il);
		}*/
}

void Perform()
{
  int i,j,k,numG,xG;
  int overX;

/*        for(j=saveS[myid]*icrb;j<saveE[myid]*icrb;j++)
		{
		  fprintf(fp,"%6d",mute[j*2]);
		  if(j%10==9) fprintf(fp,"j %8d \n",j);
		}
        fprintf(fp,"mute j %8d \n",j);*/

  if(FreBa>0) naveFreq=0;
  if(AmpBa>0) { tolAmin=0.; ntolA=0; tolAmax=0.;}
  if(fxyF>0) 
  {
	numG=numGroup;
	xG=xGroup;
  }
  else 
  {
	numG=numGroupR;
	xG=xGroupR;
  }

//  fprintf(fp,"perform numG,xG %8d %6d \n",numG,xG);
//   fflush(fp);
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0)
  {
  if(fxyF>0)
  {
    PrintTime();
    printf("start fxy filtering for bin = %5d \n",bin);
  }
  else
  {
    if(radonF==1 || radonF==3) 
	{
	  PrintTime();
      printf("start Radon inline direction filtering for bin = %5d \n",bin);
	}
    if(radonF==2) 
	{
	  PrintTime();
      printf("start Radon xline direction filtering for bin = %5d \n",bin);
	}
    if(radonF==4) 
	{
	  PrintTime();
      printf("start Radon one pass filtering for bin = %5d \n",bin);
	}
  }
  }
 //  fprintf(fp,"Perform numG,xG %7d %7d  \n",numG,xG );
 //  fflush(fp);
  if(myid>0) 
  {
  for(i=0;i<numG;i++)
  {
    xlS=0;
	for(j=0;j<xG;j++)
	{
//   fprintf(fp,"i,j %7d %7d  \n",i,j );
//   fflush(fp);
	  get_buff(numG,xG,i,j,0);
/*	  if(j==0)
	  {
        for(k=500;k<702;k++)
		{
		  fprintf(fp,"%10.1f",buff[k]);
		  if(k%10==9) fprintf(fp,"   k %8d \n",k);
		}
        fprintf(fp,"i,j,k %8d %8d %8d \n",i,j,k);
	  }*/
	  if(fxyF>0) FxyPerform(numG,xG,i,j);
	  else
	  {
	    if(radonF==1 || radonF==3)
	    {
	      sign=0;
	      lineRadon(i,j);
	    }
	    if(radonF==2)
	    {
	      sign=1;
	      lineRadon(i,j);
	    }
	    if(radonF==4) Radon3d(i,j);
	  }
	  if(j==xG-1)
	  {
	    if(radonF==0 || (fxyF==0 && radonF!=3)) ouPuts(numG,i);
	    else put_buff(numG,i,1);
	  }
	}
  }
  }

//   fprintf(fp,"done perform \n");
//   fflush(fp);
  if(!radonF) return;
  if(!fxyF && radonF!=3) return;
  if(fxyF>0) fxyDone=1;
  if(radonF==3 && fxyF==0) radonDone=1;

  MPI_Barrier(MPI_COMM_WORLD);
  if(fxyF>0) overX=overlapF;
  else overX=overlapX;
	for(k=0;k<numprocs;k++)
	{
	  if(saveS[k]>0) saveS[k]+=overX;
	  if(saveE[k]<nlines) saveE[k]-=overX;
	  if(saveE[k]<0) saveE[k]=0;
//	fprintf(fp,"(1)k,saveS[k],saveE[k],outS[k],outE[k] %5d %5d %5d %5d %5d \n",k,saveS[k],saveE[k],outS[k],outE[k]); 
//  fflush(fp);
	}
    i=saveE[myid]-saveS[myid]-overlapX*2;
	if(saveS[myid]==0 || myid==numprocs-1) i+=overlap1;
	if(saveS[myid]==0 && myid==numprocs-1) i+=overlap1;
	numGroupR=i/(numDesiX-2*overlapX);
	if(i%(numDesiX-2*overlapX)>0) numGroupR++;
//	fprintf(fp,"overlapX,numGroupR,xGroupR %5d %5d %5d \n",overlapX,numGroupR,xGroupR); 
//  fflush(fp);
/*  for(i=0;i<numprocs;i++)
  {
	fprintf(fp,"(1)i,saveS[i],saveE[i],outS[i],outE[i] %5d %5d %5d %5d %5d \n",i,saveS[i],saveE[i],outS[i],outE[i]); 
  fflush(fp);
  }*/
  if(myid==0)
  {
	if(radonF==1 || (fxyF==1 && radonF==3) )
	{
	  PrintTime();
      printf("start Radon inline direction filtering for bin = %5d \n",bin);
	}
	if(radonF==2 || (fxyF==0 && radonF==3))
	{
	  PrintTime();
      printf("start Radon xline direction filtering for bin = %5d \n",bin);
	}
    if(radonF==4) 
	{
	  PrintTime();
      printf("start Radon one pass filtering for bin = %5d \n",bin);
	}
  }

  if(myid>0)
  {
  numG=numGroupR;
  xG=xGroupR;
 // fprintf(fp,"numG,xG %8d %6d \n",numG,xG);
  for(i=0;i<numG;i++)
  {
  xlS=0;
	for(j=0;j<xG;j++)
	{
//  fprintf(fp,"i,j %8d %6d \n",i,j);
	  get_buff(numG,xG,i,j,1);
	  if(radonF==1 || (fxyF==1 && radonF==3) )
	  {
	    sign=0;
	    lineRadon(i,j);
	  }
	  if(radonF==2 || (fxyF==0 && radonF==3))
	  {
	    sign=1;
	    lineRadon(i,j);
	  }
	  if(radonF==4) Radon3d(i,j);
	  if(j==xG-1)
	  if(fxyF!=1 || radonF!=3) ouPuts(numG,i);
	  else put_buff(numG,i,0);
	}
  }
  }

  if(fxyF==0 || radonF!=3) return; 
  MPI_Barrier(MPI_COMM_WORLD);
  if(myid==0)
  {
	  PrintTime();
      printf("start Radon xline direction filtering for bin = %5d \n",bin);
  }
  if(myid==0) return;
	for(k=0;k<numprocs;k++)
	{
	  if(saveS[k]>0) saveS[k]+=overlapX;
	  if(saveE[k]<nlines) saveE[k]-=overlapX;
	  if(saveE[k]<0) saveE[k]=0;
	}
    i=saveE[myid]-saveS[myid]-overlapX*2;
	if(saveS[myid]==0 || myid==numprocs-1) i+=overlap1;
	if(saveS[myid]==0 && myid==numprocs-1) i+=overlap1;
	numGroupR=i/(numDesiX-2*overlapX);
	if(i%(numDesiX-2*overlapX)>0) numGroupR++;
  numG=numGroupR;
  xG=xGroupR;
  for(i=0;i<numG;i++)
  {
  xlS=0;
	for(j=0;j<xG;j++)
	{
	  get_buff(numG,xG,i,j,0);
	  sign=1;
	  lineRadon(i,j);
	  if(j==xG-1)
	  ouPuts(numG,i);
	}
  }
}

void FxyPerform(int numG,int xG,int k,int n)
{
  int i,j;

	  for(j=0;j<numWin;j++)
	  {
	    for(fre=0;fre<lhw[j];fre++)
	    {
		  Xij(j,n);
		  Rij();
		  Rsolve(j);
		  if(G[0]==100000000.) 
			ProgramErr("Problem in solving equation" );
	    }
	  }
	  ouTrace(numG,xG,k,n);
 //  fprintf(fp,"finish ouTrace \n" );
 //  fflush(fp);
}

void Xij(int jj,int xcG)
{
  /* jj(0->numWin-1) nstart (start point) */
  int i,j,k1,k2,num1,num;

  num1=(jj*lhwmax+fre)*2;
  if(xcG==0 || xlSo>0) num=overlap1;
  else num=0;

  for(i=0;i<numDesiF;i++)
  {
	k1=(i*ileng-xlSo+sx1-num+overlap1)*freTrace+num1;
	k2=i*numDesiF*2;
	for(j=0;j<numDesiF;j++) vmovf(&buff[k1+j*freTrace], &proc[0].xij[k2+2*j],2);
  }
}

void Rij()
{
  int i1,i2,i22,j1,j2,num;
  int numri1,numxi1,numri2,numxi2,numrj1,numxj1,numrj2,numxj2;

/*   fprintf(fp,"Rij fre %7d \n",fre);
   fflush(fp);*/
  num=numFilt/2;
  vclrf(rrij,ac*half/2);
  for(i1=0;i1<=num;i1++)
  {
	if(i1<num) i22=numFilt;
	else i22=numFilt/2;
	numri1=i1*half*numFilt;
	numxi1=i1*numDesiF;
	for(i2=0;i2<i22;i2++)
	{
	  numri2=numri1+i2*half;
	  numxi2=numxi1+i2;
	  for(j1=0;j1<numFilt;j1++)
	  {
	    numrj1=numri2+j1*numFilt;
	    numxj1=j1*numDesiF;
		for(j2=0;j2<numFilt;j2++)
		{
	      numrj2=numrj1+j2;
	      numxj2=numxj1+j2;
		  Get_rij(i2,j1,j2,numrj2,numxi2,numxj2,i1);
		} /*j2*/
	  } /*j1*/
	} /*i2*/
  } /*i1*/

/*  vmovf(&rrij[ac*(half+1)/2],&rrij[ac*(half-1)/2],ac*(half-1)/2);*/
}

void Get_rij(int i2,int j1, int j2,int numrj2,int numxi2,int numxj2,int ch11)
{
	int k;

		  if(ch11==j1 && i2==j2)
		  {
			caddr(&proc[0].xij[numxj2*2],&proc[0].xij[numxj2*2],
			      &rrij[numrj2*2],numDesiF-numFilt+1);
		  }
		  if(ch11<j1)
		  {
			  caddc(&proc[0].xij[numxj2*2],&proc[0].xij[numxi2*2],&rrij[numrj2*2],numDesiF-numFilt+1);
		  }
		  if(ch11==j1 && i2<j2)
		  {
			  caddc(&proc[0].xij[numxj2*2],&proc[0].xij[numxi2*2],&rrij[numrj2*2],numDesiF-numFilt+1);
		  }
		  if(ch11>j1)
			{
			  k=(j1*half*numFilt+j2*half+ch11*numFilt+i2)*2;
			  rrij[numrj2*2]=rrij[k];
			  rrij[numrj2*2+1]=-rrij[k+1];
			}
		  if(ch11==j1 && i2>j2)
			{
			  k=(j1*half*numFilt+j2*half+ch11*numFilt+i2)*2;
			  rrij[numrj2*2]=rrij[k];
			  rrij[numrj2*2+1]=-rrij[k+1];
			}
/*   fprintf(fp,"numrj2,numxi2,numxj2,rrij[numrj2*2],rrij[numrj2*2+1] %6d %6d %6d %10.1f %10.1f \n",numrj2,numxi2,numxj2,rrij[numrj2*2],rrij[numrj2*2+1]);
   fflush(fp);*/
}
		  
void Rsolve(int jj)
{
  int i,i1,i2,i3,j,j1,result,k1,half1,ac1;
  float rij,rik,iij,iik,norm;
  float saveLimit;
  
  half1=(half-1)/2;     /*4*/
  ac1=half1*2;         /*8*/
  saveLimit=1.;
  k1=(jj*lhwmax+fre)*ac1;
/*   fprintf(fp,"Rsolve jj,fre,k1 %7d %7d %7d \n",jj,fre,k1);
   fflush(fp);*/
  for(i=0;i<half1;i++)
  {
	G[k1+i]=rrij[i*ac+half-1];
	G[k1+i+half1]=rrij[i*ac+half];
  }
  norm=0.;
  for(i=k1;i<k1+ac1;i++) norm+=G[i];
  if(norm==0.)
  {
    vclrf(&G[k1],ac1); 
	return;
  }

  for(i=0;i<half1-1;i++)
  {
	vmovf(&rrij[i*ac+half+1],&rrij[i*(ac-2)+half-1],half-1);
	vmovf(&rrij[(i+1)*ac],&rrij[(i+1)*(ac-2)],half-1);
  }
  vmovf(&rrij[(half1-1)*ac+half+1],&rrij[(half1-1)*(ac-2)+half-1],half-1);

  for(i=0;i<half1;i++)
  {
	i1=i*(ac-2);
	i2=i*ac1;
	i3=(i+half1)*ac1;
	for(j=0;j<half1;j++)
	{
	  j1=2*j;
	  rij=rrij[i1+j1];
	  rik=rrij[i1+ac-4-j1];
	  iij=rrij[i1+j1+1];
	  iik=rrij[i1+ac-3-j1];
	  RR[i2+j]=rij+rik;
	  RR[i3+j+half1]=rij-rik;
	  RR[i2+j+half1]=iik-iij;
	  RR[i3+j]=iik+iij;
	}
  }

  if(RR[0]==0. || saveLimit==0.)
  {
    vclrf(&G[k1],ac1); 
	return;
  }

  norm=RR[0]*0.005;
  for(i=0;i<ac1;i++)
  {
    RR[i*ac1+i]+=norm;
  }

  if(fre==0) saveLimit=RR[0];

  if(RR[0]/saveLimit>0.0000001) 
  {
    result=ludcmp(RR,ac1,idx,vv);
    if(result)
    {
	  failSolve++;
      vclrf(&G[k1],ac1); 
	  return;
    }
    lubksb(RR,ac1,idx,&G[k1]);

    norm=0.;
    for(i=0;i<half1;i++) norm=norm+sqrt(G[k1+i]*G[k1+i]+G[k1+i+half1]*G[k1+i+half1]);
    norm=1.0/norm;
    for(i=0;i<ac1;i++) trace2[i]=G[k1+i]*norm;

    for(i=0;i<half1;i++)
    {
	  G[k1+2*i]=trace2[i];
	  G[k1+2*i+1]=trace2[i+half1];
    }
  }
  else vclrf(&G[k1],ac1); 
}

void Freq(int ii,int xx,int ww,int ff,float *a)
{
  /* ii inline, xx=xline, ww=windos, ff=freq */
  int i,i1,i2,j,j1,k1,numi,numx,num;
  float b[100],c[2]; 

  numx=numFilt/2;
  numi=numFilt/2;
  i1=ileng*(ii-numx)+xx-numi;
  num=(ww*lhwmax+ff)*2;
  k1=(ww*lhwmax+ff)*(half-1);
  for(i=0;i<=numx;i++)
  {
	if(i==numx) j1=numi;
	else j1=numFilt;
	for(j=0;j<j1;j++)
	{
	  i2=(i1+i*ileng+j)*freTrace;
	  vmovf(&buff[i2+num],&b[(i*numFilt+j)*2],2);
	}
  }
  caddi(b,&G[k1],a,half/2);
  c[0]=a[0];
  c[1]=a[1];
  for(i=0;i<=numx;i++)
  {
	if(i==numx) j1=numi;
	else j1=numFilt;
	for(j=0;j<j1;j++)
	{
	  i2=(i1+(numFilt-i-1)*ileng+numFilt-1-j)*freTrace;
	  vmovf(&buff[i2+num],&b[(i*numFilt+j)*2],2);
	}
  }
  caddir(b,&G[k1],a,half/2);
  a[0]+=c[0];
  a[1]+=c[1];
}

void ouTrace(int numG,int xG,int numcG,int xcG)
{
  int i,j,k,m,ii,jj,i1,j1,num,num1;

  if(xcG==0 || xlSo>0) num=overlap1;
  else num=0;
  if(xcG==0 || xlSo==0) num1=overlap1;
  else num1=0;
  if(numG==numcG+1) ii=numDesiF-overlap1;
  else ii=numDesiF-overlapF;
  if(xG==xcG+1) jj=numDesiF-overlap1;
  else jj=numDesiF-overlapF;
  for(i=outLine;i<ii;i++)
  {
	for(j=ouStart;j<jj;j++)
    {
	  for(k=0;k<numWin;k++)
		for(m=0;m<lhw[k];m++) 
		  Freq(i,j-xlSo+sx1+overlap1-num,k,m,&buff1[((i-outLine)*icrb+sx1+j
		         -num+overlap1-num1)*freTrace+2*(k*lhwmax+m)]);
	}
  }
}

void XijLine(procStruct  *ps)
{
	/* jj(0->numWin-1) kk:line nstart (start point) i(line number) */
  int ii,i,j,num1,i2,i1,i3,j1,num;
  float *AA;

  if(ps->xc==0 || xlSo>0) num=overlap1;
  else num=0;
/*        for(j=0;j<32;j++)
		{
		  fprintf(fp1,"%10.1f",ps->xij[j]);
		  if(j%10==9) fprintf(fp1,"   j %8d \n",j);
		}
        fprintf(fp1,"ps->il, j  %8d %8d \n",ps->il,j);
   fflush(fp1);*/
  i2=ps->free+low[ps->win]-low1;
  num1=ps->win*lhwmax2+ps->free*2;

  vclrf(ps->xij,numDesi*2);
  if(method==1)
  {
    AA=&LLsave[i2*numDesiR2];
	if(sign==0)
	{
	  ii=ps->il*ileng;
	  for(i=0;i<numDesiI;i++)
	  {
        i1=i+i;
	    i3=i*numDesiI;
	    for(j=0;j<numDesiI;j++)
	    {
		  j1=(ii+j-xlSo+sx1-num+overlap1)*freTrace+num1;
	      ps->xij[i1]+=buff[j1]*AA[i3+j];
	      ps->xij[i1+1]+=buff[j1+1]*AA[i3+j];
	    }
	  }
    }
    else
    {
	  ii=ps->il-xlSo+sx1-num+overlap1;
	  for(i=0;i<numDesiI;i++)
	  {
        i1=i+i;
	    i3=i*numDesiI;
	    for(j=0;j<numDesiI;j++)
	    {
		  j1=(j*ileng+ii)*freTrace+num1;
	      ps->xij[i1]+=buff[j1]*AA[i3+j];
	      ps->xij[i1+1]+=buff[j1+1]*AA[i3+j];
	    }
	  }
    }
  }
  else
  {
    if(sign==0)
    {
	  ii=ps->il*ileng-xlSo+sx1-num+overlap1;
	  for(j=0;j<numDesiI;j++)
        vmovf(&buff[(ii+j)*freTrace+num1], &ps->xij[j+j],2);
    }
    else
    {
	  for(j=0;j<numDesiX;j++)
        vmovf(&buff[(j*ileng+ps->il-xlSo+sx1-num+overlap1)*freTrace+num1],&ps->xij[j+j],2);
/*  fprintf(fp,"(1)enter XijLine free,sign,ps->il %7d %7d %7d \n",ps->free,sign,ps->il);
   fflush(fp);*/
/*	if(ps->il==10)
	{
        for(j=0;j<32;j++)
		{
		  fprintf(fp1,"%10.1f",ps->xij[j]);
		  if(j%10==9) fprintf(fp1,"   j %8d \n",j);
		}
        fprintf(fp1,"ps->il, j  %8d %8d \n",ps->il,j);
   fflush(fp1);
    }*/
    }
  }
}

void reXijLine(procStruct  *ps)
{
	/* jj(0->numWin-1) kk:line nstart (start point) i(line number) */
  int j,num1,ii,k1,k2,jj,num,num2;

  num1=ps->win*lhwmax2+ps->free*2;
  if(ps->xc==0 || xlSo>0) num=overlap1;
  else num=0;
  if(ps->xc==0 || xlSo==0) num2=overlap1;
  else num2=0;

 if(sign==0)
  {
    if(ps->xc==xGroupR-1) jj=numDesiI-overlap1;
    else jj=numDesiI-overlapI;
    ii=icrb*(ps->il-outLine);
	for(j=ouStart;j<jj;j++)
      vmovf(&ps->xij[j*2],&buff1[(ii+sx1+j-num+overlap1-num2)*freTrace+num1],2);
  }
  else
  {
    if(ps->nc==numGroupR-1 && outE[myid]==nlines) jj=numDesiX-overlap1;
    else jj=numDesiX-overlapX;
	for(j=outLine;j<jj;j++)
      vmovf(&ps->xij[j*2],&buff1[((j-outLine)*icrb+ps->il+sx1-num+overlap1-num2)*freTrace+num1],2);
  }
//   fprintf(fp,"reXij outLine,jj,sx1,num2,ps->il,num %7d %7d %7d %7d %7d %7d \n",outLine,jj,sx1,num2,ps->il,num);
//   fprintf(fp,"reXij ps->il,ps->nc,ps->xc %7d %7d %7d \n",ps->il,ps->nc,ps->xc);
//   fflush(fp);
}

void lineRadon(int numcG,int xcG)
{
  int i,i1,i2,j,k,ii,jj;

/* pthread used  */
  pthread_attr_t  joinAttr;				/* default is wrong so we must set */
  procStruct  *ps;

  pthread_attr_init( &joinAttr );
  pthread_attr_setscope(&joinAttr, PTHREAD_SCOPE_SYSTEM); 
  pthread_attr_setdetachstate( &joinAttr, PTHREAD_CREATE_DETACHED );

   /*     for(k=500;k<702;k++)
		{
		  fprintf(fp,"%10.1f",buff[k]);
		  if(k%10==9) fprintf(fp,"   k %8d \n",k);
		}
        fprintf(fp,"i,j,k %8d %8d %8d \n",i,j,k);*/
  if(sign==0)
  {
	i1=outLine;
    if(numGroupR==numcG+1) i2=numDesiX-overlap1;
    else i2=numDesiX-overlapX;
  }
  else
  {
	i1=ouStart;
    if(xGroupR==xcG+1) i2=numDesiI-overlap1;
    else i2=numDesiI-overlapI;
  }

/*   fprintf(fp,"lineRadon i1,i2,outLine,numcG,xcG,sign,numWin %7d %7d %7d %7d %7d %7d %7d \n",i1,i2,outLine,numcG,xcG,sign,numWin);
   fflush(fp);*/
  k=0;
  for(i=i1;i<i2;i++)
  {
//    pthread_join( workThreads[k], NULL );	 
	proc[k].il=i;
	proc[k].nc=numcG;
	proc[k].xc=xcG;
	proc[k].threadID=k;
	procThread();
/*  fprintf(fp,"first i,il,k %7d %7d %7d  \n",i,proc[k].il,k);
   fflush(fp);*/
//	k++;
//	if(k==numThreads) k=0;
  }
//  for(i=0;i<numThreads;i++) pthread_join( workThreads[i], NULL );
//  for(i=0;i<numThreads;i++) pthread_join( workThreads[i], NULL );
}

void strongD(procStruct  *ps,int np,int ni)
{
  int k,i2,nd,ndip,*diplim;
  float a1,a2,stDip;

  if(sign==0) 
  {
	ndip=ndipI;
	stDip=stDipI[ps->win];
	diplim=diplimI;
  }
  else 
  {
	ndip=ndipX;
	stDip=stDipX[ps->win];
	diplim=diplimX;
  }

  nd=ni%ndip;

  //fprintf(fp,"strongD stDip,ni,np,sign,ndip,diplim[0],ps->win %7.2f %7d %7d %7d %7d %7d %7d \n",stDip,ni,np,sign,ndip,diplim[0],ps->win);
  // fflush(fp);
  i2=low[ps->win]-low1;
  a1=2.+stDip;
  vclrf(ps->trace2,lenft+10);
  for(k=0;k<lhw[ps->win];k++)
	if(nd<diplim[k+i2]) vmovf(&ps->radon[k*np+ni],&ps->trace2[(low[ps->win]+k)*2],2);
  ReversePrimeFFT( ps->trace2, lenft );
  for(k=0;k<lenft;k++)
  {
    a2=fabs(ps->trace2[k]);
    if(a2>0.0000001)
	  ps->trace2[k]=pow(a2,a1)/ps->trace2[k];
  }
  ForwardPrimeFFT( ps->trace2, lenft );
  for(k=0;k<lhw[ps->win];k++)
  {
	if(nd>=diplim[k+i2]) vclrf(&ps->radon[k*np+ni],2);
	else vmovf(&ps->trace2[(low[ps->win]+k)*2],&ps->radon[k*np+ni],2);
  }
}

void  procThread()		/* routine threads run */
{
  int i,j,k,np,ndip;
  float stDip;
  procStruct  *ps;

  ps =proc;

//   fprintf(fp,"enter procThread ps->threadID,ps->il,sign,,stDipX %5d %5d %5d %5d %7.3f \n",ps->threadID,ps->il,sign,ndipX,stDipX );
//   fflush(fp);
  for(ps->win=0;ps->win<numWin;ps->win++)
  {
  if(sign==0) 
  {
	ndip=ndipI;
	stDip=stDipI[ps->win];
  }
  else 
  {
	ndip=ndipX;
	stDip=stDipX[ps->win];
  }
  np=ndip*4;
/*  fprintf(fp,"(1)ps->win,resoluI[ps->win],sign %7d %7.3f %5d \n",ps->win,resoluI[ps->win],sign);
   fflush(fp);
   fprintf(fp,"ps->threadID,ps->il,stDip %7d %7d %7.3f \n",ps->threadID,ps->il,stDip );
   fflush(fp);*/
	vclrf(ps->radon,ndip*lhw[ps->win]*4);
	for(ps->free=0;ps->free<lhw[ps->win];ps->free++) 
	{
	  XijLine(ps);
	  radov(ps);
	}
	if(stDip>0.)
	{
      for(j=0;j<ndip;j++)
      {
		strongD(ps,np,j+j);
		if(j>0) strongD(ps,np,(j+ndip)*2);
	  }
	}
    for(ps->free=0;ps->free<lhw[ps->win];ps->free++)
    {
/*  fprintf(fp,"(2)ps->win,ps->free %7d %7d \n",ps->win,ps->free);
   fflush(fp);*/
      iradon(ps);
      iradov(ps);
      reXijLine(ps);
    }
  }
//  return(ps);
}

void Radon3d(int numcG,int xcG)
{
  int i,k;

/* pthread used  */
  pthread_attr_t  joinAttr;				/* default is wrong so we must set */
  procStruct  *ps;

  pthread_attr_init( &joinAttr );
  pthread_attr_setscope(&joinAttr, PTHREAD_SCOPE_SYSTEM); 
  pthread_attr_setdetachstate( &joinAttr, PTHREAD_CREATE_DETACHED );

	k=0;
  for(i=0;i<numWin;i++)
  {
//      pthread_join( workThreads[k], NULL );	 
	  proc[k].win=i;
 /*  fprintf(fp,"win %6d \n",proc[k].win);
   fflush(fp);*/
	  proc[k].threadID=k;
	  proc[k].nc=numcG;
	  proc[k].xc=xcG;
	  procThread2( );
//	  k++;
//	  if(k==numThreads) k=0;
  }
//  for(i=0;i<numThreads;i++) pthread_join( workThreads[i], NULL );
}

void Xij3d(procStruct  *ps)
{
	/* jj(0->numWin-1) kk:line nstart (start point) i(line number) */
  int i2,i,i4,ii,ii4,j,jj4,ll,llh,num1,jj,j4,num;

/*   fprintf(fp,"Xij3d free %6d \n",ps->free);
   fflush(fp);*/
  if(ps->xc==0 || xlSo>0) num=overlap1;
  else num=0;
  i2=ps->free+low[ps->win]-low1;
  num1=ps->win*lhwmax2+ps->free*2;
  llh=i2*numDesiR2*numDesiR2;

  vclrf(ps->xij,numDesiR2*2);
  if(method>0)
  {
  for(i=0;i<numDesiI;i++)
  {
	ii=i*numDesiI;
    for(i4=0;i4<numDesiI;i4++)
    {
	  ii4=ii+i4;
	  for(j=0;j<numDesiI;j++)
	  {
        jj=j*ileng-xlSo+sx1-num+overlap1;
		ll=llh+(ii4)*numDesiR2+j*numDesiI;
	    for(j4=0;j4<numDesiI;j4++)
		{
		  jj4=(jj+j4)*freTrace+num1;
          ps->xij[ii4+ii4]+=buff[jj4]*LLsave[ll+j4];
          ps->xij[ii4+ii4+1]+=buff[jj4+1]*LLsave[ll+j4];
		}
	  }
	}
  }
  }
  else
  {
	  for(j=0;j<numDesiI;j++)
	  {
        jj=j*ileng-xlSo+sx1-num+overlap1;
	    for(j4=0;j4<numDesiI;j4++)
		{
		  jj4=(jj+j4)*freTrace+num1;
		  vmovf(&buff[jj4],&ps->xij[(j*numDesiI+j4)*2],2);
		}
	  }
  }

}

void procThread2()		/* routine threads run */
{
  int i,j,k,i1,i2,npp,k1,k2,nix,j3,j4,j5,j6,j2,i4,j0,i0;
  procStruct  *ps;
  float a2;

  ps = proc;

  vclrf(ps->radon,lhw[ps->win]*ndipI*ndipI*8);
  for(ps->free=0;ps->free<lhw[ps->win];ps->free++) 
  {
    Xij3d(ps);
	radov3d(ps);
  }
  if(stDipI[ps->win]>0.)
  {
    nix=ndipI*ndipI*2;
    npp=ndipI*ndipI*8;
    for(j=0;j<ndipI;j++)
    {
	  j2=j*ndipI*2;
      for(j0=0;j0<ndipI;j0++)
	  {
	    i4=(int)(sqrt((float)(j*j+j0*j0)));
	    if(i4>=ndipI) break;
	    j3=2*(j2+ndipI+j0);
	    j4=2*(j2+nix+j0);
	    j5=2*(j2+j0);
	    j6=2*(j2+nix+ndipI+j0);
		strongD(ps,npp,j3);
		strongD(ps,npp,j4);
		strongD(ps,npp,j5);
		strongD(ps,npp,j6);
	  }
	}
  }

  for(ps->free=0;ps->free<lhw[ps->win];ps->free++)
  {
    iradon3d(ps);
    iradov3d(ps);
  }
//  return(ps);
}

   /* radov.c */
   void radov(procStruct *ps)
   {
     int i,i1,i2,j,j1,j2,m1,nd,ndip,*diplim,numD;
     float a1,a2,ac,as,ar,ai,dp;

/*  fprintf(fp,"enter radov free %7d \n",ps->free);
   fflush(fp);*/
  if(sign==0) 
  {
	ndip=ndipI;
	diplim=diplimI;
	dp=dpI;
	numD=numDesiI;
  }
  else 
  {
	ndip=ndipX;
	diplim=diplimX;
	dp=dpX;
	numD=numDesiX;
  }
//  fprintf(fp,"radov free,sign,ndip,numD,diplim[0],dp %7d %7d %7d %7d %7d %8.5f \n",ps->free,sign,ndip,numD,diplim[0],dp);
//   fflush(fp);
	 i2=ps->free+low[ps->win];
     a1=2.*3.1415926/lenft*i2*dp;
	 m1=ps->free*ndip*2;
        
	 nd=diplim[i2-low1];

   /* start each trace */                            
     for(i=0;i<numD;i++)
     {                   
   /* start each dip */  
	   i1=i+i;
       for(j=0;j<nd;j++)
       { 
		 j1=2*(m1+j);
		 j2=2*(m1+ndip+j);
		 a2=j*i*a1;
         ac=cos(a2);
         as=sin(a2);
         ar=ps->xij[i1]*ac-ps->xij[i1+1]*as;
         ai=ps->xij[i1]*as+ps->xij[i1+1]*ac;
         ps->radon[j1]+=ar;
         ps->radon[j1+1]+=ai;
		 if(j>0)
		 {
           ar=ps->xij[i1]*ac+ps->xij[i1+1]*as;
           ai=-ps->xij[i1]*as+ps->xij[i1+1]*ac;
           ps->radon[j2]+=ar;
           ps->radon[j2+1]+=ai;
		 }
       }      
     }
   }

   /* iradon.c */
   void iradon(procStruct *ps)
   {
     int i,i2,ndip;
     float te,a1,resolu;
     
/*  fprintf(fp,"enter iradon free %7d \n",ps->free);
   fflush(fp);*/
  if(sign==0) 
  {
	ndip=ndipI;
	resolu=resoluI[ps->win];
  }
  else 
  {
	ndip=ndipX;
	resolu=resoluX[ps->win];
  }
//  fprintf(fp,"iradon free,sign,ndip,resolu %7d %7d %7d %8.5f \n",ps->free,sign,ndip,resolu);
//   fflush(fp);
	 if(method==1 && resolu==0.) return;
   /* winfil=4 omeg =flat */ 
	 i2=ps->free+low[ps->win];
     te=(float)i2/(float)lenft;

	 if(method==1) a1=resolu;
	 else a1=1.0+resolu;
     te=pow(te,a1);
   
   /* radon*omega */     
     for(i=4*ndip*ps->free;i<4*ndip*(ps->free+1);i++)  
     { 
	   ps->radon[i]*=te;
	 }
   }

   /* iradov.c */
   void iradov(procStruct *ps)
   {
     int i,i1,i2,j,nd,m1,j1,j2,ndip,*diplim,numD;
     float a1,a2,ac,as,ar,ai,dp;
                                                     
/*  fprintf(fp,"enter iradov free %7d \n",ps->free);
   fflush(fp);*/
  if(sign==0) 
  {
	ndip=ndipI;
	diplim=diplimI;
	dp=dpI;
	numD=numDesiI;
  }
  else 
  {
	ndip=ndipX;
	diplim=diplimX;
	dp=dpX;
	numD=numDesiX;
  }
//  fprintf(fp,"iradov free,sign,ndip,numD,diplim[0],dp %7d %7d %7d %7d %7d %8.5f \n",ps->free,sign,ndip,numD,diplim[0],dp);
//   fflush(fp);
	 i2=ps->free+low[ps->win];
     a1=2.*3.1415926/lenft*i2*dp;
	 m1=ps->free*ndip*2;
	 nd=diplim[i2-low1];

	 vclrf(ps->xij,numD*2);
	 for(i=0;i<numD;i++)
     {                   
	   i1=i+i;
       for(j=0;j<nd;j++)
       {                      
		 j1=2*(m1+j);
		 j2=2*(m1+ndip+j);
		 a2=j*i*a1;
         ac=cos(a2);
         as=sin(a2);
         ar=ps->radon[j1]*ac+ps->radon[j1+1]*as;
         ai=-ps->radon[j1]*as+ps->radon[j1+1]*ac;
         ps->xij[i1]+=ar;
         ps->xij[i1+1]+=ai;
		 if(j>0)
		 {
           ar=ps->radon[j2]*ac-ps->radon[j2+1]*as;
           ai=ps->radon[j2]*as+ps->radon[j2+1]*ac;
           ps->xij[i1]+=ar;
           ps->xij[i1+1]+=ai;
         }      
       }  
     }
   }

/* radov3d.c */
void radov3d(procStruct  *ps)
{
  int i,i1,i2,j,j1,j2,j3,j4,j5,j6,k,k1,i0,j0,nix;
  float a1,a2,ac,as,ar,ai,aii,ajj,*AA;
  int lhit,i4,k2;

/*  fprintf(fp,"radov3d fre %5d  \n",ps->free);
   fflush(fp);*/
  i2=ps->free+low[ps->win];
  a1=2.*3.1415926/lenft*i2*dpI;
  nix=ndipI*ndipI*2;
  k2=diplimI[i2-low1];
  
  AA=&ps->radon[ps->free*ndipI*ndipI*8];

  /* start each trace */                            
  for(j=0;j<k2;j++)
  {
	j2=j*ndipI*2;
    for(j0=0;j0<k2;j0++)
	{
	  i4=(int)(sqrt((float)(j*j+j0*j0)));
	  if(i4>=k2) break;
	  j3=2*(j2+ndipI+j0);
	  j4=2*(j2+nix+j0);
	  j5=2*(j2+j0);
	  j6=2*(j2+nix+ndipI+j0);
      for(i=0;i<numDesiI;i++)
      {  
        j1=i*numDesiI;
	    aii=(i-(numDesiI-1)/2)*j;
        for(i0=0;i0<numDesiI;i0++)
        {                   
		  i1=(j1+i0)*2;
		  ajj=(i0-(numDesiI-1)/2)*j0;
          a2=(aii-ajj)*a1;
          ac=cos(a2);
          as=sin(a2);
          ar=ps->xij[i1]*ac-ps->xij[i1+1]*as;
          ai=ps->xij[i1]*as+ps->xij[i1+1]*ac;
          AA[j3]+=ar;
          AA[j3+1]+=ai;
 		  if(j>0 || j0>0)
		  {
            ar=ps->xij[i1]*ac+ps->xij[i1+1]*as;
            ai=-ps->xij[i1]*as+ps->xij[i1+1]*ac;
            AA[j4]+=ar;
            AA[j4+1]+=ai;
		  }
		  if(j>0 || j0>0)
		  {
            a2=(aii+ajj)*a1;
            ac=cos(a2);
            as=sin(a2);
            ar=ps->xij[i1]*ac-ps->xij[i1+1]*as;
            ai=ps->xij[i1]*as+ps->xij[i1+1]*ac;
            AA[j5]+=ar;
            AA[j5+1]+=ai;
            ar=ps->xij[i1]*ac+ps->xij[i1+1]*as;
            ai=-ps->xij[i1]*as+ps->xij[i1+1]*ac;
            AA[j6]+=ar;
            AA[j6+1]+=ai;
          }
		}
      }
	}
  }
}

   /* iradon3d.c */
   void iradon3d(procStruct  *ps)
   {
     int i,i2;
     float te,a1;
                   
/*   fprintf(fp,"iradon3d free %6d \n",ps->free);
   fflush(fp);*/
	 if(method==1 && resoluI[ps->win]==0.) return;
   /* winfil=4 omeg =flat */ 
	 i2=ps->free+low[ps->win];
     te=(float)i2/(float)lenft;

	 if(method==1) a1=resoluI[ps->win];
	 else a1=2.0+resoluI[ps->win];
     te=pow(te,a1);
   
   /* radon*omega */     
     for(i=8*ndipI*ndipI*ps->free;i<8*ndipI*ndipI*(ps->free+1);i++)  
     { 
	   ps->radon[i]=ps->radon[i]*te;
	 }
   }

/* iradov3d.c */
void iradov3d(procStruct  *ps)
{
  int i,i0,i1,i2,j,j0,j1,j2,j3,j4,j5,j6,nix,num1,k;
  float a1,a2,ac,as,ar,ai,aii,ajj,*AA;
  int lhit,i4,k2,m1,num,num2;
                
/*  fprintf(fp,"iradov3d fre %5d  \n",ps->free);
   fflush(fp);*/
  if(ps->xc==0 || xlSo>0) num=overlap1;
  else num=0;
  if(ps->xc==0 || xlSo==0) num2=overlap1;
  else num2=0;
  i2=ps->free+low[ps->win];
  a1=2.*3.1415926/lenft*i2*dpI;
  num1=ps->win*lhwmax2+ps->free*2;
  k2=diplimI[i2-low1];
  m1=ps->free*ndipI*ndipI*8;
  nix=ndipI*ndipI*2;
  AA=&ps->radon[ps->free*ndipI*ndipI*8];
        
  for(j=0;j<k2;j++)
  {
	j2=j*ndipI*2;
    for(j0=0;j0<k2;j0++)
	{
	  k=(int)(sqrt((float)(j*j+j0*j0)));
	  if(k>=k2) break;
	  j3=2*(j2+ndipI+j0);
	  j4=2*(j2+nix+j0);
	  j5=2*(j2+j0);
	  j6=2*(j2+nix+ndipI+j0);
      for(i=outLine;i<numDesiI-overlap1;i++)
      {                   
        j1=(i-outLine)*icrb;
	    aii=(i-(numDesiI-1)/2)*j;
	    for(i0=ouStart;i0<numDesiI-overlap1;i0++)
        {   
		  i1=(j1+i0+sx1-num+overlap1-num2)*freTrace+num1;
		  ajj=(i0-(numDesiI-1)/2)*j0;
          a2=(aii-ajj)*a1;
          ac=cos(a2);
          as=sin(a2);
          ar=AA[j3]*ac+AA[j3+1]*as;
          ai=-AA[j3]*as+AA[j3+1]*ac;
          buff1[i1]+=ar;
          buff1[i1+1]+=ai;
	      if(j>0 || j0>0)
	      {
            ar=AA[j4]*ac-AA[j4+1]*as;
            ai=AA[j4]*as+AA[j4+1]*ac;
            buff1[i1]+=ar;
            buff1[i1+1]+=ai;
	      }

		  if(j>0 || j0>0)
		  {
            a2=(aii+ajj)*a1;
            ac=cos(a2);
            as=sin(a2);
            ar=AA[j5]*ac+AA[j5+1]*as;
            ai=-AA[j5]*as+AA[j5+1]*ac;
            buff1[i1]+=ar;
            buff1[i1+1]+=ai;
            ar=AA[j6]*ac-AA[j6+1]*as;
            ai=AA[j6]*as+AA[j6+1]*ac;
            buff1[i1]+=ar;
            buff1[i1+1]+=ai;
		  }
        }
      }
	}
  }
}


/*get_LLH1 */
void get_LLH1(procStruct *ps)
{
  int i,j,k,k1,i2,i1,j1,i3;
  int ii1,jj1,np;
  float a1,a2,a3,*AA,*g,*BB;

  a1=2.*3.1415926/lenft;
  i2=ps->free+low1;
  a1*=dpI*i2*0.5;
  AA=&LLsave[ps->free*numDesiR2*numDesiR2];
  g=ps->xij;
  BB=ps->radon;
  np=ndipI+ndipI-1;

/*  fprintf(fp,"get_LLH1 ps->free,i2,np,ndip,dp %7d %7d %7d %7d %8.4f  \n",ps->free,i2,np,ndip,dp);
  fflush(fp);*/
  g[0]=np;
  for(i=1;i<numDesiI;i++)
  {
	a2=i*a1;
	g[i]=sin(a2*np)/sin(a2);
  }

/*  fprintf(fp,"g %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n",g[0],g[1],g[2],g[3],g[4],g[5],g[6]);
  fflush(fp);*/
/* first inline squre LLH total numDesign */
  for(i=0;i<numDesiI;i++)
  {
	i1=i*numDesiI;
    for(j=0;j<numDesiI;j++) /* inside LLH */
    {
      a3=g[i]*g[j];
	  j1=i1+j;
	  AA[j1]=a3;
	  if(j>0)
	  {
	    jj1=j*numDesiR2+i1;
	    AA[jj1]=a3;
	  }
      for(k=j+1;k<numDesiI;k++)
	  {
		k1=numDesiR2*(k-j)+k-j;
		AA[j1+k1]=a3;
		if(j>0) AA[jj1+k1]=a3;
	  }
	}
  }

/* copy each toeplitz */
  i3=numDesiR2*numDesiI;
  for(i=0;i<numDesiI;i++)
  {
	i1=i*numDesiI;
	ii1=i*i3;
	if(i>0)
	{
	  for(k=0;k<numDesiI;k++)
	  {
		k1=k*numDesiR2;
		vmovf(&AA[i1+k1],&AA[ii1+k1],numDesiI);
	  }
	}
    for(j=i+1;j<numDesiI;j++)
	{
	  j1=j*(i3+numDesiI);
	  for(k=0;k<numDesiI;k++)
	  {
		k1=k*numDesiR2;
		vmovf(&AA[i1+k1],&AA[j1+k1-ii1],numDesiI);
		if(i>0)	vmovf(&AA[i1+k1],&AA[j1+k1-i1],numDesiI);
	  }
	}
  }

/*  fprintf(fp,"before inverse\n");
  fflush(fp);*/
  for(i=0;i<numDesiR2;i++)
	  AA[i*(numDesiR2+1)]*=(1.0+numDesiI*0.08/(float)i2)*(1.0+numDesiI*0.08/(float)i2);
  ludcmp(AA,numDesiR2,idx,vv);
  for(i=0;i<numDesiR2;i++)
  {
	vclrf(g,numDesiR2);
	g[i]=1.;
    lubksb(AA,numDesiR2,idx,g);
	for(j=0;j<numDesiR2;j++) BB[j*numDesiR2+i]=g[j];
  }
  vmovf(BB,AA,numDesiR2*numDesiR2);

/*  fprintf(fp,"after inverse\n");
  fflush(fp);*/
}

/* get_LLH */
void get_LLH(procStruct *ps)
{
  int i,j,k,k1,i1,i2,np;
  float a1,a2,a[60],r[60],*g,*f,*AA;

  AA=&LLsave[ps->free*numDesiR2];
  a1=2.*3.1415926/lenft;
  i2=ps->free+low1;
  a1*=dpI*i2*0.5;
  g=ps->radon;
  f=ps->xij;
  np=ndipI+ndipI-1;

/*  fprintf(fp,"get_LLH ps->free,i2,np,ndip,dp %7d %7d %7d %7d %8.4f  \n",ps->free,i2,np,ndip,dp);
  fflush(fp);*/
  r[0]=np*(1+0.08*numDesiI/(float)i2);
  for(i=1;i<numDesiI;i++)
  {
	a2=i*a1;
	r[i]=sin(a2*np)/sin(a2);
  }

/*  fprintf(fp,"r %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f \n",r[0],r[1],r[2],r[3],r[4],r[5],r[6]);
  fflush(fp);*/
  for(i=0;i<numDesiI;i++)
  {
    vclrf(g,numDesiI);
    g[i]=1.;
    stoepf(numDesiI,r,g,f,a);
	for(j=0;j<numDesiI;j++) AA[j*numDesiI+i]=f[j];
  }
}

void procThreadLL()		/* routine threads run */
{
  int i,j,k;
  procStruct  *ps;

  ps = proc;

  if(radonF==4)
  {
    get_LLH1(ps);
  }
  else
  {
    get_LLH(ps);
  }
//  return(ps);
}

void get_LLsave()
{
  int i,j,k;
  pthread_attr_t  joinAttr;				/* default is wrong so we must set */
  procStruct  *ps;

  pthread_attr_init( &joinAttr );
  pthread_attr_setscope(&joinAttr, PTHREAD_SCOPE_SYSTEM); 
  pthread_attr_setdetachstate( &joinAttr, PTHREAD_CREATE_DETACHED );

  k=0;
  for(i=0;i<lhwmax;i++)
  {
 //   pthread_join( workThreads[k], NULL );	 
	proc[k].threadID=k;
    proc[k].free=i;
/*  fprintf(fp,"i,k %7d %7d  \n",i,k);
  fflush(fp);*/
	procThreadLL();
//	k++;
//	if(k==numThreads) k=0;
  }
//  for(i=0;i<numThreads;i++) pthread_join( workThreads[i], NULL );
}

void caddr(float *trace,float *trace1,float *a,int j)
{
  int i,i1,k,k1;

  /* real of conjugate of t1*t */
  a[0]=a[1]=0.;
  for(i=0;i<numDesiF-numFilt+1;i++)
  {
    i1=i*numDesiF;
    for(k=0;k<j;k++)
	{
	  k1=2*(i1+k);
      a[0]=a[0]+trace[k1]*trace1[k1]+trace[k1+1]*trace1[k1+1];
	}
  }
}
   
void caddc(float *trace,float *trace1,float *a,int j)
{
  int i,i1,k,k1;

  /* complex of conjugate of t1*t */
  a[0]=a[1]=0.;
  for(i=0;i<numDesiF-numFilt+1;i++)
  {
    i1=i*numDesiF;
    for(k=0;k<j;k++)
    {
	  k1=2*(i1+k);
      a[0]=a[0]+trace[k1]*trace1[k1]+trace[k1+1]*trace1[k1+1];
      a[1]=a[1]-trace[k1]*trace1[k1+1]+trace1[k1]*trace[k1+1];
	}
  }
}
   
void caddi(float *trace,float *trace1,float *a,int j)
{
  int i,i1;

  /* complex of t1*t */
  a[0]=a[1]=0.;
  for(i=0;i<j;i++)
  {
	i1=2*i;
    a[0]=a[0]+trace[i1]*trace1[i1]-trace[i1+1]*trace1[i1+1];
    a[1]=a[1]+trace[i1]*trace1[i1+1]+trace1[i1]*trace[i1+1];
  }
}
   
void caddir(float *trace,float *trace1,float *a,int j)
{
  int i,i1;

  /* complex of conjugate of t1*t */
  a[0]=a[1]=0.;
  for(i=0;i<j;i++)
  {
	i1=2*i;
    a[0]=a[0]+trace[i1]*trace1[i1]+trace[i1+1]*trace1[i1+1];
    a[1]=a[1]-trace[i1]*trace1[i1+1]+trace1[i1]*trace[i1+1];
  }
}
   
void vclrf(float *trace,int j)
{
  int i;

  for(i=0;i<j;i++)
    trace[i]=0.0;
}
   
void vclri(int *trace,int j)
{
  int i;

  for(i=0;i<j;i++)
    trace[i]=0;
}
   
void vmovf(float *trace,float *trace1,int j)
{
  int i;

  for(i=0;i<j;i++)
    trace1[i]=trace[i];
}
   
void vmovi(int *trace,int *trace1,int j)
{
  int i;

  for(i=0;i<j;i++)
    trace1[i]=trace[i];
}
     
void vabsf(float *trace,int j,int k)
{
  int i,i2;
  float a1,a2;

  a1=0.;
  for(i=j;i<=k;i++) a1+=fabs(trace[i]);
  a2=absR[cdp];
  i2=1;
  if(a1>0.)
  {
	if(xl>0 && absR[cdp-1]>0.)
	{
	  i2++;
	  a2+=absR[cdp-1];
	}
	if(xl<icrb-1 && absR[cdp+1]>0.)
	{
	  i2++;
	  a2+=absR[cdp+1];
	}
	if(il>0 && absR[cdp-icrb]>0.)
	{
	  i2++;
	  a2+=absR[cdp-icrb];
	}
	if(il<nlines-1 && absR[cdp+icrb]>0.)
	{
	  i2++;
	  a2+=absR[cdp+icrb];
	}
	a1=a2/(a1*(float)i2);
	for(i=j;i<=k;i++) trace[i]*=a1;
  }
}
     
void  get_bin( float *trace, int seq, int num,int bin )
{
int  i, n,i1,i2;

//   fprintf(fp,"get_bin seq,num,bin %7d %7d %7d   \n",seq,num,bin);
//   fflush(fp);
   i = seq / tracesPerBin;             /* calc file */
   seq -= i * tracesPerBin;            /* trace in file */
//   fprintf(fp,"get_bin seq,num,bin,i %7d %7d %7d  %7d  \n",seq,num,bin,i);
//   fflush(fp);
   if(seq+num>tracesPerBin)
   {
	 i1=tracesPerBin-seq;
	 i2=num-i1;
     fseek( binF[bin][i], seq*shtotal*4, SEEK_SET );
     n = fread( trace, sizeof( float ), shtotal*i1, binF[bin][i] );
     if( n != shtotal*i1 )
       ProgramErr( "get_bin -- error reading bin i1" );
     fseek( binF[bin][i+1], 0, SEEK_SET );
     n = fread( &trace[shtotal*i1], sizeof( float ), shtotal*i2, binF[bin][i+1] );
     if( n != shtotal*i2 )
       ProgramErr( "get_bin -- error reading bin i2" );
   }
   else
   {
     fseek( binF[bin][i], seq*shtotal*4, SEEK_SET );
     n = fread( trace, sizeof( float ), shtotal*num, binF[bin][i] );
     if( n != shtotal*num )
       ProgramErr( "get_bin -- error reading bin num" );
   }
/*   for(i=0;i<num;i++)
   {
	 fprintf(fp,"%6.1f",trace[i*shtotal]);
	 if(i%10==9) fprintf(fp,"  %7d \n",i);
   }
   fprintf(fp,"  %7d \n",i);*/
}

void  put_bin( float *trace, int seq, int num,int bin )
{
int  i, n,i1,i2;

//   fprintf(fp,"put_bin seq,num,bin %7d %7d %7d   \n",seq,num,bin);
//   fflush(fp);
   i = seq / tracesPerBin;             /* calc file */
   seq -= i * tracesPerBin;            /* trace in file */
   if(seq+num>tracesPerBin)
   {
	 i1=tracesPerBin-seq;
	 i2=num-i1;
     fseek( binF[bin][i], seq*shtotal*4, SEEK_SET );
     n = fwrite( trace, sizeof( float ), shtotal*i1, binF[bin][i] );
     if( n != shtotal*i1 )
       ProgramErr( "put_bin -- error whiting bin i1" );
     fseek( binF[bin][i+1], 0, SEEK_SET );
     n = fwrite( &trace[shtotal*i1], sizeof( float ), shtotal*i2, binF[bin][i+1] );
     if( n != shtotal*i2 )
       ProgramErr( "put_bin -- error whiting bin i2" );
   }
   else
   {
     fseek( binF[bin][i], seq*shtotal*4, SEEK_SET );
     n = fwrite( trace, sizeof( float ), shtotal*num, binF[bin][i] );
     if( n != shtotal*num )
       ProgramErr( "put_bin -- error writing bin num" );
   }
/*   for(i=0;i<num;i++)
   {
	 fprintf(fp,"%6.1f",trace[i*shtotal]);
	 if(i%10==9) fprintf(fp,"  %7d \n",i);
   }
   fprintf(fp,"  %7d \n",i);*/
}

void  get_ampl( float *trace, int seq, int num )
{
int  i, n,i1,i2;

//	  fprintf(fp,"get_ampl seq,num %5d %5d \n",seq,num );
//   fflush(fp);
   i = seq / tracesPerAmpl;             /* calc file */
   seq -= i * tracesPerAmpl;            /* trace in file */
   if(seq+num>tracesPerAmpl)
   {
	 i1=tracesPerAmpl-seq;
	 i2=num-i1;
     fseek( amplF[i], seq*sampleAs*4, SEEK_SET );
     n = fread( trace, sizeof( float ), sampleAs*i1, amplF[i] );
     if( n != sampleAs*i1 )
       ProgramErr( "get_ampl -- error reading file i1" );
     fseek( amplF[i+1], 0, SEEK_SET );
     n = fread( &trace[sampleAs*i1], sizeof( float ), sampleAs*i2, amplF[i+1] );
     if( n != sampleAs*i2 )
       ProgramErr( "get_ampl -- error reading file i2" );
   }
   else
   {
     fseek( amplF[i], seq*sampleAs*4, SEEK_SET );
     n = fread( trace, sizeof( float ), sampleAs*num, amplF[i] );
     if( n != sampleAs*num )
       ProgramErr( "get_ampl -- error reading file num" );
   }
}

void  put_ampl( float *trace, int seq, int num )
{
int  i, n,i1,i2;

   i = seq / tracesPerAmpl;             /* calc file */
   seq -= i * tracesPerAmpl;            /* trace in file */
//	  fprintf(fp,"put_ampl seq,num %5d %5d \n",seq,num );
//   fflush(fp);
   if(seq+num>tracesPerAmpl)
   {
	 i1=tracesPerAmpl-seq;
	 i2=num-i1;
     fseek( amplF[i], seq*sampleAs*4, SEEK_SET );
     n = fwrite( trace, sizeof( float ), sampleAs*i1, amplF[i] );
     if( n != sampleAs*i1 )
       ProgramErr( "put_ampl -- error whiting file i1" );
     fseek( amplF[i+1], 0, SEEK_SET );
     n = fwrite( &trace[sampleAs*i1], sizeof( float ), sampleAs*i2, amplF[i+1] );
     if( n != sampleAs*i2 )
       ProgramErr( "put_ampl -- error whiting file i2" );
   }
   else
   {
     fseek( amplF[i], seq*sampleAs*4, SEEK_SET );
     n = fwrite( trace, sizeof( float ), sampleAs*num, amplF[i] );
     if( n != sampleAs*num )
       ProgramErr( "put_ampl -- error writing file num" );
   }
}

void  get_head( float *trace, int seq, int num )
{
int  i, n,i1,i2;

//   fprintf(fp,"get_head seq,num,headN[0] %7d %7d %s \n",seq,num,headN[0]);
//   fflush(fp);
   i = seq / tracesPerHead;             /* calc file */
   seq -= i * tracesPerHead;            /* trace in file */
   if(seq+num>tracesPerHead)
   {
	 i1=tracesPerHead-seq;
	 i2=num-i1;
     fseek( headF[i], seq*headers*4, SEEK_SET );
     n = fread( trace, sizeof( float ), headers*i1, headF[i] );
     if( n != headers*i1 )
       ProgramErr( "get_head -- error reading file i1" );
     fseek( headF[i+1], 0, SEEK_SET );
     n = fread( &trace[headers*i1], sizeof( float ), headers*i2, headF[i+1] );
     if( n != headers*i2 )
       ProgramErr( "get_head -- error reading file i2" );
   }
   else
   {
     fseek( headF[i], seq*headers*4, SEEK_SET );
     n = fread( trace, sizeof( float ), headers*num, headF[i] );
     if( n != headers*num )
       ProgramErr( "get_head -- error reading file num" );
   }
//   fprintf(fp,"(1)get_head seq,num,headN[0] %7d %7d %s \n",seq,num,headN[0]);
//   fflush(fp);
}

void  put_head( float *trace, int seq, int num )
{
int  i, n,i1,i2;

//   fprintf(fp,"put_head seq,num,headN[0] %7d %7d %s\n",seq,num,headN[0]);
//   fflush(fp);
   i = seq / tracesPerHead;             /* calc file */
   seq -= i * tracesPerHead;            /* trace in file */
   if(seq+num>tracesPerHead)
   {
	 i1=tracesPerHead-seq;
	 i2=num-i1;
     fseek( headF[i], seq*headers*4, SEEK_SET );
     n = fwrite( trace, sizeof( float ), headers*i1, headF[i] );
     if( n != headers*i1 )
       ProgramErr( "put_head -- error whiting file i1" );
     fseek( headF[i+1], 0, SEEK_SET );
     n = fwrite( &trace[headers*i1], sizeof( float ), headers*i2, headF[i+1] );
     if( n != headers*i2 )
       ProgramErr( "put_head -- error whiting file i2" );
   }
   else
   {
     fseek( headF[i], seq*headers*4, SEEK_SET );
     n = fwrite( trace, sizeof( float ), headers*num, headF[i] );
     if( n != headers*num )
       ProgramErr( "put_head -- error writing file num" );
   }
}

void  get_trace( float *trace, int seq, int num,int kk )
{
int  i, n,i1,i2;

 //  fprintf(fp,"get_trace seq,num,kk %7d %7d  %7d  \n",seq,num,kk);
 //  fflush(fp);
   i = seq / tracesPerFile;             /* calc file */
   seq -= i * tracesPerFile;            /* trace in file */
   if(seq+num>tracesPerFile)
   {
	 i1=tracesPerFile-seq;
	 i2=num-i1;
	 if(kk==0)
	 {
       fseek( tracF[i], seq*freTrace*4, SEEK_SET );
       n = fread( trace, sizeof( float ), freTrace*i1, tracF[i] );
       if( n != freTrace*i1 )
         ProgramErr( "get_trace trac -- error reading file i1" );
       fseek( tracF[i+1], 0, SEEK_SET );
       n = fread( &trace[freTrace*i1], sizeof( float ), freTrace*i2, tracF[i+1] );
       if( n != freTrace*i2 )
         ProgramErr( "get_trace trac -- error reading file i2" );
	 }
	 else
	 {
       fseek( trarF[i], seq*freTrace*4, SEEK_SET );
       n = fread( trace, sizeof( float ), freTrace*i1, trarF[i] );
       if( n != freTrace*i1 )
         ProgramErr( "get_trace trar -- error reading file i1" );
       fseek( trarF[i+1], 0, SEEK_SET );
       n = fread( &trace[freTrace*i1], sizeof( float ), freTrace*i2, trarF[i+1] );
       if( n != freTrace*i2 )
         ProgramErr( "get_trace trar -- error reading file i2" );
	 }
   }
   else
   {
	 if(kk==0)
	 {
       fseek( tracF[i], seq*freTrace*4, SEEK_SET );
       n = fread( trace, sizeof( float ), freTrace*num, tracF[i] );
       if( n != freTrace*num )
         ProgramErr( "get_trace trac -- error reading file num" );
	 }
	 else
	 {
       fseek( trarF[i], seq*freTrace*4, SEEK_SET );
       n = fread( trace, sizeof( float ), freTrace*num, trarF[i] );
       if( n != freTrace*num )
         ProgramErr( "get_trace trar -- error reading file num" );
	 }
   }
}

void  put_trace( float *trace, int seq, int num,int kk)
{
int  i, n,i1,i2;

//   fprintf(fp,"put_trace seq,num,kk %7d %7d  %7d  \n",seq,num,kk);
//   fflush(fp);
   i = seq / tracesPerFile;             /* calc file */
   seq -= i * tracesPerFile;            /* trace in file */
   if(seq+num>tracesPerFile)
   {
	 i1=tracesPerFile-seq;
	 i2=num-i1;
	 if(kk==0)
	 {
	   fseek( tracF[i], seq*freTrace*4, SEEK_SET );
       n = fwrite( trace, sizeof( float ), freTrace*i1, tracF[i] );
       if( n != freTrace*i1 )
         ProgramErr( "put_trace trac -- error writing file i1" );
       fseek( tracF[i+1], 0, SEEK_SET );
       n = fwrite( &trace[freTrace*i1], sizeof( float ), freTrace*i2, tracF[i+1] );
       if( n != freTrace*i2 )
         ProgramErr( "put_trace trac-- error writing file i2" );
	 }
	 else
	 {
	   fseek( trarF[i], seq*freTrace*4, SEEK_SET );
       n = fwrite( trace, sizeof( float ), freTrace*i1, trarF[i] );
       if( n != freTrace*i1 )
         ProgramErr( "put_trace trar -- error writing file i1" );
       fseek( trarF[i+1], 0, SEEK_SET );
       n = fwrite( &trace[freTrace*i1], sizeof( float ), freTrace*i2, trarF[i+1] );
       if( n != freTrace*i2 )
         ProgramErr( "put_trace trar -- error writing file i2" );
	 }
   }
   else
   {
	 if(kk==0)
	 {
       fseek( tracF[i], seq*freTrace*4, SEEK_SET );
       n = fwrite( trace, sizeof( float ), freTrace*num, tracF[i] );
       if( n != freTrace*num )
         ProgramErr( "put_trace trac -- error writing file num" );
	 }
	 else
	 {
       fseek( trarF[i], seq*freTrace*4, SEEK_SET );
       n = fwrite( trace, sizeof( float ), freTrace*num, trarF[i] );
       if( n != freTrace*num )
         ProgramErr( "put_trace trar -- error writing file num" );
	 }
   }
}

void  OpenTracFile( void )
{
int  i,n;
char  tempStr[8];

//   strcpy( tracN[numFiles], "/scratch/ssa101/exscratch1/" ); //Crusher
//   strcpy( tracN[numFiles], "/scratch/ds01/exscr1/" );  
   strcpy( tracN[numFiles], "/scratch/" );  
//   strcpy( tracN[numFiles], "/scratch/ds4/chengtest/" );  
   strcat( tracN[numFiles], scratchName[0] );
   sprintf( tempStr, "%d", myid );
   strcat( tracN[numFiles], tempStr );
   strcat( tracN[numFiles], ".trac" );
   sprintf( tempStr, "%d", numFiles );
   strcat( tracN[numFiles], tempStr );

   tracF[numFiles] = fopen( tracN[numFiles], "wb+" );
   if( tracF[numFiles] == 0 )
      ProgramErr( "error opening trac file" );

   if((fxyF && radonF) || radonF==3)
   {
//   strcpy( trarN[numFiles], "/scratch/ssa101/exscratch1/" ); //Crusher
//   strcpy( trarN[numFiles], "/scratch/ds01/exscr1/" );  
   strcpy( trarN[numFiles], "/scratch/" );  
//   strcpy( trarN[numFiles], "/scratch/ds4/chengtest/" );  
   strcat( trarN[numFiles], scratchName[0] );
   sprintf( tempStr, "%d", myid );
   strcat( trarN[numFiles], tempStr );
   strcat( trarN[numFiles], ".trar" );
   sprintf( tempStr, "%d", numFiles );
   strcat( trarN[numFiles], tempStr );

   trarF[numFiles] = fopen( trarN[numFiles], "wb+" );
   if( trarF[numFiles] == 0 )
      ProgramErr( "error opening trar file" );

   }
   if(myid==1)
   printf("open trace file number = %7d \n",numFiles);
   numFiles++;
}

void  OpenAmplFile( void )
{
int  i,n;
char  tempStr[8];

   strcpy( amplN[numAmpls], "/scratch/" ); 
   strcat( amplN[numAmpls], scratchName[0] );
   sprintf( tempStr, "%d", myid );
   strcat( amplN[numAmpls], tempStr );
   strcat( amplN[numAmpls], ".ampl" );
   sprintf( tempStr, "%d", numAmpls );
   strcat( amplN[numAmpls], tempStr );
   amplF[numAmpls] = fopen( amplN[numAmpls], "wb+" );

   if( amplF[numAmpls] == 0 )
      ProgramErr( "error opening ampl file" );
   if(myid==1)
   printf("open ampl file number = %7d \n",numAmpls);
   numAmpls++;
}

void  OpenHeadFile( void )
{
int  i,n;
char  tempStr[8];

   strcpy( headN[numHeads], "/scratch/" );  
   strcat( headN[numHeads], scratchName[0] );
   sprintf( tempStr, "%d", myid );
   strcat( headN[numHeads], tempStr );
   strcat( headN[numHeads], ".head" );
   sprintf( tempStr, "%d", numHeads );
   strcat( headN[numHeads], tempStr );
   headF[numHeads] = fopen( headN[numHeads], "wb+" );

   if( headF[numHeads] == 0 )
      ProgramErr( "error opening head file" );
   if(myid==1)
   printf("open header file number = %7d \n",numHeads);
   numHeads++;
}

void  OpenBinFile()
{
int  i,n,num;
char  tempStr[8];

for(num=0;num<nbin;num++)
{
/*fprintf(fp,"numFiles %5d \n",numFiles);
	  fflush(fp);*/

//   if(num<10) strcpy( binN[num][numBins], "/scratch/ssa101/exscratch1/" ); //Crusher
//   else strcpy( binN[num][numBins], "/scratch/ssa102/exscratch2/" ); //Crusher
//   strcpy( binN[num][numBins], "/scratch/ds4/chengtest/" ); //Crusher
//   strcpy( binN[num][numBins], "/scratch/ds01/exscr1/" );  
   strcpy( binN[num][numBins], "/scratch/" );  
   strcat( binN[num][numBins], scratchName[0] );
   sprintf( tempStr, "%d", myid );
   strcat( binN[num][numBins], tempStr );
   strcat( binN[num][numBins], ".bin" );
   i=num+1;
   if(i<10) strcat( binN[num][numBins], "0" );
   sprintf( tempStr, "%d", i );
   strcat( binN[num][numBins], tempStr );
   sprintf( tempStr, "%d", numBins );
   strcat( binN[num][numBins], tempStr );
   binF[num][numBins] = fopen( binN[num][numBins], "wb+" );

   if( binF[num][numBins] == 0 )
      ProgramErr( "error opening bin file" );
}
   if(myid==1)
   printf("open bin file number = %7d \n",numBins);
numBins++;
}

void stoepf (int n, float *r, float *g, float *f, float *a)
/*****************************************************************************
Solve a symmetric Toeplitz linear system of equations Rf=g for f
(float version) 
******************************************************************************
Input:
n		dimension of system
r		array[n] of top row of Toeplitz matrix
g		array[n] of right-hand-side column vector

Output:
f		array[n] of solution (left-hand-side) column vector
a		array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
******************************************************************************
Notes:
This routine does NOT solve the case when the main diagonal is zero, it
just silently returns.

The left column of the Toeplitz matrix is assumed to be equal to the top
row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int i,j;
	float v,e,c,w,bot;

	if (r[0] == 0.0) return;

	a[0] = 1.0;
	v = r[0];
	f[0] = g[0]/r[0];

	for (j=1; j<n; j++) {
		
		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
		a[j] = 0.0;
		f[j] = 0.0;
		for (i=0,e=0.0; i<j; i++)
			e += a[i]*r[j-i];
		c = e/v;
		v -= c*e;
		for (i=0; i<=j/2; i++) {
			bot = a[j-i]-c*a[i];
			a[i] -= c*a[j-i];
			a[j-i] = bot;
		}

		/* use a and v above to get f[i], i = 0,1,2,...,j */
		for (i=0,w=0.0; i<j; i++)
			w += f[i]*r[j-i];
		c = (w-g[j])/v;
		for (i=0; i<=j; i++)
			f[i] -= c*a[j-i];
	}
}

void GetSocketData( char *theData, int length )
{
   int  bytesRead = 0;
   int  n;

   while( length > 0 )
      {
      n = read( socketFile, theData+bytesRead, length );
      length -= n;
      bytesRead += n;
      }
}

void ProgramErr( char *errString )
{
	int result;
   printf( "%s\n", errString );

    result = -1;
    write( socketFile, ( char* )&result, 4 );
	exit( -1 );
}

void  PrintTime( void )
{
time_t  systime;
char  *s;

	systime = time( NULL );
	s = ctime( &systime );
	s[20] = 0;
	printf( "   %s ", s );
}

