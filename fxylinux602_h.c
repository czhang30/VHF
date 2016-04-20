/* program fxy.c */
#include  "cpromax.h"
#include  "cglobal.h"
#include  "cSocketTool.h" 

#include  <time.h>

/*
// internet socket stuff
*/
#define u_short unsigned short
#define u_long unsigned long
#define ulong unsigned long
#define u_char unsigned char
#define u_int unsigned int

#include  <signal.h>
#include  <sys/types.h>
#include  <sys/socket.h>
#include  <netdb.h>
#include  <netinet/in.h>
#include  <arpa/inet.h>

int  serverFd, clientFd,sockPort, result;
/*unsigned long clientLen;*/
int clientLen;
struct sockaddr_in  clientAddr, serverAddr;

float  *trace;            /* input trace */
int  *ithdr;                      /* integer trace header */
float  *rthdr;                    /* float trace header */
float  sampleRate;                /* sample rate of data */
int  samples;                     /* number of samples in data */
int  headers;                     /* number of entries in trace header */
float *traceIn;                    /* work trace */

int  indexCdp;                    /* cdp header */
int  indexTrc_type;                    /* trc_type header */
int  indexpstmBin;                /* location of pstm offset bin number header */
int  indexEnd_ens;                /* end os ensamples flag */
int  indexSeqno;                  

int  minCDP;                      /* min CDP in data */
int  maxCDP;                      /* max CDP in data */
int  minIlin;                     /* first inline */
int  maxIlin;                     /* last inline */
int  minXlin;                     /* first xline */
int  maxXlin;                     /* last xline */
int  ntrcdp;                      /* max CDP fold */

int  ncrbi;                       /* number of cdps inline */
int  icrb;                        /* number of cdps inline to fxy */
int  nlines;                      /* number of cdps xline to fxy */
int idmap;                        /* icrb*nlines */

int   lmin;                       /* minimum inline no. to fxy */
int   lmax;                       /* maximum inline no. to fxy */
int   xmin;                       /* minimum xline no. to fxy */
int   xmax;                       /* maximum xline no. to fxy */
int   tleng;                      /* length of operator time window default = 1000 ms */
int   taper;                      /* time window taper, default = 200 ms */
int  radonF;                      /* high dip filter option */
int  numDesiI,numDesiX;                  /* number of lines in filter design default=25 */
float hiradon;                   /* higher radon transform accurace */
int  dipI,dipX;                         /* dip ms/trace */
int  method;
int timeseq[7];
float stDipI[7],stDipX[7];
float resoluI[7],resoluX[7];
float qmax[7],qmaxm[7],qminp[7],qmin[7];  
int  numDesiF;                    /* number of lines in filter design default=25 */
int  numFilt;                     /* number of lines in filter default=5 */
int  fxyF;                        /* fxy decon option */
int  overlapF,overlapI,overlapX;  /* overlap between filter calculating */
int OuPut,nbin;
int dbaseParm;
int AmpBa,ixSame; /* 2004,11 */
float AmpParm;
int FreBa; /* 2004,12 */
float FreParm;

FILE *fp; 

/* 02,06,11 */
int idtyp;                        /* input data type */
int bin;
int shtotal,shtotal1,shtotal4;
int traceCount;
int numGroup;

void  GetSocketData( char *theData, long length );
void  ConnectRemote();
void  GetParameters();
void  PutGlobals();
void  PreFxy();
void binStart();
void  Fxy();
void  traceDisk(int kk);
void vclrf(float *trace,int i);
void vclri(int *trace,int i);
void vmovf(float *trace,float *trace1,int i);
void vmovi(int *trace,int *trace1,int i);
void Perform();

int  main( int ac, char **av )
{
int  i;
   initStandAlone( ac, av, "SOCKET_TOOL" );
   stConnectToServer();
   stSetToolName( "fxypower" );

/* check input dataset */
   indexCdp = stHdrIndex( "cdp");
   indexTrc_type = stHdrIndex( "trc_type");

   GetParameters();
   PreFxy();
   ConnectRemote();
   printf( "server -- passing globals\n" );
   PutGlobals();
   stEndInitialization();
   Fxy();

   stCloseSocketLink();
}             

void  GetParameters( void )
{
  int i,j,lmin1,lmax1,xmin1,xmax1,dbaseParm1,numDesi;
  char *token,*token1;
  char sep0[]="-(),/";
  char sep1[]="0123456789.";
  char *stimeseq;                       /* string with time position */
  char *sqmaxm,*sqmin,*sqminp,*sqmax;                        /* string with smooth */
  char *sresoluI,*sresoluX;                       /* string with degrees after migration in depth domain */
  char *sstDipI,*sstDipX;                     /* string with migration interval */

/*   fp=fopen("/home/cheng/72/advance/port/src/exe/fxylinux/arg","w"); */

  uParGetInt( "sockPort", &sockPort );  /* intenet code */
  uParGetInt( "dbaseParm", &dbaseParm1 );          
  if(dbaseParm1)
  {
    i = stHdrExists( "CDP" );
    if( i!=1 ) stErrFatal( "error reading CDP header" );
    else indexCdp = stHdrIndex( "cdp");
    if( globalRuntime->idomain != ITX )
      stErrFatal( "Domain must be time-space" );

    idtyp = globalRuntime->idtyp;         /* 0=unstacked 7=stacked */
    sampleRate = globalRuntime->samprat;
    samples    = globalRuntime->numsmp;
    headers    = globalRuntime->nth;
    minIlin = globalGeom->minilin;        /* first inline */
    maxIlin = globalGeom->maxilin;        /* last inline */
    minXlin = globalGeom->minxlin;        /* first xline */
    maxXlin = globalGeom->maxxlin;        /* last xline */
    minCDP  = globalGeom->mincdp;         /* first cdp */
    maxCDP  = globalGeom->maxcdp;         /* last cdp */
    ntrcdp  = globalGeom->ntrcdp;         /* maximum fold of CDP's */

/* get fxy parameters input */
    uParGetInt( "lmin", &lmin );          /* minimum inline no. to fxy default=minXlin */
    uParGetInt( "lmax", &lmax );          /* maximum inline no. to fxy default=maxIlin */
    uParGetInt( "xmin", &xmin );          /* minimum xline no. to fxy default=minXlin */
    uParGetInt( "xmax", &xmax );          /* maximum xline no. to fxy default=minXlin */
  }
  else
  {
	idtyp=1;
    uParGetInt( "samples", &samples );   
    uParGetFloat( "sampleRate", &sampleRate );         
    uParGetInt( "headers", &headers );   
    uParGetInt( "minCDP", &minCDP );         
    uParGetInt( "maxCDP", &maxCDP );         
    uParGetInt( "minIlin", &minIlin );   
    uParGetInt( "maxIlin", &maxIlin );   
    uParGetInt( "minXlin", &minXlin );   
    uParGetInt( "maxXlin", &maxXlin );   
    uParGetInt( "lmin1", &lmin1 );          /* minimum inline no. to fxy default=minXlin */
    uParGetInt( "lmax1", &lmax1 );          /* maximum inline no. to fxy default=maxIlin */
    uParGetInt( "xmin1", &xmin1 );          /* minimum xline no. to fxy default=minXlin */
    uParGetInt( "xmax1", &xmax1 );          /* maximum xline no. to fxy default=minXlin */
	lmin=lmin1;
	lmax=lmax1;
	xmin=xmin1;
	xmax=xmax1;
  }
  if(lmin<minIlin) lmin=minIlin; 
  if(lmax>maxIlin) lmax=maxIlin; 
  if(xmin<minXlin) xmin=minXlin; 
  if(xmax>maxXlin) xmax=maxXlin; 

  uParGetInt( "tleng", &tleng );        /* length of operator time window default = 1000 ms */
  if(tleng<200) tleng=200;
  uParGetInt( "taper", &taper );        /* time window taper, default = 200 ms */
  if(taper>tleng/2) taper=tleng/2;
  uParGetString( "timeseq", &stimeseq );     /* ramp parameters */
  uParGetString( "qmin", &sqmin );     /* ramp parameters */
  uParGetString( "qminp", &sqminp );     /* ramp parameters */
  uParGetString( "qmaxm", &sqmaxm );     /* ramp parameters */
  uParGetString( "qmax", &sqmax );     /* ramp parameters */
/*   timeseq[0]   = ( int )atoi( strtok( stimeseq, "-" ) );
   for(i=1;i<6;i++) timeseq[i] = ( int )atoi( strtok( NULL, "-" ) );
   for(i=1;i<6;i++) if(timeseq[i]==0) break;
   numGroup=i;
   qmin[0]   = ( float )atof( strtok( sqmin, "-" ) );
   for(i=1;i<numGroup;i++) qmin[i] = ( float )atof( strtok( NULL, "-" ) );
   qminp[0] = ( float )atof( strtok( sqminp, "-" ) );
   if(qminp[0]<qmin[0]) qminp[0]=qmin[0]; 
   for(i=1;i<numGroup;i++) 
   {
	 qminp[i] = ( float )atof( strtok( NULL, "-" ) );
     if(qminp[i]<qmin[i]) qminp[i]=qmin[i];
   }
   qmax[0]   = ( float )atof( strtok( sqmax, "-" ) );
   if(qmax[0]>490./sampleRate) qmax[0]=490./sampleRate;
   for(i=1;i<numGroup;i++) 
   {
	 qmax[i] = ( float )atof( strtok( NULL, "-" ) );
     if(qmax[i]>490./sampleRate) qmax[i]=490./sampleRate;
   }
   qmaxm[0]   = ( float )atof( strtok( sqmaxm, "-" ) );
   if(qmaxm[0]>qmax[0]) qmaxm[0]>qmax[0];
   for(i=1;i<numGroup;i++) 
   {
	 qmaxm[i] = ( float )atof( strtok( NULL, "-" ) );
     if(qmaxm[i]>qmax[i]) qmaxm[i]>qmax[i];
   }*/

   vclri(timeseq,6);
   i=0;
   token=strtok(stimeseq,sep0);
   while(token!=NULL)
   {
     timeseq[i]   = ( int )atoi( token );
	 i++;
   	 token=strtok(NULL,sep0);
   }
   numGroup=i;
   vclrf(qmin,6);
   i=0;
   token=strtok(sqmin,sep0);
   while(token!=NULL)
   {
     qmin[i]   = ( float )atof( token );
	 i++;
   	 token=strtok(NULL,sep0);
   }
   vclrf(qminp,6);
   i=0;
   token=strtok(sqminp,sep0);
   while(token!=NULL)
   {
     qminp[i]   = ( float )atof( token );
	 i++;
   	 token=strtok(NULL,sep0);
   }
   vclrf(qmax,6);
   i=0;
   token=strtok(sqmax,sep0);
   while(token!=NULL)
   {
     qmax[i]   = ( float )atof( token );
	 i++;
   	 token=strtok(NULL,sep0);
   }
   vclrf(qmaxm,6);
   i=0;
   token=strtok(sqmaxm,sep0);
   while(token!=NULL)
   {
     qmaxm[i]   = ( float )atof( token );
	 i++;
   	 token=strtok(NULL,sep0);
   }

  numDesiI=numDesiX=numDesiF=9;
  numFilt=3;
  dipI=dipX=4;
  overlapF=overlapI=overlapX=3;
  uParGetInt( "fxyF", &fxyF );          /* doing fxy filter? yes or no */
  if(fxyF>0)
  {
    uParGetInt( "numDesiF", &numDesiF );/* number of lines in filter design default=25 */
    uParGetInt( "numFilt", &numFilt );/* number of inlines in filter default=5 */
    uParGetInt( "overlapF", &overlapF );
    if(numDesiF<9) numDesiF=9;
    if(numDesiF>61) numDesiF=61;
    numDesiF=numDesiF/2*2+1;
    if(numFilt<3) numFilt=3;
    if(numFilt>7) numFilt=7;
    numFilt=numFilt/2*2+1;
    if(overlapF<numFilt/2+1) overlapF=numFilt/2+1;
    if(overlapF<3) overlapF=3;
    if((numDesiF-3)/2<overlapF) numDesiF=overlapF*2+3;
  }
  hiradon=1.;
  uParGetInt( "radonF", &radonF );       /* doing inline dip filter? yes or no */
  if(radonF>0) 
  {
    uParGetFloat( "hiradon", &hiradon );    /* higher radon transform accurace */
    if(hiradon<1.) hiradon=1.;
	if(radonF==3) uParGetInt( "ixSame", &ixSame);
    uParGetInt( "method", &method );
	if(radonF!=2)
	{
      uParGetInt( "numDesiI", &numDesiI );/* number of lines in filter design default=25 */
	  uParGetInt( "dipI", &dipI ); /* dip inline direction (ms/trace) */
	  uParGetString( "resoluI", &sresoluI );
	  uParGetString( "stDipI", &sstDipI );
   vclrf(resoluI,6);
   i=0;
   token=strtok(sresoluI,sep0);
   while(token!=NULL)
   {
     resoluI[i]   = ( float )atof( token );
	 i++;
   	 token=strtok(NULL,sep0);
   }
   vclrf(stDipI,6);
   i=0;
   token=strtok(sstDipI,sep0);
   while(token!=NULL)
   {
     stDipI[i]   = ( float )atof( token );
	 i++;
   	 token=strtok(NULL,sep0);
   }
/*      resoluI[0]   = ( float )atof( strtok( sresoluI, "-" ) );
      for(i=1;i<6;i++) resoluI[i] = ( float )atof( strtok( NULL, "-" ) );
      stDipI[0]   = ( float )atof( strtok( sstDipI, "-" ) );
      for(i=1;i<6;i++) stDipI[i] = ( float )atof( strtok( NULL, "-" ) );*/

      uParGetInt( "overlapI", &overlapI );
	  if(radonF==4) overlapI=3;
/*      if(numDesiI<9) numDesiI=9;*/
      if(numDesiI>61) numDesiI=61;
      numDesiI=numDesiI/2*2+1;
      if(dipI<4) dipI=4;
      if(overlapI<3) overlapI=3;
      if(overlapI*2+1>numDesiI) overlapI=(numDesiI-1)/2;
	}
	if(radonF==2 || (radonF==3 && ixSame==0))
	{
      uParGetInt( "numDesiX", &numDesiX );/* number of lines in filter design default=25 */
	  uParGetInt( "dipX", &dipX ); /* dip inline direction (ms/trace) */
	  uParGetString( "resoluX", &sresoluX );
	  uParGetString( "stDipX", &sstDipX );
   vclrf(resoluX,6);
   i=0;
   token=strtok(sresoluX,sep0);
   while(token!=NULL)
   {
     resoluX[i]   = ( float )atof( token );
	 i++;
   	 token=strtok(NULL,sep0);
   }
   vclrf(stDipX,6);
   i=0;
   token=strtok(sstDipX,sep0);
   while(token!=NULL)
   {
     stDipX[i]   = ( float )atof( token );
	 i++;
   	 token=strtok(NULL,sep0);
   }
/*      resoluX[0]   = ( float )atof( strtok( sresoluX, "-" ) );
      for(i=1;i<6;i++) resoluX[i] = ( float )atof( strtok( NULL, "-" ) );
      stDipX[0]   = ( float )atof( strtok( sstDipX, "-" ) );
      for(i=1;i<6;i++) stDipX[i] = ( float )atof( strtok( NULL, "-" ) );*/
      uParGetInt( "overlapX", &overlapX );
/*        if(numDesiX<9) numDesiX=9;*/
        if(numDesiX>61) numDesiX=61;
        numDesiX=numDesiX/2*2+1;
        if(dipX<6) dipX=6;
        if(overlapX<3) overlapX=3;
        if(overlapX*2+1>numDesiX) overlapX=(numDesiX-1)/2;
	}
	  if(method==1 || radonF==4 || (radonF==3 && ixSame==1))
	  {
        numDesiX=numDesiI;
		dipX=dipI;
		vmovf(resoluI,resoluX,numGroup);
		vmovf(stDipI,stDipX,numGroup);
		overlapX=overlapI;
	  }
  }

  uParGetInt( "AmpBa", &AmpBa );      
  if(AmpBa) uParGetFloat( "AmpParm", &AmpParm ); 
  else AmpParm=0.;
  uParGetInt( "FreBa", &FreBa );      
  if(FreBa>0) uParGetFloat( "FreParm", &FreParm ); 
  else FreParm=0.;
  uParGetInt( "OuPut", &OuPut );      
  nbin=1;
  if(OuPut==0) uParGetInt( "nbin", &nbin );
  if(idtyp) nbin=1;

  numDesi=numDesiF;
  if(numDesiF<numDesiI) numDesi=numDesiI;
  if(numDesi<numDesiX) numDesi=numDesiX;
   nlines = lmax-lmin+1;
   icrb = xmax-xmin+1;
   ncrbi = maxXlin-minXlin+1;
   idmap = icrb*nlines;

	  if(radonF==1)
	  {
        if(nlines<numDesiI) numDesiX=9;
		else numDesiX=numDesiI;
		vmovf(resoluI,resoluX,numGroup);
		vmovf(stDipI,stDipX,numGroup);
        dipX=dipI;
        overlapX=3;
	  }
	  if(radonF==2)
	  {
        if(icrb<numDesiX) numDesiI=9;
		else numDesiI=numDesiX;
		vmovf(resoluX,resoluI,numGroup);
		vmovf(stDipX,stDipI,numGroup);
        dipI=dipX;
        overlapI=3;
	  }

  if(!idtyp) 
  {
    i = stHdrExists( "pstmBin" );
    if( i != 1 )
     indexpstmBin = stHdrAdd( "pstmBin", "Offset bin number used in PSTM", 1, HDRINT );
    else 
     indexpstmBin=stHdrIndex( "pstmBin" );
    indexSeqno  = stHdrIndex( "seqno" );
    indexEnd_ens = stHdrIndex( "end_ens" );                /* end os ensamples flag */
    globalRuntime->itrno_valid=0;           /* trace no. novalid after this process */
	if(OuPut==1)
	{
      globalRuntime->ipkey=indexpstmBin+1;
      globalRuntime->iskey=indexCdp+1;
      globalRuntime->maxdtr=idmap;  /* making more trace ouput in each CDP */
	}
	if(OuPut==0)
	{
      globalRuntime->ipkey=indexCdp+1;
      globalRuntime->iskey=indexpstmBin+1;
      globalRuntime->maxdtr=nbin;  /* making more trace ouput in each CDP */
	}
  }
	  
   printf( "\n Linux Socket VHF 6.0 job parameters\n" );
   printf( "--------------\n\n" );
   if(dbaseParm1) printf( "the database available \n" );
   else printf( "the database is not available \n" );
   printf( "minimum inline no. to fxy  %7d\n", lmin );
   printf( "maximum inline no. to fxy  %7d\n", lmax );
   printf( "minimum xline no. to fxy  %7d\n", xmin );
   printf( "maximum xline no. to fxy  %7d\n", xmax );
   printf( "length of operator time window %7d\n", tleng );
   printf( "time window taper %7d\n", taper );
   printf( "total put values = %7d in each sequence \n",numGroup);
   printf( "time sequence \n");
   for(i=0;i<numGroup;i++) printf( " %7d ", timeseq[i] );
   printf( " \n");
   printf( "qmin minimum frequency sequence \n");
   for(i=0;i<numGroup;i++) printf( " %6.1f ", qmin[i] );
   printf( " \n");
   printf( "larger than minimum frequency sequence for taper \n");
   for(i=0;i<numGroup;i++) printf( " %6.1f ", qminp[i] );
   printf( " \n");
   printf( "maximum frequency sequence \n");
   for(i=0;i<numGroup;i++) printf( " %6.1f ", qmax[i] );
   printf( " \n");
   printf( "less than maximum frequency sequence for taper \n");
   for(i=0;i<numGroup;i++) printf( " %6.1f ", qmaxm[i] );
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
	 if(radonF==3) 
	 {
	   printf( "you use two pass \n");
	   if(ixSame) printf( "the parms are the same in both inline and xline \n");
	   else printf( "the parms are the different in inline and xline \n");
	 }
	 if(radonF!=2)
	 {
       printf( "the following parms using inline direction \n");
       printf( "number of lines for radon filter design %7d\n", numDesiI );
	   printf( "dip(ms/trace) for line radon filter used %7d\n", dipI );
       printf( "overlap between Radon group = %5d \n",overlapI);
       printf( "increasing resolutoin coefficient sequence for inline direction \n");
       for(i=0;i<numGroup;i++) printf( " %6.2f ", resoluI[i] );
       printf( " \n");
       printf( "enhancing stronger dips sequence for inline direction \n");
       for(i=0;i<numGroup;i++) printf( " %6.2f ", stDipI[i] );
       printf( " \n");
	 }
	 if(radonF==2 || (radonF==3 && ixSame==0))
	 {
       printf( "the following parms using xline direction \n");
       printf( "number of lines for radon filter design %7d\n", numDesiX );
	   printf( "dip(ms/trace) for line radon filter used %7d\n", dipX );
       printf( "overlap between Radon group = %5d \n",overlapX);
       printf( "increasing resolutoin coefficient sequence for xline direction \n");
       for(i=0;i<numGroup;i++) printf( " %6.2f ", resoluX[i] );
       printf( " \n");
       printf( "enhancing stronger dips sequence for xline direction \n");
       for(i=0;i<numGroup;i++) printf( " %6.2f ", stDipX[i] );
       printf( " \n");
	 }
   }
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
   if(AmpBa) printf( "amplitude balacing performing uisng parm %8.2f \n",AmpParm);
   else printf( "amplitude balacing not performing \n");
   if(FreBa==1) printf( "low frequency balacing performing, the parm = %8.2f \n",FreParm);
   else printf( "low frequency balacing not performing \n");

   printf( "\n Linux Socket VHF 6.0 geometry parameters\n" );
   printf( "--------------\n\n" );
   printf( "first CDP               %7ld\n", minCDP );
   printf( "last CDP                %7ld\n", maxCDP );
   printf( "max cdp fold            %7d\n",  ntrcdp );
   printf( "first inline            %7ld\n", minIlin );
   printf( "last inline             %7ld\n", maxIlin );
   printf( "first xline             %7ld\n", minXlin );
   printf( "last xline              %7ld\n", maxXlin ); 
   printf( "sampleRate              %8.2f\n",sampleRate );
   printf( "samples                 %7d\n",  samples );
   printf( "headers                 %7d\n",  headers );

   printf("you process total inline number %7d \n",nlines);
   printf("you process total crossline number %7d \n",icrb);
   printf("you process total CRB number %7d \n",idmap); 

/*      uErrStop("fxy stop"); */                     
 }

void  PreFxy( void )
{

/* alloc arry */
  shtotal=samples+headers+1;
  shtotal1=shtotal*10;
  shtotal4=shtotal*40;
  traceIn = alloc1float( shtotal1 );
  trace = alloc1float( samples+10 );
  vclrf(trace,samples+10);
  ithdr = alloc1int( headers+3 );
  rthdr = ( float * )ithdr;

}

void  PutGlobals( void )
{
  int i;

/*  fprintf(fp,"overlapI,numDesiI,dipI,resoluI,strongDipI %7d %7d %7d %7.2f %7.2f \n",overlapI,numDesiI,dipI,resoluI,stDipI);
  fprintf(fp,"overlapX,numDesiX,dipX,resoluX,strongDipX %7d %7d %7d %7.2f %7.2f \n",overlapX,numDesiX,dipX,resoluX,stDipX);
  fflush(fp);*/

/*   fprintf(fp,"put global \n");
   fflush(fp);*/
   write( clientFd, ( char* )&sampleRate, 4 );
   write( clientFd, ( char* )&samples, 4 );
   write( clientFd, ( char* )&headers, 4 );
   write( clientFd, ( char* )&shtotal, 4 );
   write( clientFd, ( char* )&minCDP, 4 );
   write( clientFd, ( char* )&minIlin, 4 );
   write( clientFd, ( char* )&minXlin, 4 );
   write( clientFd, ( char* )&lmin, 4 );
   write( clientFd, ( char* )&xmin, 4 );

   write( clientFd, ( char* )&idtyp, 4 );
   write( clientFd, ( char* )&icrb, 4 );
   write( clientFd, ( char* )&nlines, 4 );
   write( clientFd, ( char* )&idmap, 4 );
   write( clientFd, ( char* )&ncrbi, 4 );

   write( clientFd, ( char* )&tleng, 4 );
   write( clientFd, ( char* )&taper, 4 );
   write( clientFd, ( char* )&fxyF, 4 );
   write( clientFd, ( char* )&numDesiF, 4 );
   write( clientFd, ( char* )&numFilt, 4 );
   write( clientFd, ( char* )&overlapF, 4 );

   write( clientFd, ( char* )&radonF, 4 );
   write( clientFd, ( char* )&hiradon, 4 );
   write( clientFd, ( char* )&method, 4 );
   write( clientFd, ( char* )&numDesiI, 4 );
   write( clientFd, ( char* )&dipI, 4 );
   write( clientFd, ( char* )&overlapI, 4 );
   write( clientFd, ( char* )&numDesiX, 4 );
   write( clientFd, ( char* )&dipX, 4 );
   write( clientFd, ( char* )&overlapX, 4 );

   write( clientFd, ( char* )&numGroup, 4 );
   write( clientFd, ( char* )&timeseq[0], 4*numGroup );
   write( clientFd, ( char* )&qmin[0], 4*numGroup );
   write( clientFd, ( char* )&qminp[0], 4*numGroup );
   write( clientFd, ( char* )&qmax[0], 4*numGroup );
   write( clientFd, ( char* )&qmaxm[0], 4*numGroup );
   write( clientFd, ( char* )&stDipI[0], 4*numGroup  );
   write( clientFd, ( char* )&resoluI[0], 4*numGroup  );
   write( clientFd, ( char* )&stDipX[0], 4*numGroup  );
   write( clientFd, ( char* )&resoluX[0], 4*numGroup  );

   write( clientFd, ( char* )&AmpBa, 4 );
   write( clientFd, ( char* )&AmpParm, 4 );
   write( clientFd, ( char* )&FreBa, 4 );
   write( clientFd, ( char* )&FreParm, 4 );
   write( clientFd, ( char* )&OuPut, 4 );
   write( clientFd, ( char* )&nbin, 4 );

   write( clientFd, ( char* )&indexCdp, 4 );
   write( clientFd, ( char* )&indexTrc_type, 4 );
   if(!idtyp)
   {
     write( clientFd, ( char* )&indexSeqno, 4 );
     write( clientFd, ( char* )&indexEnd_ens, 4 );
     write( clientFd, ( char* )&indexpstmBin, 4 );
   }
}

/* binStart */
void binStart()
{

/* fprintf(fp,"binStart numline,lmin,lmax %6d %6d %6d \n",numline,lmin,lmax);
       fflush(fp);*/
  if(!idtyp) bin=ithdr[indexpstmBin];
  else bin=0;
  traceCount=0;

/* fprintf(fp,"binStart,bin %7d \n",bin);
       fflush(fp);*/
}

void  Fxy()
{
  int i,k,cdp,il,xl;

  if(idtyp) binStart();
  while( 1 )                                 
  {
    if( ! stGetTrace( trace, ithdr ) ) 
	{
	  if(!idtyp) 
	  {
		traceIn[0]=-3.;
/* fprintf(fp,"traceIn[0]==-3 %7.1f \n",traceIn[0]);
       fflush(fp);*/
	    write( clientFd, ( char* )traceIn, shtotal4 );
	  }
	  if(OuPut==0 || idtyp) Perform();
	  break;
	}
    else
    {      
      cdp=ithdr[indexCdp]-minCDP;                   /* no of cdp */
      il=cdp/ncrbi+minIlin;
      xl=cdp%ncrbi+minXlin;
      if(ithdr[indexSeqno]==1 && !idtyp) binStart();
	  if(ithdr[indexEnd_ens]==1 && !idtyp) 
	  {
		if(il<lmin || il>lmax || xl>xmax || xl<xmin) k=0;
		else 
		{
          il=il-lmin;
          xl=xl-xmin;  
          cdp=il*icrb+xl;
		  k=-1;
		}
		traceDisk(k);
		if(OuPut==1 && !idtyp) Perform();
		continue;
      }
	  if(il<lmin || il>lmax || xl>xmax || xl<xmin) continue;
      il=il-lmin;
      xl=xl-xmin;  
      cdp=il*icrb+xl;
/*   fprintf(fp,"fxy cdp,il,xl %7d %7d %7d  \n",cdp,il,xl);
   fflush(fp);*/
	  traceDisk(1);               /* traces input to disk */
    }
  }
}

void traceDisk(int kk)
{
  int i,j,k,k1,ii;
  static int first=0; 
  static float norm; 

  k=traceCount*shtotal;
  if(kk!=0)
  {
	vmovf(rthdr,&traceIn[k+1],headers);
	vmovf(trace,&traceIn[k+headers+1],samples);
    traceIn[k]=(float)kk;
	traceCount++;
  }
  if(traceCount==10 || kk<1)
  {
	if(kk==0 && traceCount>0) traceIn[k]=-1.;  
	if(traceCount==0) traceIn[0]=-2.;
/*   fprintf(fp,"traceDisk traceCount,kk %7d %7d  \n",traceCount,kk);
   fflush(fp);*/
	write( clientFd, ( char* )traceIn, shtotal4 );
	traceCount=0;
  }
}

void Perform()
{
  int i,j,k,nn;

  if(idtyp)
  {
	if(traceCount==0) traceIn[0]=-2.;
	else traceIn[(traceCount-1)*shtotal]=-1.;
	write( clientFd, ( char* )traceIn, shtotal4 );
  }
/*  nn=0;
   fprintf(fp,"nn %5d \n",nn);
   fflush(fp);*/
  while(1)
  {
	vclrf(traceIn,shtotal1);
    GetSocketData( ( char* )traceIn, shtotal4);
/*	nn++;*/
	if(traceIn[0]==-2.) break;
	for(i=0;i<10;i++)
	{
      k=i*shtotal;
	  vmovf(&traceIn[k+1],rthdr,headers);
      stPutTrace( &traceIn[k+headers+1], ithdr);
/*   fprintf(fp,"traceIn[k],indexCdp,indexseqno %7.1f %5d %5d \n",traceIn[k],ithdr[indexCdp],ithdr[indexSeqno]);
   fflush(fp);*/
	  if(traceIn[k]==-1.) break;
	}
/*   fprintf(fp,"nn,i %5d %5d \n",nn,i);
   fflush(fp);*/
	if(traceIn[k]==-1. && i<10) break;
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
     
void  GetSocketData( char *theData, long length )
{
   long  bytesRead = 0;
   long  n;

   while( length > 0 )
      {
      n = read( clientFd, theData+bytesRead, length );
      length -= n;
      bytesRead += n;
      }
}

void  ConnectRemote( void )
{
    signal( SIGCHLD, SIG_IGN );

    if( ( serverFd = socket( AF_INET, SOCK_STREAM, 0 ) ) < 0 )
            {
            perror( "generate error" );
            exit( 1 );
            }
    printf( "socket created\n" );

    memset( &serverAddr, 0, sizeof ( serverAddr ) );
    serverAddr.sin_family = AF_INET;
    serverAddr.sin_addr.s_addr = htonl( INADDR_ANY );
    serverAddr.sin_port = htons( sockPort );

    if( bind( serverFd, ( struct sockaddr * )&serverAddr, sizeof( serverAddr ) ) < 0 )
            {
            perror( "bind error" );
            close( serverFd );
            exit( 2 );
            }
    printf( "bind complete\n" );

    if( listen( serverFd, 1 ) < 0 )
            {
            perror( "listen error" );
            exit( 3 );
            }
    printf( "listen complete\n" );

    clientLen = sizeof( clientAddr );
    if( ( clientFd = accept( serverFd, ( struct sockaddr * )&clientAddr, &clientLen ) ) < 0 )
            {
            perror( "accept error" );
            close( serverFd );
            exit( 4 );
            }
}

