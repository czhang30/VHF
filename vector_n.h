#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <time.h>
#include  <unistd.h>

#include  <signal.h>
#include  <sys/types.h>
#include  <sys/socket.h>
#include  <netdb.h>
#include  <netinet/in.h>
#include  <arpa/inet.h>
#include  <fcntl.h>
#include  <sys/errno.h>
#include  <sys/ipc.h>
#include  <sys/shm.h>
#include  <errno.h>

/* include ProMAX prototypes and globals */

void pfacc0( long isign, long n, float *z );
long rpfN( long nmin );
long cpfN( long nmin );
long RealPrimeFactorN( long nmin, long nmax );
long ComplexPrimeFactorN( long nmin, long nmax );
void ForwardPrimeFFT( float *z, long n );
void ReversePrimeFFT( float *z, long n );
void pfacc0( long isign, long n, float *z );
float  **FloatMatrix( int nrl, int nrh, int ncl, int nch );
int  Power2( int samples );
float  *FloatVector( int nl, int nh );
int  *IntVector( int nl, int nh );

