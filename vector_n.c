#include  "vector_n.h"

#define NTAB 240
#define NFAX 10
#define  TINY 1.0e-20

#define PI  ( 3.14159265358979323846 )

#define P120 0.120536680
#define P142 0.142314838
#define P173 0.173648178
#define P222 0.222520934
#define P239 0.239315664
#define P281 0.281732557
#define P342 0.342020143
#define P354 0.354604887
#define P382 0.382683432
#define P415 0.415415013
#define P433 0.433883739
#define P464 0.464723172
#define P540 0.540640817
#define P559 0.559016994
#define P568 0.568064747
#define P587 0.587785252
#define P623 0.623489802
#define P642 0.642787610
#define P654 0.654860734
#define P663 0.663122658
#define P707 0.707106781
#define P748 0.748510748
#define P755 0.755749574
#define P766 0.766044443
#define P781 0.781831482
#define P822 0.822983866
#define P841 0.841253533
#define P866 0.866025404
#define P885 0.885456026
#define P900 0.900968868
#define P909 0.909631995
#define P923 0.923879533
#define P935 0.935016243
#define P939 0.939692621
#define P951 0.951056516
#define P959 0.959492974
#define P970 0.970941817
#define P974 0.974927912
#define P984 0.984807753
#define P989 0.989821442
#define P992 0.992708874

static struct
	{
	long n;
	float c;
	} nctab[NTAB] = {
{      1, 0.00000447 },
{      2, 0.00000524 },
{      3, 0.00000542 },
{      4, 0.00000555 },
{      5, 0.00000578 },
{      6, 0.00000734 },
{      7, 0.00000623 },
{      8, 0.00000621 },
{      9, 0.00000675 },
{     10, 0.00000849 },
{     11, 0.00000894 },
{     12, 0.00000912 },
{     13, 0.00001278 },
{     14, 0.00000992 },
{     15, 0.00001006 },
{     16, 0.00000927 },
{     18, 0.00001150 },
{     20, 0.00001142 },
{     21, 0.00001228 },
{     22, 0.00001644 },
{     24, 0.00001217 },
{     26, 0.00002484 },
{     28, 0.00001373 },
{     30, 0.00001525 },
{     33, 0.00002122 },
{     35, 0.00001653 },
{     36, 0.00001700 },
{     39, 0.00003445 },
{     40, 0.00001728 },
{     42, 0.00002135 },
{     44, 0.00002561 },
{     45, 0.00001975 },
{     48, 0.00002533 },
{     52, 0.00004304 },
{     55, 0.00003292 },
{     56, 0.00002232 },
{     60, 0.00002683 },
{     63, 0.00002660 },
{     65, 0.00005281 },
{     66, 0.00003933 },
{     70, 0.00003179 },
{     72, 0.00002810 },
{     77, 0.00004385 },
{     78, 0.00006717 },
{     80, 0.00003822 },
{     84, 0.00003430 },
{     88, 0.00004730 },
{     90, 0.00003738 },
{     91, 0.00007141 },
{     99, 0.00005596 },
{    104, 0.00007901 },
{    105, 0.00004317 },
{    110, 0.00006585 },
{    112, 0.00005040 },
{    117, 0.00009263 },
{    120, 0.00004883 },
{    126, 0.00005233 },
{    130, 0.00010454 },
{    132, 0.00007352 },
{    140, 0.00005697 },
{    143, 0.00013565 },
{    144, 0.00006560 },
{    154, 0.00008791 },
{    156, 0.00012719 },
{    165, 0.00009488 },
{    168, 0.00006743 },
{    176, 0.00010731 },
{    180, 0.00007250 },
{    182, 0.00014611 },
{    195, 0.00015979 },
{    198, 0.00011825 },
{    208, 0.00017536 },
{    210, 0.00009567 },
{    220, 0.00012467 },
{    231, 0.00013668 },
{    234, 0.00019435 },
{    240, 0.00011833 },
{    252, 0.00010185 },
{    260, 0.00020584 },
{    264, 0.00014740 },
{    273, 0.00022219 },
{    280, 0.00010928 },
{    286, 0.00027964 },
{    308, 0.00017497 },
{    312, 0.00024961 },
{    315, 0.00013650 },
{    330, 0.00020718 },
{    336, 0.00016660 },
{    360, 0.00014801 },
{    364, 0.00029525 },
{    385, 0.00022833 },
{    390, 0.00033693 },
{    396, 0.00023321 },
{    420, 0.00018959 },
{    429, 0.00042543 },
{    440, 0.00025178 },
{    455, 0.00037477 },
{    462, 0.00029453 },
{    468, 0.00038417 },
{    495, 0.00029701 },
{    504, 0.00021141 },
{    520, 0.00041746 },
{    528, 0.00034909 },
{    546, 0.00047775 },
{    560, 0.00028079 },
{    572, 0.00056647 },
{    585, 0.00049646 },
{    616, 0.00035592 },
{    624, 0.00056422 },
{    630, 0.00030628 },
{    660, 0.00041804 },
{    693, 0.00042831 },
{    715, 0.00072498 },
{    720, 0.00037201 },
{    728, 0.00059863 },
{    770, 0.00050051 },
{    780, 0.00068036 },
{    792, 0.00047195 },
{    819, 0.00070298 },
{    840, 0.00038971 },
{    858, 0.00090686 },
{    880, 0.00059890 },
{    910, 0.00081273 },
{    924, 0.00060065 },
{    936, 0.00078653 },
{    990, 0.00066282 },
{   1001, 0.00103437 },
{   1008, 0.00053763 },
{   1040, 0.00096204 },
{   1092, 0.00097834 },
{   1144, 0.00116165 },
{   1155, 0.00077654 },
{   1170, 0.00107579 },
{   1232, 0.00086242 },
{   1260, 0.00063187 },
{   1287, 0.00136386 },
{   1320, 0.00087092 },
{   1365, 0.00126090 },
{   1386, 0.00096058 },
{   1430, 0.00157132 },
{   1456, 0.00138956 },
{   1540, 0.00104381 },
{   1560, 0.00142789 },
{   1584, 0.00116187 },
{   1638, 0.00156293 },
{   1680, 0.00101045 },
{   1716, 0.00190932 },
{   1820, 0.00172937 },
{   1848, 0.00129647 },
{   1872, 0.00188508 },
{   1980, 0.00143839 },
{   2002, 0.00233567 },
{   2145, 0.00256525 },
{   2184, 0.00220544 },
{   2288, 0.00286041 },
{   2310, 0.00199268 },
{   2340, 0.00249935 },
{   2520, 0.00165966 },
{   2574, 0.00328185 },
{   2640, 0.00250880 },
{   2730, 0.00324685 },
{   2772, 0.00246072 },
{   2860, 0.00376331 },
{   3003, 0.00389973 },
{   3080, 0.00282769 },
{   3120, 0.00386549 },
{   3276, 0.00389132 },
{   3432, 0.00463863 },
{   3465, 0.00322068 },
{   3640, 0.00430320 },
{   3696, 0.00385603 },
{   3960, 0.00387239 },
{   4004, 0.00550178 },
{   4095, 0.00488055 },
{   4290, 0.00634835 },
{   4368, 0.00572348 },
{   4620, 0.00506549 },
{   4680, 0.00579823 },
{   5005, 0.00688648 },
{   5040, 0.00475490 },
{   5148, 0.00725322 },
{   5460, 0.00739223 },
{   5544, 0.00572243 },
{   5720, 0.00813454 },
{   6006, 0.00915082 },
{   6160, 0.00688741 },
{   6435, 0.00899499 },
{   6552, 0.00841530 },
{   6864, 0.01047864 },
{   6930, 0.00797688 },
{   7280, 0.00994134 },
{   7920, 0.00910392 },
{   8008, 0.01164887 },
{   8190, 0.01152095 },
{   8580, 0.01341202 },
{   9009, 0.01296296 },
{   9240, 0.01101695 },
{   9360, 0.01303952 },
{  10010, 0.01572065 },
{  10296, 0.01526346 },
{  10920, 0.01571038 },
{  11088, 0.01310844 },
{  11440, 0.01795483 },
{  12012, 0.01923848 },
{  12870, 0.02061874 },
{  13104, 0.01874544 },
{  13860, 0.01699460 },
{  15015, 0.02351713 },
{  16016, 0.02560160 },
{  16380, 0.02413479 },
{  17160, 0.02824737 },
{  18018, 0.02972973 },
{  18480, 0.02492284 },
{  20020, 0.03319398 },
{  20592, 0.03379152 },
{  21840, 0.03561436 },
{  24024, 0.04072959 },
{  25740, 0.04356223 },
{  27720, 0.03649691 },
{  30030, 0.05422948 },
{  32760, 0.05186703 },
{  34320, 0.06369732 },
{  36036, 0.06516064 },
{  40040, 0.07561521 },
{  45045, 0.08809524 },
{  48048, 0.09973118 },
{  51480, 0.10560345 },
{  55440, 0.09884259 },
{  60060, 0.13863636 },
{  65520, 0.14276557 },
{  72072, 0.16074297 },
{  80080, 0.18254505 },
{  90090, 0.23093434 },
{ 102960, 0.24238506 },
{ 120120, 0.32568027 },
{ 144144, 0.34959350 },
{ 180180, 0.54267677 },
{ 240240, 0.68750000 },
{ 360360, 1.15729167 },
{ 720720, 2.39687500 },
};

long  cpfN( long nmin )
{
long i;

	for( i = 0; i < NTAB-1 && nctab[i].n < nmin; ++i ) ;
	return( nctab[i].n );
}

long ComplexPrimeFactorN( long nmin, long nmax )
{
long  i, j;

	for( i = 0; i < NTAB-1 && nctab[i].n < nmin; ++i ) ;
	for( j = i+1; j < NTAB-1 && nctab[j].n <= nmax; ++j )
		{
		if( nctab[j].c < nctab[i].c )
			i = j;
		}
	return( nctab[i].n );
}

long rpfN( long nmin )
{
    return( 2 * cpfN( ( nmin + 1 ) / 2 ) );
}

long RealPrimeFactorN( long nmin, long nmax )
{
	return( 2 * ComplexPrimeFactorN( ( nmin + 1 ) / 2, ( nmax + 1 ) / 2 ) );
}

void  ReversePrimeFFT( float *z, long  n )
{
	double temp,theta,wtemp,wpr,wpi,wr,wi;
	long i,no2,ir,ii,jr,ji;
	double sumr,sumi,difr,difi,tempr,tempi;

	z[1] = z[0]-z[n];
	z[0] = z[0]+z[n];
    temp = 1 / ( double )n;
	for(   i = 0; i < n; i++ )
		z[i] *= temp;


	  theta = -2.0 * PI / ( double )n;
	  wtemp = sin( 0.5 * theta );
	  wpr = -2.0 * wtemp * wtemp;
	  wpi = sin( theta );
	  wr = 1.0 + wpr;
	  wi = wpi;

	  no2 = n / 2;
	for( ir = 2, ii = 3, jr = n-2, ji = n-1; ir <= no2; ir+=2, ii+=2, jr-=2, ji-=2 )
		{
		  sumr = z[ir] + z[jr];
		  sumi = z[ii] + z[ji];
		  difr = z[ir] - z[jr];
		  difi = z[ii] - z[ji];
		  tempr = wi * difr - wr * sumi;
		  tempi = wi * sumi + wr * difr;
		z[ir] = sumr + tempr;
		z[ii] = difi + tempi;
		z[jr] = sumr - tempr;
		z[ji] = tempi - difi;
		wtemp = wr;
		wr += wr * wpr - wi * wpi;
		wi += wi * wpr + wtemp * wpi;
		}


	pfacc0( 1, n/2, z);
}

void  ForwardPrimeFFT( float *z, long n )
{
	double temp,theta,wtemp,wpr,wpi,wr,wi;
	long i,no2,ir,ii,jr,ji;
	double sumr,sumi,difr,difi,tempr,tempi;

	for(  i = 0; i < n; i++) z[i] *= 0.5;

	pfacc0( -1, n/2, z);

	z[n] = 2.0 * ( z[0] - z[1] );
	z[0] = 2.0 * ( z[0] + z[1] );
	z[n+1] = 0.0;
	z[1] = 0.0;


	  theta = -2.0 * PI / (double)n;
	  wtemp = sin( 0.5 * theta );
	  wpr = -2.0 * wtemp * wtemp;
	  wpi = sin( theta );
	  wr = 1.0 + wpr;
	  wi = wpi;

     no2 = n / 2;
	for (  ir = 2, ii = 3, jr = n-2, ji = n-1; ir <= no2; ir+=2, ii+=2, jr-=2, ji-=2 )
		{
		  sumr = z[ir] + z[jr];
		  sumi = z[ii] + z[ji];
		  difr = z[ir] - z[jr];
		  difi = z[ii] - z[ji];
		  tempr = wi * difr + wr * sumi;
		  tempi = wi * sumi - wr * difr;
		z[ir] = sumr + tempr;
		z[ii] = difi + tempi;
		z[jr] = sumr - tempr;
		z[ji] = tempi - difi;
		wtemp = wr;
		wr += wr * wpr - wi * wpi;
		wi += wi * wpr + wtemp * wpi;
		}
}

void  pfacc0( long isign, long n, float *z )
{
static long  kfax[] = { 16, 13, 11, 9, 8, 7, 5, 4, 3, 2 };

/* keep track of n left after dividing by factors */
	long  nleft = n;
	long  mm = 0;
	long  mu = 0;
    long jfax,ifac,ndiv,m,jinc,jmax,j0,j1,l,jt,j2,j3,j4,j5,j6,j7,j8;
    long j9,j10,j11,j12,j13,j14,j15,jfac;
	double t1r,t1i,c1,y1r,y1i,y2r,y2i,t2r,t2i,y3r,y3i,c2,c3;
	double t3r,t3i,y4r,y4i,t4r,t4i,t5r,t5i,t6r,t6i,t7r,t7i;
	double t8r,t8i,t9r,t9i,y5r,y5i,y6r,y6i;
	double t10r,t10i,t11r,t11i,t12r,t12i,y7r,y7i;
	double t13r,t13i,t14r,t14i,t15r,t15i,y8r,y8i;
	double t16r,t16i,t17r,t17i,t18r,t18i;
	double t19r,t19i,t20r,t20i,t21r,t21i;
	double t22r,t22i,t23r,t23i;
	double t24r,t24i,t25r,t25i,t26r,t26i;
	double t27r,t27i,t28r,t28i,t29r,t29i;
	double t30r,t30i,t31r,t31i;
	double t32r,t32i,t33r,t33i;
	double t34r,t34i,t35r,t35i,t36r,t36i;
	double t37r,t37i,t38r,t38i,t39r,t39i;
	double t40r,t40i,t41r,t41i,t42r,t42i;
	double y9r,y9i,y10r,y10i;
	double y11r,y11i,y12r,y12i;
	double y13r,y13i,y14r,y14i;
	double y15r,y15i;
			double  c4, c5, c6,c7,c8,c9,c10;
			double  c11,c12;
/* begin loop over possible factors (from biggest to smallest) */
	for(  jfax = 0; jfax < NFAX; jfax++ )
		{
	     ifac = kfax[jfax];
		  ndiv = nleft / ifac;
		if( ndiv * ifac != nleft)
			continue;								
 
		nleft = ndiv;
		  m = n / ifac;
 
		for(   jfac = 1; jfac <= ifac; jfac++ )
			{
			mu = jfac;
			mm = jfac * m;
			if( mm % ifac == 1 )  break;		
			}
 
		if( isign < 0 )  mu = ifac - mu;
		  jinc = 2 * mm;
		  jmax = 2 * n;
		  j0 = 0;
		  j1 = j0 + jinc;

		if( ifac == 2 )						       
			{
			for(   l = 0; l < m; l++ )
				{
				  t1r = z[j0] - z[j1];
				  t1i = z[j0+1] - z[j1+1];
				z[j0] = z[j0] + z[j1];
				z[j0+1] = z[j0+1] + z[j1+1];
				z[j1] = t1r;
				z[j1+1] = t1i;
			  jt = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}

		 j2 = j1 + jinc;
		if( j2 >= jmax )  j2 = j2 - jmax;

		if( ifac == 3 )							
			{
			if( mu == 1 )
				c1 = P866;
			else
				c1 = -P866;
			for(   l = 0; l < m; l++ )
				{
				  t1r = z[j1] + z[j2];
				  t1i = z[j1+1] + z[j2+1];
				  y1r = z[j0] - 0.5 * t1r;
				  y1i = z[j0+1] - 0.5 * t1i;
				  y2r = c1 * ( z[j1] - z[j2] );
				  y2i = c1 * ( z[j1+1] - z[j2+1] );
				z[j0] = z[j0] + t1r;
				z[j0+1] = z[j0+1] + t1i;
				z[j1] = y1r - y2i;
				z[j1+1] = y1i + y2r;
				z[j2] = y1r + y2i;
				z[j2+1] = y1i - y2r;
				  jt = j2 + 2;
				j2 = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}

		  j3 = j2 + jinc;
		if( j3 >= jmax )  j3 = j3 - jmax;

		if( ifac == 4 )							
			{
			  c1;
			if( mu == 1 )
				c1 = 1.0;
			else
				c1 = -1.0;
			for(   l = 0; l < m; l++ )
				{
				  t1r = z[j0] + z[j2];
				  t1i = z[j0+1] + z[j2+1];
				  t2r = z[j1] + z[j3];
				  t2i = z[j1+1] + z[j3+1];
				  y1r = z[j0] - z[j2];
				  y1i = z[j0+1] - z[j2+1];
				  y3r = c1 * ( z[j1] - z[j3] );
				  y3i = c1 * ( z[j1+1] - z[j3+1] );
				z[j0] = t1r + t2r;
				z[j0+1] = t1i + t2i;
				z[j1] = y1r - y3i;
				z[j1+1] = y1i + y3r;
				z[j2] = t1r - t2r;
				z[j2+1] = t1i - t2i;
				z[j3] = y1r + y3i;
				z[j3+1] = y1i - y3r;
			   jt = j3 + 2;
				j3 = j2 + 2;
				j2 = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}

		  j4 = j3 + jinc;
		if( j4 >= jmax )  j4 = j4 - jmax;

		if( ifac == 5 )							
			{
			if( mu == 1 )
				{
				c1 = P559;
				c2 = P951;
				c3 = P587;
				}
			else if( mu == 2 )
				{
				c1 = -P559;
				c2 = P587;
				c3 = -P951;
				}
			else if( mu == 3 )
				{
				c1 = -P559;
				c2 = -P587;
				c3 = P951;
				}
			else
				{ 
				c1 = P559;
				c2 = -P951;
				c3 = -P587;
				}
			for(  l = 0; l < m; l++ )
				{
				  t1r = z[j1] + z[j4];
				  t1i = z[j1+1] + z[j4+1];
				  t2r = z[j2] + z[j3];
				  t2i = z[j2+1] + z[j3+1];
				  t3r = z[j1] - z[j4];
				  t3i = z[j1+1] - z[j4+1];
				  t4r = z[j2] - z[j3];
				  t4i = z[j2+1] - z[j3+1];
				  t5r = t1r + t2r;
				  t5i = t1i + t2i;
				  t6r = c1 * ( t1r - t2r );
				  t6i = c1 * ( t1i - t2i );
				  t7r = z[j0] - 0.25 * t5r;
				  t7i = z[j0+1] - 0.25 * t5i;
				  y1r = t7r + t6r;
				  y1i = t7i + t6i;
				  y2r = t7r - t6r;
				  y2i = t7i - t6i;
				  y3r = c3 * t3r - c2 * t4r;
				  y3i = c3 * t3i - c2 * t4i;
				  y4r = c2 * t3r + c3 * t4r;
				  y4i = c2 * t3i + c3 * t4i;
				z[j0] = z[j0] + t5r;
				z[j0+1] = z[j0+1] + t5i;
				z[j1] = y1r - y4i;
				z[j1+1] = y1i + y4r;
				z[j2] = y2r - y3i;
				z[j2+1] = y2i + y3r;
				z[j3] = y2r + y3i;
				z[j3+1] = y2i - y3r;
				z[j4] = y1r +y4i;
				z[j4+1] = y1i - y4r;
				  jt = j4 + 2;
				j4 = j3 + 2;
				j3 = j2 + 2;
				j2 = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}

		  j5 = j4 + jinc;
		if( j5 >= jmax )  j5 = j5 - jmax;
		  j6 = j5 + jinc;
		if( j6 >= jmax )  j6 = j6 - jmax;

		if( ifac == 7 )
			{
			if( mu == 1)
				{
				c1 = P623;
				c2 = -P222;
				c3 = -P900;
				c4 = P781;
				c5 = P974;
				c6 = P433;
				}
			else if( mu == 2 )
				{
				c1 = -P222;
				c2 = -P900;
				c3 = P623;
				c4 = P974;
				c5 = -P433;
				c6 = -P781;
				}
			else if( mu == 3 )
				{
				c1 = -P900;
				c2 = P623;
				c3 = -P222;
				c4 = P433;
				c5 = -P781;
				c6 = P974;
				}
			else if( mu == 4 )
				{
				c1 = -P900;
				c2 = P623;
				c3 = -P222;
				c4 = -P433;
				c5 = P781;
				c6 = -P974;
				}
			else if( mu == 5 )
				{
				c1 = -P222;
				c2 = -P900;
				c3 = P623;
				c4 = -P974;
				c5 = P433;
				c6 = P781;
				}
			else
				{
				c1 = P623;
				c2 = -P222;
				c3 = -P900;
				c4 = -P781;
				c5 = -P974;
				c6 = -P433;
				}
			for(   l = 0; l < m; l++ )
				{
				  t1r = z[j1] + z[j6];
				  t1i = z[j1+1] + z[j6+1];
				  t2r = z[j2] + z[j5];
				  t2i = z[j2+1] + z[j5+1];
				  t3r = z[j3] + z[j4];
				  t3i = z[j3+1] + z[j4+1];
				  t4r = z[j1] - z[j6];
				  t4i = z[j1+1] - z[j6+1];
				  t5r = z[j2] - z[j5];
				  t5i = z[j2+1] - z[j5+1];
				  t6r = z[j3] - z[j4];
				  t6i = z[j3+1] - z[j4+1];
				  t7r = z[j0] - 0.5 * t3r;
				  t7i = z[j0+1] - 0.5 * t3i;
				  t8r = t1r - t3r;
				  t8i = t1i - t3i;
				  t9r = t2r - t3r;
				  t9i = t2i - t3i;
				  y1r = t7r + c1 * t8r + c2 * t9r;
				  y1i = t7i + c1 * t8i + c2 * t9i;
				  y2r = t7r + c2 * t8r + c3 * t9r;
				  y2i = t7i + c2 * t8i + c3 * t9i;
				  y3r = t7r + c3 * t8r + c1 * t9r;
				  y3i = t7i + c3 * t8i + c1 * t9i;
				  y4r = c6 * t4r - c4 * t5r + c5 * t6r;
				  y4i = c6 * t4i - c4 * t5i + c5 * t6i;
				  y5r = c5 * t4r - c6 * t5r - c4 * t6r;
				  y5i = c5 * t4i - c6 * t5i - c4 * t6i;
				  y6r = c4 * t4r + c5 * t5r + c6 * t6r;
				  y6i = c4 * t4i + c5 * t5i + c6 * t6i;
				z[j0] = z[j0] + t1r + t2r + t3r;
				z[j0+1] = z[j0+1] + t1i + t2i + t3i;
				z[j1] = y1r - y6i;
				z[j1+1] = y1i + y6r;
				z[j2] = y2r - y5i;
				z[j2+1] = y2i + y5r;
				z[j3] = y3r - y4i;
				z[j3+1] = y3i + y4r;
				z[j4] = y3r + y4i;
				z[j4+1] = y3i - y4r;
				z[j5] = y2r + y5i;
				z[j5+1] = y2i - y5r;
				z[j6] = y1r + y6i;
				z[j6+1] = y1i - y6r;
				  jt = j6 + 2;
				j6 = j5 + 2;
				j5 = j4 + 2;
				j4 = j3 + 2;
				j3 = j2 + 2;
				j2 = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}

		  j7 = j6 + jinc;
		if( j7 >= jmax )  j7 = j7 - jmax;

		if( ifac == 8 )								
			{
			  c1, c2;
			if( mu==1 )
				{
				c1 = 1.0;
				c2 = P707;
				}
			else if( mu == 3 )
				{
				c1 = -1.0;
				c2 = -P707;
				}
			else if( mu == 5 )
				{
				c1 = 1.0;
				c2 = -P707;
				}
			else
				{
				c1 = -1.0;
				c2 = P707;
				}
			  c3 = c1 * c2;
			for(   l = 0; l < m; l++ )
				{
				  t1r = z[j0] + z[j4];
				  t1i = z[j0+1] + z[j4+1];
				  t2r = z[j0] - z[j4];
				  t2i = z[j0+1] - z[j4+1];
				  t3r = z[j1] + z[j5];
				  t3i = z[j1+1] + z[j5+1];
				  t4r = z[j1] - z[j5];
				  t4i = z[j1+1] - z[j5+1];
				  t5r = z[j2] + z[j6];
				  t5i = z[j2+1] + z[j6+1];
				  t6r = c1 * ( z[j2] - z[j6] );
				  t6i = c1 * ( z[j2+1] - z[j6+1] );
				  t7r = z[j3] + z[j7];
				  t7i = z[j3+1] + z[j7+1];
				  t8r = z[j3] - z[j7];
				  t8i = z[j3+1] - z[j7+1];
				  t9r = t1r + t5r;
				  t9i = t1i + t5i;
				  t10r = t3r + t7r;
				  t10i = t3i + t7i;
				  t11r = c2 * ( t4r - t8r );
				  t11i = c2 * ( t4i - t8i );
				  t12r = c3 * ( t4r + t8r );
				  t12i = c3 * ( t4i + t8i );
				  y1r = t2r + t11r;
				  y1i = t2i + t11i;
				  y2r = t1r - t5r;
				  y2i = t1i - t5i;
				  y3r = t2r - t11r;
				  y3i = t2i - t11i;
				  y5r = t12r - t6r;
				  y5i = t12i - t6i;
				  y6r = c1 * ( t3r - t7r );
				  y6i = c1 * ( t3i - t7i );
				  y7r = t12r + t6r;
				  y7i = t12i + t6i;
				z[j0] = t9r + t10r;
				z[j0+1] = t9i + t10i;
				z[j1] = y1r - y7i;
				z[j1+1] = y1i + y7r;
				z[j2] = y2r - y6i;
				z[j2+1] = y2i + y6r;
				z[j3] = y3r - y5i;
				z[j3+1] = y3i + y5r;
				z[j4] = t9r - t10r;
				z[j4+1] = t9i - t10i;
				z[j5] = y3r + y5i;
				z[j5+1] = y3i - y5r;
				z[j6] = y2r + y6i;
				z[j6+1] = y2i - y6r;
				z[j7] = y1r + y7i;
				z[j7+1] = y1i - y7r;
				  jt = j7 + 2;
				j7 = j6 + 2;
				j6 = j5 + 2;
				j5 = j4 + 2;
				j4 = j3 + 2;
				j3 = j2 + 2;
				j2 = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}

		  j8 = j7 + jinc;
		if( j8 >= jmax )  j8 = j8 - jmax;

		if( ifac==9 )										
			{
			if( mu == 1)
				{
				c1 = P866;
				c2 = P766;
				c3 = P642;
				c4 = P173;
				c5 = P984;
				}
			else if( mu == 2 )
				{
				c1 = -P866;
				c2 = P173;
				c3 = P984;
				c4 = -P939;
				c5 = P342;
				}
			else if( mu == 4 )
				{
				c1 = P866;
				c2 = -P939;
				c3 = P342;
				c4 = P766;
				c5 = -P642;
				}
			else if( mu == 5 )
				{
				c1 = -P866;
				c2 = -P939;
				c3 = -P342;
				c4 = P766;
				c5 = P642;
				}
			else if( mu == 7 )
				{
				c1 = P866;
				c2 = P173;
				c3 = -P984;
				c4 = -P939;
				c5 = -P342;
				}
			else
				{
				c1 = -P866;
				c2 = P766;
				c3 = -P642;
				c4 = P173;
				c5 = -P984;
				}
			  c6 = c1 * c2;
			  c7 = c1 * c3;
			  c8 = c1 * c4;
			  c9 = c1 * c5;
			for(   l = 0; l < m; l++)
				{
				  t1r = z[j3] + z[j6];
				  t1i = z[j3+1] + z[j6+1];
				  t2r = z[j0] - 0.5 * t1r;
				  t2i = z[j0+1] - 0.5 * t1i;
				  t3r = c1 * ( z[j3] - z[j6] );
				  t3i = c1 * ( z[j3+1] - z[j6+1] );
				  t4r = z[j0] + t1r;
				  t4i = z[j0+1] + t1i;
				  t5r = z[j4] + z[j7];
				  t5i = z[j4+1] + z[j7+1];
				  t6r = z[j1] - 0.5 * t5r;
				  t6i = z[j1+1] - 0.5 * t5i;
				  t7r = z[j4] - z[j7];
				  t7i = z[j4+1] - z[j7+1];
				  t8r = z[j1] + t5r;
				  t8i = z[j1+1] + t5i;
				  t9r = z[j2] + z[j5];
				  t9i = z[j2+1] + z[j5+1];
				  t10r = z[j8] - 0.5 * t9r;
				  t10i = z[j8+1] - 0.5 * t9i;
				  t11r = z[j2] - z[j5];
				  t11i = z[j2+1] - z[j5+1];
				  t12r = z[j8] + t9r;
				  t12i = z[j8+1] + t9i;
				  t13r = t8r + t12r;
				  t13i = t8i + t12i;
				  t14r = t6r + t10r;
				  t14i = t6i + t10i;
				  t15r = t6r - t10r;
				  t15i = t6i - t10i;
				  t16r = t7r + t11r;
				  t16i = t7i + t11i;
				  t17r = t7r - t11r;
				  t17i = t7i - t11i;
				  t18r = c2 * t14r - c7 * t17r;
				  t18i = c2 * t14i - c7 * t17i;
				  t19r = c4 * t14r + c9 * t17r;
				  t19i = c4 * t14i + c9 * t17i;
				  t20r = c3 * t15r + c6 * t16r;
				  t20i = c3 * t15i + c6 * t16i;
				  t21r = c5 * t15r - c8 * t16r;
				  t21i = c5 * t15i - c8 * t16i;
				  t22r = t18r + t19r;
				  t22i = t18i + t19i;
				  t23r = t20r  -  t21r;
				  t23i = t20i - t21i;
				  y1r = t2r + t18r;
				  y1i = t2i + t18i;
				  y2r = t2r + t19r;
				  y2i = t2i + t19i;
				  y3r = t4r - 0.5 * t13r;
				  y3i = t4i - 0.5 * t13i;
				  y4r = t2r - t22r;
				  y4i = t2i - t22i;
				  y5r = t3r - t23r;
				  y5i = t3i - t23i;
				  y6r = c1 * ( t8r - t12r );
				  y6i = c1 * ( t8i - t12i );
				  y7r = t21r - t3r;
				  y7i = t21i - t3i;
				  y8r = t3r + t20r;
				  y8i = t3i + t20i;
				z[j0] = t4r + t13r;
				z[j0+1] = t4i + t13i;
				z[j1] = y1r - y8i;
				z[j1+1] = y1i + y8r;
				z[j2] = y2r - y7i;
				z[j2+1] = y2i + y7r;
				z[j3] = y3r - y6i;
				z[j3+1] = y3i + y6r;
				z[j4] = y4r - y5i;
				z[j4+1] = y4i + y5r;
				z[j5] = y4r + y5i;
				z[j5+1] = y4i - y5r;
				z[j6] = y3r + y6i;
				z[j6+1] = y3i - y6r;
				z[j7] = y2r + y7i;
				z[j7+1] = y2i - y7r;
				z[j8] = y1r + y8i;
				z[j8+1] = y1i - y8r;
				  jt = j8 + 2;
				j8 = j7 + 2;
				j7 = j6 + 2;
				j6 = j5 + 2;
				j5 = j4 + 2;
				j4 = j3 + 2;
				j3 = j2 + 2;
				j2 = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}

		  j9 = j8 + jinc;
		if( j9 >= jmax )  j9 = j9 - jmax;
		  j10 = j9 + jinc;
		if( j10 >= jmax ) j10 = j10 - jmax;

		if( ifac==11 )										
			{

			if( mu == 1 )
				{
				c1 = P841;
				c2 = P415;
				c3 = -P142;
				c4 = -P654;
				c5 = -P959;
				c6 = P540;
				c7 = P909;
				c8 = P989;
				c9 = P755;
				c10 = P281;
				}
			else if( mu == 2 )
				{
				c1 = P415;
				c2 = -P654;
				c3 = -P959;
				c4 = -P142;
				c5 = P841;
				c6 = P909;
				c7 = P755;
				c8 = -P281;
				c9 = -P989;
				c10 = -P540;
				}
			else if( mu == 3 )
				{
				c1 = -P142;
				c2 = -P959;
				c3 = P415;
				c4 = P841;
				c5 = -P654;
				c6 = P989;
				c7 = -P281;
				c8 = -P909;
				c9 = P540;
				c10 = P755;
				}
			else if( mu == 4 )
				{
				c1 = -P654;
				c2 = -P142;
				c3 = P841;
				c4 = -P959;
				c5 = P415;
				c6 = P755;
				c7 = -P989;
				c8 = P540;
				c9 = P281;
				c10 = -P909;
				}
			else if( mu == 5 )
				{
				c1 = -P959;
				c2 = P841;
				c3 = -P654;
				c4 = P415;
				c5 = -P142;
				c6 = P281;
				c7 = -P540;
				c8 = P755;
				c9 = -P909;
				c10 = P989;
				}
			else if( mu == 6 )
				{
				c1 = -P959;
				c2 = P841;
				c3 = -P654;
				c4 = P415;
				c5 = -P142;
				c6 = -P281;
				c7 = P540;
				c8 = -P755;
				c9 = P909;
				c10 = -P989;
				}
			else if( mu == 7 )
				{
				c1 = -P654;
				c2 = -P142;
				c3 = P841;
				c4 = -P959;
				c5 = P415;
				c6 = -P755;
				c7 = P989;
				c8 = -P540;
				c9 = -P281;
				c10 = P909;
				}
			else if( mu == 8 )
				{
				c1 = -P142;
				c2 = -P959;
				c3 = P415;
				c4 = P841;
				c5 = -P654;
				c6 = -P989;
				c7 = P281;
				c8 = P909;
				c9 = -P540;
				c10 = -P755;
				}
			else if( mu == 9 )
				{
				c1 = P415;
				c2 = -P654;
				c3 = -P959;
				c4 = -P142;
				c5 = P841;
				c6 = -P909;
				c7 = -P755;
				c8 = P281;
				c9 = P989;
				c10 = P540;
				}
			else
				{
				c1 = P841;
				c2 = P415;
				c3 = -P142;
				c4 = -P654;
				c5 = -P959;
				c6 = -P540;
				c7 = -P909;
				c8 = -P989;
				c9 = -P755;
				c10 = -P281;
				}
			for(   l = 0; l < m; l ++ )
				{
				  t1r = z[j1] + z[j10];
				  t1i = z[j1+1] + z[j10+1];
				  t2r = z[j2] + z[j9];
				  t2i = z[j2+1] + z[j9+1];
				  t3r = z[j3] + z[j8];
				  t3i = z[j3+1] + z[j8+1];
				  t4r = z[j4] + z[j7];
				  t4i = z[j4+1] + z[j7+1];
				  t5r = z[j5] + z[j6];
				  t5i = z[j5+1] + z[j6+1];
				  t6r = z[j1] - z[j10];
				  t6i = z[j1+1] - z[j10+1];
				  t7r = z[j2] - z[j9];
				  t7i = z[j2+1] - z[j9+1];
				  t8r = z[j3] - z[j8];
				  t8i = z[j3+1] - z[j8+1];
				  t9r = z[j4] - z[j7];
				  t9i = z[j4+1] - z[j7+1];
				  t10r = z[j5] - z[j6];
				  t10i = z[j5+1] - z[j6+1];
				  t11r = z[j0] - 0.5 * t5r;
				  t11i = z[j0+1] - 0.5 * t5i;
				  t12r = t1r - t5r;
				  t12i = t1i - t5i;
				  t13r = t2r - t5r;
				  t13i = t2i - t5i;
				  t14r = t3r - t5r;
				  t14i = t3i - t5i;
				  t15r = t4r - t5r;
				  t15i = t4i - t5i;
				  y1r = t11r + c1 * t12r + c2 * t13r + c3 * t14r + c4 * t15r;
				  y1i = t11i + c1 * t12i + c2 * t13i + c3 * t14i + c4 * t15i;
				  y2r = t11r + c2 * t12r + c4 * t13r + c5 * t14r + c3 * t15r;
				  y2i = t11i + c2 * t12i + c4 * t13i + c5 * t14i + c3 * t15i;
				  y3r = t11r + c3 * t12r + c5 * t13r + c2 * t14r + c1 * t15r;
				  y3i = t11i + c3 * t12i + c5 * t13i + c2 * t14i + c1 * t15i;
				  y4r = t11r + c4 * t12r + c3 * t13r + c1 * t14r + c5 * t15r;
				  y4i = t11i + c4 * t12i + c3 * t13i + c1 * t14i + c5 * t15i;
				  y5r = t11r + c5 * t12r + c1 * t13r + c4 * t14r + c2 * t15r;
				  y5i = t11i + c5 * t12i + c1 * t13i + c4 * t14i + c2 * t15i;
				  y6r = c10 * t6r - c6 * t7r + c9 * t8r - c7 * t9r + c8 * t10r;
				  y6i = c10 * t6i - c6 * t7i + c9 * t8i - c7 * t9i + c8 * t10i;
				  y7r = c9 * t6r - c8 * t7r + c6 * t8r + c10 * t9r - c7 * t10r;
				  y7i = c9 * t6i - c8 * t7i + c6 * t8i + c10 * t9i - c7 * t10i;
				  y8r = c8 * t6r - c10 * t7r - c7 * t8r + c6 * t9r + c9 * t10r;
				  y8i = c8 * t6i - c10 * t7i - c7 * t8i + c6 * t9i + c9 * t10i;
				  y9r = c7 * t6r + c9 * t7r - c10 * t8r - c8 * t9r - c6 * t10r;
				  y9i = c7 * t6i + c9 * t7i - c10 * t8i - c8 * t9i - c6 * t10i;
				  y10r = c6 * t6r + c7 * t7r + c8 * t8r + c9 * t9r + c10 * t10r;
				  y10i = c6 * t6i + c7 * t7i + c8 * t8i + c9 * t9i + c10 * t10i;
				z[j0] = z[j0] + t1r + t2r + t3r + t4r + t5r;
				z[j0+1] = z[j0+1] + t1i + t2i + t3i + t4i + t5i;
				z[j1] = y1r - y10i;
				z[j1+1] = y1i + y10r;
				z[j2] = y2r - y9i;
				z[j2+1] = y2i + y9r;
				z[j3] = y3r - y8i;
				z[j3+1] = y3i + y8r;
				z[j4] = y4r - y7i;
				z[j4+1] = y4i + y7r;
				z[j5] = y5r - y6i;
				z[j5+1] = y5i + y6r;
				z[j6] = y5r + y6i;
				z[j6+1] = y5i - y6r;
				z[j7] = y4r + y7i;
				z[j7+1] = y4i - y7r;
				z[j8] = y3r + y8i;
				z[j8+1] = y3i - y8r;
				z[j9] = y2r + y9i;
				z[j9+1] = y2i - y9r;
				z[j10] = y1r + y10i;
				z[j10+1] = y1i - y10r;
				  jt = j10 + 2;
				j10 = j9 + 2;
				j9 = j8 + 2;
				j8 = j7 + 2;
				j7 = j6 + 2;
				j6 = j5 + 2;
				j5 = j4 + 2;
				j4 = j3 + 2;
				j3 = j2 + 2;
				j2 = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}

		  j11 = j10 + jinc;
		if( j11 >= jmax )  j11 = j11 - jmax;
		  j12 = j11 + jinc;
		if( j12 >= jmax )  j12 = j12 - jmax;

		if( ifac == 13 )									
			{
//			  c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12;
			if( mu == 1 )
				{
				c1 = P885;
				c2 = P568;
				c3 = P120;
				c4 = -P354;
				c5 = -P748;
				c6 = -P970;
				c7 = P464;
				c8 = P822;
				c9 = P992;
				c10 = P935;
				c11 = P663;
				c12 = P239;
				}
			else if( mu == 2 )
				{
				c1 = P568;
				c2 = -P354;
				c3 = -P970;
				c4 = -P748;
				c5 = P120;
				c6 = P885;
				c7 = P822;
				c8 = P935;
				c9 = P239;
				c10 = -P663;
				c11 = -P992;
				c12 = -P464;
				}
			else if( mu == 3 )
				{
				c1 = P120;
				c2 = -P970;
				c3 = -P354;
				c4 = P885;
				c5 = P568;
				c6 = -P748;
				c7 = P992;
				c8 = P239;
				c9 = -P935;
				c10 = -P464;
				c11 = P822;
				c12 = P663;
				}
			else if( mu == 4 )
				{
				c1 = -P354;
				c2 = -P748;
				c3 = P885;
				c4 = P120;
				c5 = -P970;
				c6 = P568;
				c7 = P935;
				c8 = -P663;
				c9 = -P464;
				c10 = P992;
				c11 = -P239;
				c12 = -P822;
				}
			else if( mu == 5 )
				{
				c1 = -P748;
				c2 = P120;
				c3 = P568;
				c4 = -P970;
				c5 = P885;
				c6 = -P354;
				c7 = P663;
				c8 = -P992;
				c9 = P822;
				c10 = -P239;
				c11 = -P464;
				c12 = P935;
				}
			else if( mu == 6 )
				{
				c1 = -P970;
				c2 = P885;
				c3 = -P748;
				c4 = P568;
				c5 = -P354;
				c6 = P120;
				c7 = P239;
				c8 = -P464;
				c9 = P663;
				c10 = -P822;
				c11 = P935;
				c12 = -P992;
				}
			else if( mu == 7 )
				{
				c1 = -P970;
				c2 = P885;
				c3 = -P748;
				c4 = P568;
				c5 = -P354;
				c6 = P120;
				c7 = -P239;
				c8 = P464;
				c9 = -P663;
				c10 = P822;
				c11 = -P935;
				c12 = P992;
				}
			else if( mu == 8 )
				{
				c1 = -P748;
				c2 = P120;
				c3 = P568;
				c4 = -P970;
				c5 = P885;
				c6 = -P354;
				c7 = -P663;
				c8 = P992;
				c9 = -P822;
				c10 = P239;
				c11 = P464;
				c12 = -P935;
				}
			else if( mu == 9 )
				{
				c1 = -P354;
				c2 = -P748;
				c3 = P885;
				c4 = P120;
				c5 = -P970;
				c6 = P568;
				c7 = -P935;
				c8 = P663;
				c9 = P464;
				c10 = -P992;
				c11 = P239;
				c12 = P822;
				}
			else if( mu == 10 )
				{
				c1 = P120;
				c2 = -P970;
				c3 = -P354;
				c4 = P885;
				c5 = P568;
				c6 = -P748;
				c7 = -P992;
				c8 = -P239;
				c9 = P935;
				c10 = P464;
				c11 = -P822;
				c12 = -P663;
				}
			else if( mu == 11 )
				{
				c1 = P568;
				c2 = -P354;
				c3 = -P970;
				c4 = -P748;
				c5 = P120;
				c6 = P885;
				c7 = -P822;
				c8 = -P935;
				c9 = -P239;
				c10 = P663;
				c11 = P992;
				c12 = P464;
				}
			else
				{
				c1 = P885;
				c2 = P568;
				c3 = P120;
				c4 = -P354;
				c5 = -P748;
				c6 = -P970;
				c7 = -P464;
				c8 = -P822;
				c9 = -P992;
				c10 = -P935;
				c11 = -P663;
				c12 = -P239;
				}
			for(   l = 0; l < m; l++ )
				{
				  t1r = z[j1] + z[j12];
				  t1i = z[j1+1] + z[j12+1];
				  t2r = z[j2] + z[j11];
				  t2i = z[j2+1] + z[j11+1];
				  t3r = z[j3] + z[j10];
				  t3i = z[j3+1] + z[j10+1];
				  t4r = z[j4] + z[j9];
				  t4i = z[j4+1] + z[j9+1];
				  t5r = z[j5] + z[j8];
				  t5i = z[j5+1] + z[j8+1];
				  t6r = z[j6] + z[j7];
				  t6i = z[j6+1] + z[j7+1];
				  t7r = z[j1] - z[j12];
				  t7i = z[j1+1] - z[j12+1];
				  t8r = z[j2] - z[j11];
				  t8i = z[j2+1] - z[j11+1];
				  t9r = z[j3] - z[j10];
				  t9i = z[j3+1] - z[j10+1];
				  t10r = z[j4] - z[j9];
				  t10i = z[j4+1] - z[j9+1];
				  t11r = z[j5] - z[j8];
				  t11i = z[j5+1] - z[j8+1];
				  t12r = z[j6] - z[j7];
				  t12i = z[j6+1] - z[j7+1];
				  t13r = z[j0] - 0.5 * t6r;
				  t13i = z[j0+1] - 0.5 * t6i;
				  t14r = t1r - t6r;
				  t14i = t1i - t6i;
				  t15r = t2r - t6r;
				  t15i = t2i - t6i;
				  t16r = t3r - t6r;
				  t16i = t3i - t6i;
				  t17r = t4r - t6r;
				  t17i = t4i - t6i;
				  t18r = t5r - t6r;
				  t18i = t5i - t6i;
				  y1r = t13r + c1 * t14r + c2 * t15r + c3 * t16r + c4 * t17r + c5 * t18r;
				  y1i = t13i + c1 * t14i + c2 * t15i + c3 * t16i + c4 * t17i + c5 * t18i;
				  y2r = t13r + c2 * t14r + c4 * t15r + c6 * t16r + c5 * t17r + c3 * t18r;
				  y2i = t13i + c2 * t14i + c4 * t15i + c6 * t16i + c5 * t17i + c3 * t18i;
				  y3r = t13r + c3 * t14r + c6 * t15r + c4 * t16r + c1 * t17r + c2 * t18r;
				  y3i = t13i + c3 * t14i + c6 * t15i + c4 * t16i + c1 * t17i + c2 * t18i;
				  y4r = t13r + c4 * t14r + c5 * t15r + c1 * t16r + c3 * t17r + c6 * t18r;
				  y4i = t13i + c4 * t14i + c5 * t15i + c1 * t16i + c3 * t17i + c6 * t18i;
				  y5r = t13r + c5 * t14r + c3 * t15r + c2 * t16r + c6 * t17r + c1 * t18r;
				  y5i = t13i + c5 * t14i + c3 * t15i + c2 * t16i + c6 * t17i + c1 * t18i;
				  y6r = t13r + c6 * t14r + c1 * t15r + c5 * t16r + c2 * t17r + c4 * t18r;
				  y6i = t13i + c6 * t14i + c1 * t15i + c5 * t16i + c2 * t17i + c4 * t18i;
				  y7r = c12 * t7r - c7 * t8r + c11 * t9r - c8 * t10r + c10 * t11r - c9 * t12r;
				  y7i = c12 * t7i - c7 * t8i + c11 * t9i - c8 * t10i + c10 * t11i - c9 * t12i;
				  y8r = c11 * t7r - c9 * t8r + c8 * t9r - c12 * t10r - c7 * t11r + c10 * t12r;
				  y8i = c11 * t7i - c9 * t8i + c8 * t9i - c12 * t10i - c7 * t11i + c10 * t12i;
				  y9r = c10 * t7r - c11 * t8r - c7 * t9r + c9 * t10r - c12 * t11r - c8 * t12r;
				  y9i = c10 * t7i - c11 * t8i - c7 * t9i + c9 * t10i - c12 * t11i - c8 * t12i;
				  y10r = c9 * t7r + c12 * t8r - c10 * t9r - c7 * t10r + c8 * t11r + c11 * t12r;
				  y10i = c9 * t7i + c12 * t8i - c10 * t9i - c7 * t10i + c8 * t11i + c11 * t12i;
				  y11r = c8 * t7r + c10 * t8r + c12 * t9r - c11 * t10r - c9 * t11r - c7 * t12r;
				  y11i = c8 * t7i + c10 * t8i + c12 * t9i - c11 * t10i - c9 * t11i - c7 * t12i;
				  y12r = c7 * t7r + c8 * t8r + c9 * t9r + c10 * t10r + c11 * t11r + c12 * t12r;
				  y12i = c7 * t7i + c8 * t8i + c9 * t9i + c10 * t10i + c11 * t11i + c12 * t12i;
				z[j0] = z[j0] + t1r + t2r + t3r + t4r + t5r + t6r;
				z[j0+1] = z[j0+1] + t1i + t2i + t3i + t4i + t5i + t6i;
				z[j1] = y1r - y12i;
				z[j1+1] = y1i + y12r;
				z[j2] = y2r - y11i;
				z[j2+1] = y2i + y11r;
				z[j3] = y3r - y10i;
				z[j3+1] = y3i + y10r;
				z[j4] = y4r - y9i;
				z[j4+1] = y4i + y9r;
				z[j5] = y5r - y8i;
				z[j5+1] = y5i + y8r;
				z[j6] = y6r - y7i;
				z[j6+1] = y6i + y7r;
				z[j7] = y6r + y7i;
				z[j7+1] = y6i - y7r;
				z[j8] = y5r + y8i;
				z[j8+1] = y5i - y8r;
				z[j9] = y4r + y9i;
				z[j9+1] = y4i - y9r;
				z[j10] = y3r + y10i;
				z[j10+1] = y3i - y10r;
				z[j11] = y2r + y11i;
				z[j11+1] = y2i - y11r;
				z[j12] = y1r + y12i;
				z[j12+1] = y1i - y12r;
				  jt = j12 + 2;
				j12 = j11 + 2;
				j11 = j10 + 2;
				j10 = j9 + 2;
				j9 = j8 + 2;
				j8 = j7 + 2;
				j7 = j6 + 2;
				j6 = j5 + 2;
				j5 = j4 + 2;
				j4 = j3 + 2;
				j3 = j2 + 2;
				j2 = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}

		  j13 = j12 + jinc;
		if( j13 >= jmax )  j13 = j13 - jmax;
		  j14 = j13 + jinc;
		if( j14 >= jmax ) j14 = j14 - jmax;
		  j15 = j14 + jinc;
		if( j15 >= jmax ) j15 = j15 - jmax;

		if( ifac == 16 )									
			{
			  c1, c2, c3, c4;
			if( mu == 1 )
				{
				c1 = 1.0;
				c2 = P923;
				c3 = P382;
				c4 = P707;
				}
			else if( mu == 3 )
				{
				c1 =  -1.0;
				c2 = P382;
				c3 = P923;
				c4 = -P707;
				}
			else if( mu == 5 )
				{
				c1 = 1.0;
				c2 = -P382;
				c3 = P923;
				c4 = -P707;
				}
			else if( mu == 7 )
				{
				c1 =  -1.0;
				c2 = -P923;
				c3 = P382;
				c4 = P707;
				}
			else if( mu == 9 )
				{
				c1 = 1.0;
				c2 = -P923;
				c3 = -P382;
				c4 = P707;
				}
			else if( mu == 11 )
				{
				c1 =  -1.0;
				c2 = -P382;
				c3 = -P923;
				c4 = -P707;
				}
			else if( mu == 13 )
				{
				c1 = 1.0;
				c2 = P382;
				c3 = -P923;
				c4 = -P707;
				}
			else
				{
				c1 =  -1.0;
				c2 = P923;
				c3 = -P382;
				c4 = P707;
				}
			  c5 = c1 * c4;
			  c6 = c1 * c3;
			  c7 = c1 * c2;
			for(   l = 0; l < m; l++ )
				{
				  t1r = z[j0] + z[j8];
				  t1i = z[j0+1] + z[j8+1];
				  t2r = z[j4] + z[j12];
				  t2i = z[j4+1] + z[j12+1];
				  t3r = z[j0] - z[j8];
				  t3i = z[j0+1] - z[j8+1];
				  t4r = c1 * (z[j4] - z[j12]);
				  t4i = c1 * (z[j4+1] - z[j12+1]);
				  t5r = t1r + t2r;
				  t5i = t1i + t2i;
				  t6r = t1r - t2r;
				  t6i = t1i - t2i;
				  t7r = z[j1] + z[j9];
				  t7i = z[j1+1] + z[j9+1];
				  t8r = z[j5] + z[j13];
				  t8i = z[j5+1] + z[j13+1];
				  t9r = z[j1] - z[j9];
				  t9i = z[j1+1] - z[j9+1];
				  t10r = z[j5] - z[j13];
				  t10i = z[j5+1] - z[j13+1];
				  t11r = t7r + t8r;
				  t11i = t7i + t8i;
				  t12r = t7r - t8r;
				  t12i = t7i - t8i;
				  t13r = z[j2] + z[j10];
				  t13i = z[j2+1] + z[j10+1];
				  t14r = z[j6] + z[j14];
				  t14i = z[j6+1] + z[j14+1];
				  t15r = z[j2] - z[j10];
				  t15i = z[j2+1] - z[j10+1];
				  t16r = z[j6] - z[j14];
				  t16i = z[j6+1] - z[j14+1];
				  t17r = t13r + t14r;
				  t17i = t13i + t14i;
				  t18r = c4 * (t15r - t16r);
				  t18i = c4 * (t15i - t16i);
				  t19r = c5 * (t15r + t16r);
				  t19i = c5 * (t15i + t16i);
				  t20r = c1 * (t13r - t14r);
				  t20i = c1 * (t13i - t14i);
				  t21r = z[j3] + z[j11];
				  t21i = z[j3+1] + z[j11+1];
				  t22r = z[j7] + z[j15];
				  t22i = z[j7+1] + z[j15+1];
				  t23r = z[j3] - z[j11];
				  t23i = z[j3+1] - z[j11+1];
				  t24r = z[j7] - z[j15];
				  t24i = z[j7+1] - z[j15+1];
				  t25r = t21r + t22r;
				  t25i = t21i + t22i;
				  t26r = t21r - t22r;
				  t26i = t21i - t22i;
				  t27r = t9r + t24r;
				  t27i = t9i + t24i;
				  t28r = t10r + t23r;
				  t28i = t10i + t23i;
				  t29r = t9r - t24r;
				  t29i = t9i - t24i;
				  t30r = t10r - t23r;
				  t30i = t10i - t23i;
				  t31r = t5r + t17r;
				  t31i = t5i + t17i;
				  t32r = t11r + t25r;
				  t32i = t11i + t25i;
				  t33r = t3r + t18r;
				  t33i = t3i + t18i;
				  t34r = c2 * t29r - c6 * t30r;
				  t34i = c2 * t29i - c6 * t30i;
				  t35r = t3r - t18r;
				  t35i = t3i - t18i;
				  t36r = c7 * t27r - c3 * t28r;
				  t36i = c7 * t27i - c3 * t28i;
				  t37r = t4r + t19r;
				  t37i = t4i + t19i;
				  t38r = c3 * t27r + c7 * t28r;
				  t38i = c3 * t27i + c7 * t28i;
				  t39r = t4r - t19r;
				  t39i = t4i - t19i;
				  t40r = c6 * t29r + c2 * t30r;
				  t40i = c6 * t29i + c2 * t30i;
				  t41r = c4 * (t12r - t26r);
				  t41i = c4 * (t12i - t26i);
				  t42r = c5 * (t12r + t26r);
				  t42i = c5 * (t12i + t26i);
				  y1r = t33r + t34r;
				  y1i = t33i + t34i;
				  y2r = t6r + t41r;
				  y2i = t6i + t41i;
				  y3r = t35r + t40r;
				  y3i = t35i + t40i;
				  y4r = t5r - t17r;
				  y4i = t5i - t17i;
				  y5r = t35r - t40r;
				  y5i = t35i - t40i;
				  y6r = t6r - t41r;
				  y6i = t6i - t41i;
				  y7r = t33r - t34r;
				  y7i = t33i - t34i;
				  y9r = t38r - t37r;
				  y9i = t38i - t37i;
				  y10r = t42r - t20r;
				  y10i = t42i - t20i;
				  y11r = t36r + t39r;
				  y11i = t36i + t39i;
				  y12r = c1 * (t11r - t25r);
				  y12i = c1 * (t11i - t25i);
				  y13r = t36r - t39r;
				  y13i = t36i - t39i;
				  y14r = t42r + t20r;
				  y14i = t42i + t20i;
				  y15r = t38r + t37r;
				  y15i = t38i + t37i;
				z[j0] = t31r + t32r;
				z[j0+1] = t31i + t32i;
				z[j1] = y1r - y15i;
				z[j1+1] = y1i + y15r;
				z[j2] = y2r - y14i;
				z[j2+1] = y2i + y14r;
				z[j3] = y3r - y13i;
				z[j3+1] = y3i + y13r;
				z[j4] = y4r - y12i;
				z[j4+1] = y4i + y12r;
				z[j5] = y5r - y11i;
				z[j5+1] = y5i + y11r;
				z[j6] = y6r - y10i;
				z[j6+1] = y6i + y10r;
				z[j7] = y7r - y9i;
				z[j7+1] = y7i + y9r;
				z[j8] = t31r - t32r;
				z[j8+1] = t31i - t32i;
				z[j9] = y7r + y9i;
				z[j9+1] = y7i - y9r;
				z[j10] = y6r + y10i;
				z[j10+1] = y6i - y10r;
				z[j11] = y5r + y11i;
				z[j11+1] = y5i - y11r;
				z[j12] = y4r + y12i;
				z[j12+1] = y4i - y12r;
				z[j13] = y3r + y13i;
				z[j13+1] = y3i - y13r;
				z[j14] = y2r + y14i;
				z[j14+1] = y2i - y14r;
				z[j15] = y1r + y15i;
				z[j15+1] = y1i - y15r;
				  jt = j15 + 2;
				j15 = j14 + 2;
				j14 = j13 + 2;
				j13 = j12 + 2;
				j12 = j11 + 2;
				j11 = j10 + 2;
				j10 = j9 + 2;
				j9 = j8 + 2;
				j8 = j7 + 2;
				j7 = j6 + 2;
				j6 = j5 + 2;
				j5 = j4 + 2;
				j4 = j3 + 2;
				j3 = j2 + 2;
				j2 = j1 + 2;
				j1 = j0 + 2;
				j0 = jt;
				}
			continue;
			}
	}
}

float  *FloatVector( int nl, int nh )
{
float  *v;

	v = ( float * )malloc( ( int )( nh-nl+1 ) * 4L );
	return( v-nl );
}
int  *IntVector( int nl, int nh )
{
int  *v;

	v = ( int * )malloc( ( int )( nh-nl+1 ) * 4L );
	return( v-nl );
}
int  Power2( int samples )
{
int   tpow[15] = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768 };

	int  i = 0;
	while( tpow[i] < samples )  i++;
	return( tpow[i] );
}

int  ludcmp( float *a, int n, int *indx, float *vv )
{
int  i, imax, j, k;
float  big, dum, sum, temp;

	for( i = 0; i < n; i++ )
		{
		big = 0.0;
		for( j = 0; j < n; j++ )
			{
			if( ( temp = fabs( a[i*n+j] ) ) > big )
				big = temp;
			}
		if( big == 0.0 )  return( 1 );
		vv[i] = 1.0 / big;
		}
	for( j = 0; j < n; j++ )
		{
		for( i = 0; i < j; i++ )
			{
			sum = a[i*n+j];
			for( k = 0; k < i; k++ )
				sum -= a[i*n+k] * a[k*n+j];
			a[i*n+j] = sum;
			}
		big = 0.0;
		for( i = j; i < n; i++ )
			{
			sum = a[i*n+j];
			for( k = 0; k < j; k++ )
				sum -= a[i*n+k] * a[k*n+j];
			a[i*n+j] = sum;
			if( ( dum = vv[i] * fabs( sum ) ) >= big )
				{
				big = dum;
				imax = i;
				}
			}
		if( j != imax )
			{
			for( k = 0; k < n; k++ )
				{
				dum = a[imax*n+k];
				a[imax*n+k] = a[j*n+k];
				a[j*n+k] = dum;
				}
			vv[imax] = vv[j];
			}
		indx[j] = imax;
		if( a[j*n+j] == 0.0 )
			a[j*n+j] = TINY;

		if( j != n )
			{
			dum = 1.0 / a[j*n+j];
			for( i = j+1; i < n; i++ )
				a[i*n+j] *= dum;
			}
		}
/*	free1float( vv );*/
	return( 0 );
}
void  lubksb( float *a, int n, int *indx, float *b )
{
int  i, ii = -1, ip, j;
float  sum;

	for( i = 0; i < n; i++ )
		{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if( ii >= 0 )
			{
			for( j = ii; j <= i-1; j++ )
				sum -= a[i*n+j] * b[j];
			}
		else if( sum )
			ii = i;
		b[i] = sum;
		}
	for( i = n-1; i >= 0; i-- )
		{
		sum = b[i];
		for( j = i+1; j < n; j++ )
			sum -= a[i*n+j] * b[j];
		b[i] = sum / a[i*n+i];
		}
}
