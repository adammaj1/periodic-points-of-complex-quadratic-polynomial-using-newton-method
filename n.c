/* 

gcc -std=c99 -Wall -Wextra -pedantic -O3 n.c -lm
./a.out

"The "flowers" are beautiful to behold but totally abhorrent from the numerical point of view." Ramillies
https://math.stackexchange.com/questions/2407659/why-does-the-newton-raphson-method-not-converge-for-some-functions


-----------------------------------
https://arxiv.org/abs/1703.05847





Newton's method in practice II: The iterated refinement Newton method and near-optimal complexity for finding all roots of some polynomials of very large degrees
Marvin Randig, Dierk Schleicher, Robin Stoll
(Submitted on 16 Mar 2017 (v1), last revised 31 Dec 2017 (this version, v2))
We present a practical implementation based on Newton's method to find all roots of several families of complex polynomials of degrees exceeding one billion (109) so that the observed complexity to find all roots is between O(dlnd) and O(dln3d) (measuring complexity in terms of number of Newton iterations or computing time). All computations were performed successfully on standard desktop computers built between 2007 and 2012.



------------------------------------
cd existing_folder
git init
git remote add origin git@gitlab.com:adammajewski/periodic-points-of-complex-quadratic-polynomial-using-newton-method.git
git add .
git commit -m "Initial commit"
git push -u origin master





*/

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		// strcat



/* ----------------------- */
int period = 2;
int degree; // = 2^period
int sMax; // = 4*d = number of starting points
long double dt; //  = 1.0 / sMax

// failure: we have n > Mfail, where Mfail is the largest number of allowed iterations before the algorithm gives up
int nMax; // = 1000000; // maxima number of Newton iterations  = 10 * d = M_fail 
int maximal_n = 0;

complex long double c = 0.0+1.0*I ; // dendrit Julia set

const long double pi = 3.1415926535897932384626433832795029L;









// viewport is a square : iWidth = iHeight and ( ZxMax-ZxMin) = (ZyMax-ZyMin)
// square in dyanamic ( complex plane) -> virtual ( memory) 2D array -> virtual 1D array -> pgm file 

/* coordinate in world units  */
static const long double ZxMin = -2.1;	// slightly more than radius
static const long double ZxMax = 2.1;	//
static const long double ZyMin = -2.1;	//
static const long double ZyMax = 2.1;	//
 long double PixelWidth;	// =(ZxMax-ZxMin)/ixMax;
long double PixelHeight;	// =(ZyMax-ZyMin)/iyMax;
static long double ratio;

complex long double z0 ; // initial aproximation for Newton method; 
complex long double zp; // periodic point



// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
//unsigned int ix, iy; // var
int ixMin = 0;	// Indexes of array starts from 0 not 1
int ixMax;	//
int iWidth;	// horizontal dimension of array

int iyMin = 0;	// Indexes of array starts from 0 not 1
int iyMax;	//

int iHeight = 1600;	//  


// The size of 1D array has to be a positive constant integer 
int iSize;	// = iWidth*iHeight; 



// memmory 1D array 
unsigned char *data;
// iMax = iSize-1; // Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].

//long double ER = 2.1;
//long double ER2; //  = 2.0 * 2.0; // ER*ER
long double radius = 2.0;

// epsilon_succes = success: we have |zn+1−zn| ≤ εsuccess, where εsuccess is a predefined accuracy threshold for success; 
long double EPS = 1e-18; // to much for double so I use long double 
long double EPS2 ;// EPS*EPS

// the minimal distance between any two of them was 5.47 · 10^−11
long double m_dist = 1e-10;
long double m_dist2 ; // = m_dist*m_dist
long double MinimalDistanceBetweenRoots; // minimal distance between any two




// arrays for complex points 
complex long double * zzs; // array for found points ( length = sMax)
// distinct roots = not equal 
complex long double * zzd; // array for periodic points = distinct points , where  length = degree -1  


// counters 
//int point_errors = 0;
//int point_drawn = 0;
int distinc_points = 0; // number of distinct roots = periodic points
long double VieteSum;








/*

https://stackoverflow.com/questions/40273079/find-all-prime-factors-of-a-given-number
https://en.wikipedia.org/wiki/Trial_division
*/

int GiveAllPrimeFactors( int n )
  {
	int a; 
	printf("prime factors of %d  = ", n);
  	for (a=2; a<=n; ++a)
  	
    		while(n%a==0)
    		{
      			printf("%d\t", a);
      			n = n/a;
    		}
    		
    	printf("\n");
  	
  return 0;

}





/* -----------  array functions -------------- */


/* gives position of 2D point (ix,iy) in 1D array  ; uses also global variable iWidth */
int Give_i( int ix,  int iy)
{ int i;
  i = ix + iy*iWidth; 
   if (i>=0 && i<iSize) 
   	return i;
   	else {printf("error from Give_i: i = %d ix = %d iy = %d iSize = %d \n", i, ix, iy, iSize); return -1;}
  
  }




// plots raster point (ix,iy) 
int iDrawPoint( int ix,  int iy, unsigned char iColor, unsigned char A[])
{ 
	int i;
  	i =  Give_i(ix,iy ); // compute index of 1D array from indices of 2D array 
 	if (i>=0 && i<iSize) 
 		A[i] = iColor; // draw
 		else {printf("error from iDrawPoint\n"); return 1;} // not draw

return 0; // succes
}


// draws point to memmory array data
// uses complex type so #include <complex.h> and -lm 
// The roots are marked by  big spots with inverted color  
int dDrawBigPoint_Inverted(complex long double point, unsigned char A[] )
{

  int ix, iy; // screen coordinate = indices of virtual 2D array
  int ixCenter, iyCenter;
  unsigned int i; // index of 1D array
  long double zx = creall(point);
  long double zy = cimagl(point);
  //int r;
  unsigned char iColor;
  int iSide =ixMax/500; /* half of width or height of big pixel */
  
  
  if (zx< ZxMax && zx > ZxMin && zy > ZyMin && zy < ZyMax){
  
  	ixCenter = (int) ((zx- ZxMin)/PixelWidth); 
  	iyCenter = (int) ((ZyMax - zy)/PixelHeight); // inverse Y axis 
  	//
  	i = Give_i(ixCenter,iyCenter);
  	// find color
  	iColor = A[i]; 
  	// inverse color
  	if (iColor < 150 && iColor>100) // if color is near 255/2 then 255-color is similar
  		iColor = 254; // more different 
  		else iColor = 255 - iColor ; // typical inverse
  	//printf(" color A[i] = %d and inverse = %d\n", A[i], 255-A[i]);
  	/* mark point by big pixel */
 
  	for(iy=iyCenter-iSide; iy<=iyCenter+iSide; iy++) 
    		for(ix=ixCenter-iSide; ix<=ixCenter+iSide; ix++){ 
      			iDrawPoint(ix, iy, iColor, A);
  			
  			}
  		
        
  	}
  	
  
  
  return 0; //  
}

int dDrawBigPoint(complex long double point, unsigned char A[], unsigned char iColor )
{

  int ix, iy; // screen coordinate = indices of virtual 2D array
  int ixCenter, iyCenter;
  long double zx = creall(point);
  long double zy = cimagl(point);
 
 
  int iSide =ixMax/500; /* half of width or height of big pixel */
  
  
  if (zx< ZxMax && zx > ZxMin && zy > ZyMin && zy < ZyMax){
  
  	ixCenter = (int) ((zx- ZxMin)/PixelWidth); 
  	iyCenter = (int) ((ZyMax - zy)/PixelHeight); // inverse Y axis 
  	 	  	
  	/* mark point by big pixel */
   	for(iy=iyCenter-iSide; iy<=iyCenter+iSide; iy++) 
    		for(ix=ixCenter-iSide; ix<=ixCenter+iSide; ix++){ 
      			iDrawPoint(ix, iy, iColor, A);
  			
  			}
  		
        
  	}
  	
  
  
  return 0; // r 
}



/*
http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm
Instead of swaps in the initialisation use error calculation for both directions x and y simultaneously:
*/
void iDrawLine( int x0, int y0, int x1, int y1, unsigned char iColor, unsigned char A[]) 
{
  int x=x0; 
  int y=y0;
  int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
  int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1; 
  int err = (dx>dy ? dx : -dy)/2, e2;

  for(;;){
    iDrawPoint(x, y, iColor, A);
    if (x==x1 && y==y1) break;
    e2 = err;
    if (e2 >-dx) { err -= dy; x += sx; }
    if (e2 < dy) { err += dx; y += sy; }
  }
}




int dDrawLine(long double Zx0, long double Zy0, long double Zx1, long double Zy1, unsigned char color, unsigned char A[]) 
{

 int ix0, iy0; // screen coordinate = indices of virtual 2D array 
 int ix1, iy1; // screen coordinate = indices of virtual 2D array

   // first step of clipping
   //if (  Zx0 < ZxMax &&  Zx0 > ZxMin && Zy0 > ZyMin && Zy0 <ZyMax 
    //  && Zx1 < ZxMax &&  Zx1 > ZxMin && Zy1 > ZyMin && Zy1 <ZyMax )
   ix0= (Zx0- ZxMin)/PixelWidth; 
   iy0 = (ZyMax - Zy0)/PixelHeight; // inverse Y axis 
   ix1= (Zx1- ZxMin)/PixelWidth; 
   iy1= (ZyMax - Zy1)/PixelHeight; // inverse Y axis 
   // second step of clipping
   if (ix0 >=ixMin && ix0<=ixMax && ix0 >=ixMin && ix0<=ixMax && iy0 >=iyMin && iy0<=iyMax 
      && iy1 >=iyMin && iy1<=iyMax )
   iDrawLine(ix0,iy0,ix1,iy1,color, A) ;

return 0;
}




int cDrawLine(complex long double z0, complex long double z1, unsigned char color, unsigned char A[]){

	dDrawLine(creall(z0), cimagl(z0), creall(z1), cimagl(z1), color, A);
	return 0; 


}



unsigned char GiveColor (complex long double z, unsigned char A[]){

	int ix, iy; // screen coordinate = indices of virtual 2D array 
 	int i; //  indices of 1D array
 	unsigned char iColor;
 	
 	ix = (int) ((creall(z)- ZxMin)/PixelWidth); 
        iy = (int) ((ZyMax - cimagl(z))/PixelHeight); // inverse Y axis 
        i = Give_i(ix,iy);
        
        // find color
  	iColor = A[i]; 
  	// inverse color
  	if (iColor < 150 && iColor>100) // if color is near 255/2 then 255-color is similar
  		iColor = 254; // more different 
  		else iColor = 255 - iColor ; // typical inverse
        return iColor;
        


}









int FillArray(unsigned char A[]){

	int i;
	for (i=0; i<iSize; i++){
	
		A[i] = 255; // white
	}
	return 0;



}





// save data array to pgm file 
int SaveArray2PGMFile( unsigned char A[], int p, int s, int k, char* comment )
{
  
  FILE * fp;
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  char name [100]; /* name of file */
  snprintf(name, sizeof name, "%d_%d_%d", p,s, k); /* period, sMax */
  char *filename =strncat(name,".pgm", 4);
  
  
  
  /* save image to the pgm file  */      
  fp= fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n # %s\n %u %u\n %u\n", comment, iWidth, iHeight, MaxColorComponentValue);  /*write header to the file*/
  fwrite(A,iSize,1,fp);  /*write image data bytes to the file in one step */
  
  //
  printf("File %s saved. ", filename);
  if (comment == NULL)  printf ("empty comment \n");
                   else printf ("%s \n", comment); 
  fclose(fp);

  return 0;
}





/* --------------- complex plane --------------------------- */



long double complex Give_z0(long double InternalAngleInTurns, long double radius )
{
  
  //0 <= InternalAngleInTurns <=1
  long double a = InternalAngleInTurns *2.0*pi; // from turns to radians
  long double Cx, Cy; /* C = Cx+Cy*i */
  
  
      Cx = radius*cosl(a); 
      Cy = radius*sinl(a); 
      

  return Cx + Cy*I;
}





long double cabs2(complex long double z) {
  return (creall(z) * creall(z) + cimagl(z) * cimagl(z));
}






/* 
newton function 
 N(z) = z - (fp(z)-z)/fp'(z)) 
 used for computing periodic point 
 of complex quadratic polynomial
 f(z) = z*z +c

*/

complex long double N( complex long double c, complex long double zn , int pMax ){

  
complex long double z = zn;
complex long double d = 1.0; /* d = first derivative with respect to z */
int p;

for (p=0; p < pMax; p++){

   
   d = 2*z*d; /* first derivative with respect to z */
   z = z*z +c ; /* complex quadratic polynomial */
  
    
}

      z = zn - (z - zn)/(d - 1) ; // Newton function

return z;
}




/* 
compute periodic point 
of complex quadratic polynomial1
using Newton iteration = numerical method

*/

complex long double GivePeriodicAndDrawRay(complex long double c, complex long double z0, int period, long double eps2, unsigned char A[]){

	complex long double z = z0;
	complex long double zPrev = z0; // prevoiuus value of z
	int n ; // iteration
	unsigned char iColor;


	iColor = GiveColor(z0, A);
	for (n=0; n<nMax; n++) {
     
    		z = N( c, z, period);
   
    
  
    
    		cDrawLine(zPrev, z, iColor, A);
   		if (cabs2(z - zPrev)< eps2) break; // succes 
    
   		 zPrev = z; }
    
   	if ( n>maximal_n) maximal_n = n; // 
     
  	// printf(" n = %d\n", n); 
  	
  	dDrawBigPoint(z, A, iColor);

	return z;
}




int DrawRays(unsigned char A[]){
	int s; // number of starting point z0
  	long double t= 0.0; //
  	
  
  	complex long double z0;
  	complex long double zp;

	for (s=0; s<sMax; s++ ){
  		z0 = Give_z0(t, radius); // compute new starting point on the circle 
  		
  		 
  		zp = GivePeriodicAndDrawRay( c , z0, period, EPS2, A); // compute periodic point 
  		zzs[s] = zp; // save all found points 
  		t += dt; // next angle using globl variable dt
  	}
  	
	return 0;

}






complex long double GivePeriodic(complex long double c, complex long double z0, int period, long double eps2){

	complex long double z = z0;
	complex long double zPrev = z0; // prevoiuus value of z
	int n ; // iteration




	for (n=0; n<nMax; n++) {
     
		z = N( c, z, period);
   
    		//cDrawLine(zPrev, z, 100, data);
    		if (cabs2(z - zPrev)< eps2) break; // succes 
    
    		zPrev = z; }
    
    	if ( n>maximal_n) maximal_n = n; // 
     
  	// printf(" n = %d\n", n); 

	return z;
}



/*

check if point z is different 
from all first dMax points of zzd array

uses global :
* array zzd
* limit m_dist2


*/
int IsDistinct(complex long double z, int d_length){
	int d; 
	int dMax = d_length+1; // d starts from 0 to dMax-1
	
	long double distance2 = 0.0L;
	
	for (d=0; d<dMax; d++ ){
	 	distance2 = cabs2(z - zzd[d]);
		 if (distance2< m_dist2) // the mutual distance between pairs of periodic points: z and zzd[d]
		 	return 0; // no = not distinct = not different 
		}
		
		
	//printf("distance = %Lf\n", sqrtl(distance2));	 	
	return 1; // distinct 
		 	



}





/* 

find distinct point in zzs and moves it to zzd 

take point from zzs
check if it is distinct
if yes then move it to zzd
if no go to next point from zzs

uses global variables : 
* arrays : zzs i zzd
* counter : distinc_points

*/

int FindDistinctPoints(){


	int s = 0; // index of zzs array 
	int d = 0; // index of zzd array
	complex long double z;
	
	
	//if ( iLength != sMax) {printf(" error : length(zzs) != sMax \n"); return 1;}

	// first point without checking 
	z = zzs[s]; // take point z from zzs
	zzd[d] = z; // move to zzd 
	d +=1;
	
	// rest of points
	for (s=1; s<sMax; s++ ){
		z = zzs[s]; // take next point
		// 
		if (IsDistinct(z,d)){
			zzd[d] = z; // move to zzd 
			d +=1;
		}
			
	
	}
	
	distinc_points = d ; // update counter = number of all distinct points 
	
	
	
	
	return 0; 




}



// from screen to world coordinate ; linear mapping
// uses global cons
complex long double GiveZ( int ix, int iy){ 
	long double zx = ZxMin + ix*PixelWidth;
	long double zy = ZyMax - iy*PixelHeight; // reverse y axis
   	return (zx + zy*I);

}

/*
	input : 
		zp = periodic point 
		zzd = table of periodic points
	output :
		index of table zzd

 each periodic points has it's own basin
 find it's number'

*/
int GiveBasinNumber(complex long double zp, complex long double zzd[]){

	int d = 0; // index of zzd array
	int dMax = distinc_points;
	
	
	// 
	for (d=0; d<dMax; d++ )
		if (cabs2(zp - zzd[d])<m_dist2) // the mutual distance between pairs of periodic points: z and zzd[d])
			return d; // 
	return -1; 

}



unsigned char GiveBasinColor(complex long double z0){

	
	complex long double zp;
	int i; 
	int step = (int)((250 - 10)/distinc_points);
	
	zp = GivePeriodic( c , z0, period, EPS2); // compute periodic point using Newton method
	i = GiveBasinNumber(zp, zzd); // color is proportional to Number of Newton Basin in zzd array
	
	return 50 + i*step; // return 8 bit color = shades of gray
	

}

//
int DrawNewtonBasins(unsigned char A[]){

	int ix, iy;
	int i;
	complex long double z;
	
	
	// for all points (x,y) of the image ( 2D array)
	for(iy=0;iy<iyMax;iy++)
 		for(ix=0;ix<ixMax;ix++){ 
 			     
 			z = GiveZ(ix, iy);  // compute pixel coordinate   
 			i = Give_i(ix, iy); // compute idex of 1D array	
			A[i] = GiveBasinColor(z); // compute color and save it to 1D array         
			}
	return 0;



}



int DrawRoots(unsigned char A[]){
	int d;
	for (d=0; d<distinc_points; d++ )
 		dDrawBigPoint_Inverted( zzd[d], A);
 		
 	return 0;



}



int DrawRootsColor(unsigned char A[], unsigned char iColor){
	int d;
	for (d=0; d<distinc_points; d++ )
 		dDrawBigPoint( zzd[d], A, iColor);
 		
 	return 0;



}


// https://en.wikipedia.org/wiki/Vieta%27s_formulas

long double ComputeVieteSum(){
	complex long double sum = 0;
	int d;
	for (d=0; d<distinc_points; d++ )
		sum += zzd[d]; // zzd is a array of roots (global variable)
	return cabsl(sum);
	
	



}

/* 

https://en.wikipedia.org/wiki/Closest_pair_of_points_problem
*/
long double ComputeMinimalDIstanceBetweenRoots(){

	long double minimal_dist2 = 2.0L; // minimal_dist^2
	long double minimal_dist;
	long double dist2;
	int i,j;
	complex long double z1 = 100;
	complex long double z2 = 200;
	
	
	for (i=0; i<distinc_points-1; i++ )
		for (j=i+1; j<distinc_points; j++ ){
			dist2 = cabs2(zzd[i] - zzd[j]);
			if (dist2 < minimal_dist2) {
				minimal_dist2 = dist2;
				z1 = zzd[i];
				z2 = zzd[j];
				}
		 	}
	
	minimal_dist = sqrtl(minimal_dist2);
	printf(" minimal distnce =%Le between\nz1 = %+.18Lf ; %+.18Lf \nz2 = %+.18Lf ; %+.18Lf\n ", minimal_dist, creall(z1), cimagl(z1), creall(z2), cimagl(z2));
	//dDrawBigPoint_Inverted( z1, data);
	//dDrawBigPoint_Inverted( z2, data);
	return minimal_dist;


}






int setup()
{	
	// numerical optimisaton for cabs2
	EPS2 = EPS * EPS;
	m_dist2 = m_dist*m_dist;
	
	degree = (int) pow(2.0, period);// degree = 2^period
	sMax = 2*degree; // 4*d = number of starting points
	nMax = 10* degree + 100; // 
	dt = 1.0/ sMax;
	
	
	// 
	iWidth = iHeight; // square
	iyMax = iHeight -1;
	ixMax = iWidth - 1; 
	iSize = iWidth * iHeight;
	//
	PixelWidth = (ZxMax-ZxMin)/ixMax;
        PixelHeight = (ZyMax-ZyMin)/iyMax;
	
	ratio = ( ZxMax - ZxMin)/ (ZyMax-ZyMin);
	
	
	
	//  dynamic 1D array for the image 
	data = malloc( iSize * sizeof(unsigned char) );
  	if (data == NULL )
    		{
      			fprintf(stderr," Could not allocate memory for data\n");
      			return 1;
    		}
    		
    		
    		
    	//
    	zzs = malloc((sMax) * sizeof(complex long double));	
	if (zzs == NULL )
    		{
      			fprintf(stderr," Could not allocate memory for zzs \n");
      			return 1;
    		}
    		
    	//
    	zzd = malloc(degree * sizeof(complex long double));	
	if (zzd == NULL )
    		{
      			fprintf(stderr," Could not allocate memory for zzd \n");
      			return 1;
    		}
	
	
	FillArray(data); // make all points white
	return 0;
}




int info(){
	int d; 

	printf (" parameter c from fc(z) = z^2+c is c = %Lf ; %Lf \n",  creall(c), cimagl(c)); 
	printf ("\n");
	printf (" period = %d  \n",  period); 
	
	
	printf (" degree of polynomial = %d  = 2^period\n", degree); 
	
	GiveAllPrimeFactors(period);
	printf (" number of roots = number of periodic points = degree of polynomial = %d  \n", degree); 
	printf (" number of starting points sMax = %d\n", sMax);
	if ( sMax<degree) printf("sMax < degree !!! it should be equal or greater \n");
	if (distinc_points < degree) 
		printf(" only %d from %d distinct points are found, increase sMax \n ", distinc_points, degree);
		else printf("succes : all %d distinct points are found !!\n", distinc_points);
	
	printf ("\n");
	printf (" dt = %Le\n", dt);
	printf (" radius of the circle around all periodic points = %Lf\n", radius);
	printf (" \n");
	printf (" maximal allowed number of Newton iterations nMax = %d  =  10*degree + 100, see setup\n", nMax);
	printf (" maximal used number of Newton iterations maximal_n = %d \n", maximal_n);
	if (nMax == maximal_n) printf(" possible error : nMax == maximal_n ; increase nMax \n");
	
	printf ("stopping criterion for the Newton iteration is epsilon_stop = %Le\n", EPS);
	
  	
  	
  	//if (point_errors>0) printf ("point_errors = %d\n", point_errors);
  	//printf ("point_drawn = %d\n", point_drawn);
  	
  	
  	// print only 
  	/*
	printf("\n\n all points:\n");
	for (s=0; s<sMax; s++ ) 
		printf(" s = %d z = %.18Lf; %.18Lf \n", s, creall(zzs[s]), cimagl(zzs[s]));	
	printf("\n\n\n");
	*/
	//  only distinct points from zzd array
	printf("\n\n periodic points are: \n");
	for (d=0; d<distinc_points; d++ )
		printf(" d = %2d z = %+.18Lf; %+.18Lf \n", d, creall(zzd[d]), cimagl(zzd[d]));
	printf("\n\n\n");
  	
  	
  	VieteSum = ComputeVieteSum();
  	printf(" the sum of all roots should be zero by Viete’s formula (this sum should be the negative of the degree d − 1 coefficient)\n");
  	printf(" Viete sum = %.18Le ( it should be zero )\n", VieteSum);
  	
  	MinimalDistanceBetweenRoots = ComputeMinimalDIstanceBetweenRoots();
  	
  	
	return 0;
}


int end(){
	info();
	// free memory
        free(data);
        free(zzs);
        free(zzd);
        
        
	return 0;
}


/*===================== main ======================================================= */ 

int main() {

 
 
 
  
  	
  	setup();
  
  
  	// create rays image 
  	DrawRays(data);
    
  	FindDistinctPoints();
  			
  	SaveArray2PGMFile( data, period, sMax, 11,  "Newton iterations (rays)");
  
  	//
  	if (degree < 255 ) { // there are only 255 colors !!
  		printf("Drawing Newton Basins \n");
  		DrawNewtonBasins(data);
   		SaveArray2PGMFile( data, period, sMax, 12,  "Newton Basins");
   		
   		DrawRays(data);
   		SaveArray2PGMFile( data, period, sMax, 13,  "Newton Basins and rays");
   		
  	  	
  		}
  
  	// 
  	FillArray(data); // make all points white = 255
  	DrawRootsColor(data, 0); // draw black roots on white background
  	SaveArray2PGMFile( data, period, sMax, 14,  "only roots");
  
  
  
  
  	//
  	end();
  
 
  
  
  
  return 0;
}
