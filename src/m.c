/* 

gcc -std=c99 -Wall -Wextra -pedantic -O3 n.c -lm
./a.out



--------- text output -------------
 period = 12  
 degree of polynomial = 4096  
 number of starting points sMax = 400
 dt = 2.500000e-03
 radius of the circle = 2.000000
 maximal allowed number of Newton iterations nMax = 40960 ( 10*degree) 
 maximal used number of Newton iterations maximal_n = 3337  
epsilon = 1.000000e-06
 c = 0.000000 ; 1.000000 
nan_errors = 0
point_errors = 0
point_drawn = 400


-----------------------------------
https://arxiv.org/abs/1703.05847





Newton's method in practice II: The iterated refinement Newton method and near-optimal complexity for finding all roots of some polynomials of very large degrees
Marvin Randig, Dierk Schleicher, Robin Stoll
(Submitted on 16 Mar 2017 (v1), last revised 31 Dec 2017 (this version, v2))
We present a practical implementation based on Newton's method to find all roots of several families of complex polynomials of degrees exceeding one billion (109) so that the observed complexity to find all roots is between O(dlnd) and O(dln3d) (measuring complexity in terms of number of Newton iterations or computing time). All computations were performed successfully on standard desktop computers built between 2007 and 2012.




Figure 2. The Newton iteration for 400 starting points on a circle
of radius r = 2 (here for a polynomial of degree 4096, so we do not
have enough starting points to find all roots; the polynomial shown
here describes periodic points of period dividing 12 of z 7→ z
2 + i).
The apparent lines connect orbits under the Newton dynamics; colors
indicate the number of iterations until an approximate root is
found. The behavior of the iterations outside of the disk containing
all roots is very parallel and “wasteful”, but required to carry
over the control from the circle of starting points to the interesting
dynamics on D. Observe that even if we had a very precise bound
on the smallest disk containing all roots, this would not help much
as most roots are away from the boundary of this dis

----------------------
convert 12.pgm -resize 600x600 12.png

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

int iHeight = 2000;	//  


// The size of 1D array has to be a positive constant integer 
int iSize;	// = iWidth*iHeight; 



// memmory 1D array 
unsigned char *data;
// iMax = iSize-1; // Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].

//long double ER = 2.1;
//long double ER2; //  = 2.0 * 2.0; // ER*ER
long double radius = 2.0;
long double EPS = 1e-6; // to much for double so I use long double 
long double EPS2 ;// EPS*EPS




int period = 12;
int degree; // = 2^period
int sMax; // = 4*d = number of starting points
long double dt; //  = 1.0 / sMax
int nMax; // = 1000000; // maxima number of Newton iterations  = 10 * d
int maximal_n = 0;

complex long double c = 0.0+1.0*I ; // dendrit Julia set

const long double pi = 3.1415926535897932384626433832795029L;



int nan_errors = 0;
int point_errors = 0;
int point_drawn = 0;






















/* -----------  array functions -------------- */


/* gives position of 2D point (ix,iy) in 1D array  ; uses also global variable iWidth */
int Give_i( int ix,  int iy)
{ int i;
  i = ix + iy*iWidth; 
   if (i>0 && i<iSize) 
   	return i;
   	else {printf("error from Give_i\n"); return -1;}
  
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
int dDrawPoint(complex long double point,unsigned char iColor, unsigned char A[] )
{

   int ix, iy; // screen coordinate = indices of virtual 2D array
  //unsigned int i; // index of 1D array
  long double zx = creall(point);
  long double zy = cimagl(point);
  int r;
  
  if (zx< ZxMax && zx > ZxMin && zy > ZyMin && zy < ZyMax){
  
  	ix = (int) ((zx- ZxMin)/PixelWidth); 
  	iy = (int) ((ZyMax - zy)/PixelHeight); // inverse Y axis 
  	r = iDrawPoint(ix, iy, iColor, A);
  	if (r==0) return 0; // succes
        
  	}
  	
  printf(" bad point from dDrawPoint : z = %Le;%Le  \n", zx, zy); 
  point_errors += 1; // not draw 
  return 1; // error 
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













int FillArray(unsigned char A[]){

	int i;
	for (i=0; i<iSize; i++){
	
		A[i] = 255; // white
	}
	return 0;



}





// save data array to pgm file 
int SaveArray2PGMFile( unsigned char A[], long double k, char* comment )
{
  
  FILE * fp;
  const unsigned int MaxColorComponentValue=255; /* color component is coded from 0 to 255 ;  it is 8 bit color file */
  char name [100]; /* name of file */
  snprintf(name, sizeof name, "%.0Lf", k); /*  */
  char *filename =strncat(name,".pgm", 4);
  
  
  
  /* save image to the pgm file  */      
  fp= fopen(filename,"wb"); /*create new file,give it a name and open it in binary mode  */
  fprintf(fp,"P5\n # %s\n %u %u\n %u\n", comment, iWidth, iHeight, MaxColorComponentValue);  /*write header to the file*/
  fwrite(A,iSize,1,fp);  /*write image data bytes to the file in one step */
  
  //
  printf("File %s saved. ", filename);
  if (comment == NULL)  printf ("empty comment \n");
                   else printf (" comment = %s \n", comment); 
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

complex long double GivePeriodic(complex long double c, complex long double z0, int period, long double eps2){

complex long double z = z0;
complex long double zPrev = z0; // prevoiuus value of z
int n ; // iteration




for (n=0; n<nMax; n++) {
     
    z = N( c, z, period);
   
    
  //  if (isnan(creall(z)) || isnan(cimagl(z))) 
   //	{ // printf("nan\n"); 
   //	  nan_errors+=1; 
   //	  break; }
    
    cDrawLine(zPrev, z, 100, data);
    if (cabs2(z - zPrev)< eps2) break; // succes 
    
    zPrev = z; }
    
    if ( n>maximal_n) maximal_n = n; // 
     
  // printf(" n = %d\n", n); 

return z;
}























int setup()
{	
	//ER2  = ER * ER; // ER*ER
	EPS2 = EPS * EPS;
	degree = (int) pow(2.0, period);// degree = 2^period
	sMax = 400; //4*degree; // number of starting points
	nMax = 10* degree;
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
	
	
	
	//  dynamic 1D array 
	data = malloc( iSize * sizeof(unsigned char) );
  	if (data == NULL )
    		{
      			fprintf(stderr," Could not allocate memory\n");
      			return 1;
    		}
	
	
	FillArray(data); // make all points white
	return 0;
}




int info(){
	printf (" period = %d  \n",  period); 
	printf (" degree of polynomial = %d  \n", degree); 
	printf (" number of starting points sMax = %d\n", sMax);
	printf (" dt = %Le\n", dt);
	printf (" radius of the circle = %Lf\n", radius);
	//printf (" escape radius ER = %Lf\n", ER);
	printf (" maximal allowed number of Newton iterations nMax = %d ( 10*degree) \n", nMax);
	printf (" maximal used number of Newton iterations maximal_n = %d  \n", maximal_n);
	
	printf ("stopping criterion epsilon_stop = %Le\n", EPS);
	printf (" c = %Lf ; %Lf \n",  creall(c), cimagl(c)); 
  	
  	// printf ("nan_errors = %d\n", nan_errors);
  	printf ("point_errors = %d\n", point_errors);
  	printf ("point_drawn = %d\n", point_drawn);
	return 0;
}


int end(){
	
        free(data);// free memory
        info();
	return 0;
}


/*===================== main ======================================================= */ 

int main() {

 
 
 
  
  int s; // number of starting point z0
  long double t= 0.0; //
  int r; // result of dDrawPoint  
  
  complex long double z0;
  
  setup();
  
  
  // create image 
  for (s=0; s<sMax; s++ ){
  		z0 = Give_z0(t, radius); // compute new starting point 
  		
  		//dDrawPoint(z0, 0,data); // draw 
  		zp = GivePeriodic( c , z0, period, EPS2); // compute periodic point 
  		//printf ("s = %d  t = %f z0 = %Lf ; %Lf  zp = %Lf;%Lf \n", s, t,  creall(z0), cimagl(z0), creall(zp), cimagl(zp));
  		r = dDrawPoint(zp, 0,data); // 
  		if ( r==0) point_drawn +=1;
  		//r=1;
  		t += dt;
  	}
  	
  	
  // 	
  SaveArray2PGMFile( data, period, "");
  end();
  
 
  
  
  
  return 0;
}
