/*

gcc m.c -Wall -lm
./a.out >2.txt

*/
#include <complex.h>
#include <math.h>
#include <stdio.h>


// parameter c of the function fc(z) = z^2+c is c = -0.7500000000000000 ; 0.0000000000000000 

const long double pi = 3.1415926535897932384626433832795029L;
long double EPS2 = 1e-18L*1e-18L; // 
complex double c = -0.75;
complex double z = -0.5;//1.5; //-0.5;






//https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
int sign(long double x){
	if (x > 0.0) return 1;
	if (x < 0.0) return -1;
	return 0;

}

int DifferentSign(long double x, long double y){
	if (sign(x)!=sign(y)) return 1;
	return 0;


}




long double complex Give_z0(long double InternalAngleInTurns, long double radius )
{
  
  //0 <= InternalAngleInTurns <=1
  long double a = InternalAngleInTurns *2.0*pi; // from turns to radians
  long double Cx, Cy; /* C = Cx+Cy*i */
  
  
      Cx = radius*cosl(a); 
      Cy = radius*sinl(a); 
      

  return Cx + Cy*I;
}
/*
multiplicity of the root

https://sites.math.washington.edu/~morrow/336_14/fta.pdf
https://stackoverflow.com/questions/60801319/how-to-estimate-multiplicity-of-the-polynomial-roots
https://math.stackexchange.com/questions/2716268/what-is-the-intuition-for-the-multiplicity-of-a-root-of-a-polynomial-equation)

*/

int GiveMultiplicity(complex long double zr, int pMax){

	int s; // number of starting point z0
	int sMax = 20*pMax; // it should be greater then 2*pMax
  	long double t= 0.0; // angle of circle around zr, measured using carg function in radians, range [-pi, pi] 
  	long double dt = 1.0 / sMax; // t step 
  	long double radius = 0.001; // radius should be smaller then minimal distance between roots 
 
 	int p; 
  	
  	long double arg_old = 0.0;
  	long double  arg_new = 0.0;
  	
  	int change = 0;
  	complex long double z;
  	complex long double z0;
  	//complex long double zp;
  	
  
  	printf("#t carg\n"); // for gnuplot
	//
	for (s=0; s<sMax; ++s){
  		z0 = zr + Give_z0(t, radius); // z =  point on the circle around root zr
  		
  		// compute zp = f^p(z)
  		z = z0;
 		for (p=0; p < pMax; ++p){z = z*z + c ;} /* complex quadratic polynomial */
  		// turn (zp-z0)	
  		z = z - z0; // equation for periodic_points of f for period p 
 		
 		arg_new = carg(z); 
 		printf("%.16Lf %.16Lf\n", t, arg_new);// for the  gnuplot 
 		if (DifferentSign(arg_new, arg_old)) {change+=1;}
 		arg_old = arg_new;
			
  		
  		
  		//printf("z0 = %.16f %.16f  zp = %.16f %.16f\n", creal(z0), cimag(z0), creal(zp), cimag(zp));
  		t += dt; // next angle using globl variable dt
  	}
  	//printf("change = %d\n",change);
  	return change/2;
  	
}



int main(){

	printf("multiplicity = %d\n", GiveMultiplicity(z,2));

	return 0;
}
