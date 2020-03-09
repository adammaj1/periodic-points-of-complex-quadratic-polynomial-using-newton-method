
Introduction



Constant values:
* c is a parameter of f function
* p is a period


Roots $`\{ z : F(z) =  0 \}`$ are a periodic points of function $`f^p`$

Functions and derivatives:
* $`f`$ 
* $`F`$ is a function for computing periodic points of 
* $`N`$ is a Newton function = function used for Newton iteration [(Newton method) ](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/periodic_points#Newton_method)


Iterations: 
* quadratic iteration : $`z_{k+1} = f(z_k) `$
* Newton iteration : $`z_{n+1} = N(z_n) `$



Function f is [the complex quadratic polynomial](https://en.wikipedia.org/wiki/Complex_quadratic_polynomial) 

$`f(z) = z^2 + c`$

[First derivative of function f with respect to z](https://en.wikipedia.org/wiki/Complex_quadratic_polynomial#First_derivative_with_respect_to_z) is denoted by 

$`z' = f'(z) = d`$


[Iterated function](https://en.wikipedia.org/wiki/Complex_quadratic_polynomial#Notation)

$`z_1 = f^1(z_0) = f(z_0) `$  

$`z_n = f^n(z_0) =  f^1(f^{n-1}(z_0)) `$


[Periodic points of f ](https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/periodic_points) are roots of the equation : 

$`f^p(z) =  z `$

it can be converted to standard form: 

$`f^p(z) - z = 0`$

New function is intoduced:

$`F(z) = f^p(z) - z `$

It is used as a basic function in Newton method. 




# Iteration of quadratic polynomial
Both function and it's derivative are computed together by iteration : 


$` z'_0 = 1`$  
$`z_0`$   
$`z'_p = 2*z_{p-1}*z'_{p-1}`$  
$`z_p = z_{p-1}^2 + c`$  


# Newton method ( iteration)
Now one can iterate Newton function N 


```math
z_{n+1} = z_n - \frac {F(z_n)}{ F'(z_n)} =  N(z_n)
```

All of it is computed in c function: 
```c
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
```
## Number of periodic points and cycles

Number of periodic points d for period p and it's divisors = degree of function F defining periodic points 

$`d = 2^p`$


Number e of periodic points for exact period p  

$` e \le d `$ 

Number of cycles  

$`r = e / p `$





# Vieta's formula

One cech check if all roots are found using [Vieta's formulas](https://en.wikipedia.org/wiki/Vieta%27s_formulas)


It is computed using function: 
```c
long double ComputeVieteSum(){
	complex long double sum = 0;
	int d;
	for (d=0; d<distinc_points; d++ )
		sum += zzd[d]; // zzd is a array of roots (global variable)
	return cabsl(sum);
}
```



# Stability

[Stability](https://en.wikipedia.org/wiki/Periodic_points_of_complex_quadratic_mappings#Stability_of_periodic_points_(orbit)_-_multiplier) of periodic point / cycle
