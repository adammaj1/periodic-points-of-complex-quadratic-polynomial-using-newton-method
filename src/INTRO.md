[README](README.md)


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
Newton method = Newton-Raphson method = normal Newton-Raphson method

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

## Schroeder method = Newton first type method 


[Modified for multiple roots](https://en.wikipedia.org/wiki/Newton%27s_method#Slow_convergence_for_roots_of_multiplicity_greater_than_1)


```math
z_{n+1} = z_n - m \frac {F(z_n)}{ F'(z_n)} 
```

where:
* m is a multiplicity of the root ( positive integer  )


# Number of periodic points and cycles

Number of periodic points d for period p and it's divisors = degree of function F defining periodic points 

$`d = 2^p`$


Number e of periodic points for exact period p  

$` e \le d `$ 

Number of cycles  

$`r = \frac{e}{p} `$





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

[Stability](https://en.wikipedia.org/wiki/Periodic_points_of_complex_quadratic_mappings#Stability_of_periodic_points_(orbit)_-_multiplier) of periodic point / cycle = magnitude of the multiplier


# Multiple roots and multiplicity of polynomial root
* [math.stackexchange question: estimating-the-multiplicity-of-a-root-numerically](https://math.stackexchange.com/questions/698858/estimating-the-multiplicity-of-a-root-numerically)
* [Dynamical control of Newton’s method for multiple roots of polynomials by S. Graillat, F. Jézéquel and M. S. Ibrahim](https://hal.archives-ouvertes.fr/hal-01363961/document)
* [math.stackexchange question: is-modified-newtons-raphson-method-redundant](https://math.stackexchange.com/questions/3089817/is-modified-newtons-raphson-method-redundant)

## How to compute (aproximate) multiplicity of the root numerically?

Multiplicity estimates
* [A study of accelerated Newton methods for multiple polynomial roots by Aurél Galántai and Aurél Galántai](https://www.researchgate.net/publication/225620593_A_study_of_accelerated_Newton_methods_for_multiple_polynomial_roots)
* [Multiplicity estimating algorithm for zeros of a complex polynomial and its applications by Tsuyako Miyakoda ](https://core.ac.uk/download/pdf/82716417.pdf)
* [Nonlinear Polynomial Systems: Multiple Roots and their Multiplicities by K. H. Ko, T Sakkalis](http://deslab.mit.edu/DesignLab/new_deslab/pubs/smi04.pdf)


### graphical method

[math.stackexchange question : what-is-the-intuition-for-the-multiplicity-of-a-root-of-a-polynomial-equation](https://math.stackexchange.com/questions/2716268/what-is-the-intuition-for-the-multiplicity-of-a-root-of-a-polynomial-equation)



### definition

[Math Tutor - multiplicity](http://math.feld.cvut.cz/mt/txtb/3/txe4ba3c.htm)

Let c be a root of a function f. We say that it is a root of multiplicity k if f (c) = 0, f ′(c) = 0, ... f (k−1)(c) = 0, but f (k)(c) ≠ 0.



###  derivative ratios
[ROOTS OF POLYNOMIALS BY RATIO OF SUCCESSIVE DERIVATIVES by James E. Crouse and Charles W. Putt](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19720015971.pdf)

the ratios of successive polynomial derivatives give the following very useful information: 
* the multiplicity of a root; that is, the number of roots at the root point approached
* the closeness of a trial point to the root approached
*  a good approximation as to where the next nearest root is when at a root.


The first derivative ratio  = $` P'(x)/P(x) `$
The second derivative ratio = $` P''(x)/P'(x) `$


### SE 

[SE](https://math.stackexchange.com/questions/698858/estimating-the-multiplicity-of-a-root-numerically)

if $`f(\cdot)`$ has a root at $`x`$, it holds $`f(x)=0`$. Furthermore, if you want to calculate the multiplicity you have to find the minimum $`m`$ s.t.:
$`f^{(m)}=0 `$

So, you can compute the derivatives, and if $`|f^{(m)}|<\varepsilon`$, where $`\varepsilon`$ represents a tolerance variable, thus $`m`$ is the number you are looking for.



### Grobner
 
  
 In the univariate case, where a is a root of f(x) of multiplicity d iff 
 
 ∂n ∂Xn (f)(α)=0 ∀n, 0 ≤ n<d.  

see : [ON MULTIPLICITIES IN POLYNOMIAL SYSTEM SOLVING by M. G. MARINARI, H. M. MOLLER, AND T. MORA](https://www.ams.org/journals/tran/1996-348-08/S0002-9947-96-01671-6/S0002-9947-96-01671-6.pdf)

### Ostrowski method
[Ostrowski 1973](https://books.google.nl/books?id=L_dqjIjOGBcC&pg=PA349&lpg=PA349&dq=Estimating%20the%20multiplicity%20of%20a%20root%20newton&source=bl&ots=yEtuRQ0PPR&sig=xNYnTGAzhctwCXzLvAWNvsGmDHs&hl=en&sa=X&ei=SaoVU-iJDof-ygPE24HYDg&ved=0CDcQ6AEwAjgU#v=onepage&q=Estimating%20the%20multiplicity%20of%20a%20root%20newton&f=false)


```math
m(x_1) = \lbrack  \frac{1}{2} + \frac{x_1 - x_2}{x_3 - 2x_2 + x_1}   \rbrack

```




where:
* m is a multiplicity of root x1
* x1 is a root
* x2 is the result of first Newton iteration
* x3 is the result of second Newton iteration

```math
\begin{matrix}
x_2 = N(x_1)\\
x_3 = N(x_2)
\end{matrix}
```

 

### Sanyasiraju VSS Yedida method

[Sanyasiraju VSS Yedida](https://mat.iitm.ac.in/home/sryedida/public_html/) [method](https://mat.iitm.ac.in/home/sryedida/public_html/caimna/transcendental/iteration%20methods/accelerating%20the%20convergence/mrac.html)


```math
m(x_1) =  \frac{1}{1-c}\\
```

where :   


```math
c =  \frac{x_3 - x_2}{x_2- x_1}
```


# Compare
* multipiplicity, multiplier and multiple root





# Newton Fractals by Chris Harshaw: 
* [blog](http://www.chrisharshaw.com/?p=34)
* [videos ](https://vimeo.com/user47350684)
* [python code](https://github.com/crharshaw/blog_materials)
  
# root finder  

online
* [Polynomial Root finder by Henrik Vestermark](http://www.hvks.com/Numerical/websolver.php)

offline



# See also
* [LA MÉTHODE DE NEWTON ET SON FRACTAL](http://images.math.cnrs.fr/La-methode-de-Newton-et-son-fractal.html?lang=fr) by Tan Lei
* [Newton Fractal](https://www.mitchr.me/SS/newton/index.html) by Mitch Richling
* [Newton's method for periodic points](http://mathr.co.uk/blog/2018-11-17_newtons_method_for_periodic_points.html) by Claude Heiland-Allen
* [Newton Fractals](https://cindyjs.org/gallery/main/NewtonFractal/) - interactive version ( using CindyJS)
* [Finding Roots in the Complex Plane by Ricky Reusser](https://observablehq.com/@rreusser/finding-roots-in-the-complex-plane)


