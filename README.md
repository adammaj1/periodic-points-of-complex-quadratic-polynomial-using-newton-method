[How to find numerically all roots of a polynomial, especially the complex ones, by using the Newton-Raphson Method](https://math.stackexchange.com/questions/998333/finding-the-all-roots-of-a-polynomial-by-using-newton-raphson-method)


# Introduction



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




## Iteration of quadratic polynomial
Both function and it's derivative are computed together by iteration : 


$` z'_0 = 1`$  
$`z_0`$   
$`z'_p = 2*z_{p-1}*z'_{p-1}`$  
$`z_p = z_{p-1}^2 + c`$  


## Newton method ( iteration)
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





## Vieta's formula

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





# Period 12
>>>
The Newton iteration for 400 starting points on a circle of radius r = 2 (here for a polynomial of degree 4096, so we do not have enough starting points to find all roots; the polynomial shown
here describes periodic points of period dividing 12 of f(z) = z^2 + i). The apparent lines connect orbits under the Newton dynamics; 
>>>
   
   
Quote from the paper: [Newton's method in practice II: "The iterated refinement Newton method and near-optimal complexity for finding all roots of some polynomials of very large degrees" by Marvin Randig, Dierk Schleicher, Robin Stoll](https://arxiv.org/abs/1703.05847)



Here is my image: 

![12](./images/12.png) 

Text output of c program 
```
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
```


Code:
* [m.c](./src/m.c) - c code ( 1 file program which creates pgm file)


Conversion from pgm to png using Image Magick: 

```bash
convert 12.pgm -resize 600x600 12.png
```

# Period 1
![1](./images/1_4_15.png)

``` 
 parameter c of the function fc(z) = z^2+c is c = 0.0000000000000000 ; 1.0000000000000000 

 period = 1  
 degree of polynomial = 2  = 2^period
prime factors of 1  = 
 number of roots = number of periodic points = degree of polynomial = 2  
 number of starting points sMax = 4
succes : all 2 distinct points are found !!

 dt = 2.500000e-01
 radius of the circle around all periodic points = 2.000000e+00
 
 maximal allowed number of Newton iterations nMax = 120  =  10*degree + 100, see setup
 maximal used number of Newton iterations maximal_n = 29 
stopping criterion for the Newton iteration is epsilon_stop = 1.000000e-18
m_dist = 1.000000000000000036e-10
 minimal distnce in zzd =1.414214e+00 between roots
 

 periodic points are: 
 z = +1.300242590220120419; -0.624810533843826587 exact period = 1 stability = 2.885147480892119399
 z = -0.300242590220120419; +0.624810533843826587 exact period = 1 stability = 1.386410929247594584





 attracting limit cycle with exact period : 



 the sum of all roots should be zero by Viete’s formula (this sum should be the negative of the degree d − 1 coefficient)
 Viete sum = 1.000000000000000000e+00 ( it should be zero )```
```


Summary:
* Number of the roots = 2
  * 2 fixed points 



# Period 2





Newton basins  

![2_8_12.png](./images/2_8_12.png) 



Newton basins, rays and periodic points  

![2_8_13.png](./images/2_8_13.png) 


Newton basins with levels sets   

![2_8_14.png](./images/2_8_14.png) 


Newton basins with level sets, rays and periodic points  

![2_8_15.png](./images/2_8_15.png) 



File names are p_sMax_n.png where:
* p is a period 
* sMax is a number of starting points 
* n is an arbitrary number of a picture


Text output of [the program](./src/n.c) 

```
parameter c of the function fc(z) = z^2+c is c = 0.0000000000000000 ; 1.0000000000000000 

 period = 2  
 degree of polynomial = 4  = 2^period
prime factors of 2  = 2	
 number of roots = number of periodic points = degree of polynomial = 4  
 number of starting points sMax = 8
succes : all 4 distinct points are found !!

 dt = 1.250000e-01
 radius of the circle around all periodic points = 2.000000e+00
 
 maximal allowed number of Newton iterations nMax = 140  =  10*degree + 100, see setup
 maximal used number of Newton iterations maximal_n = 42 
stopping criterion for the Newton iteration is epsilon_stop = 1.000000e-18
m_dist = 1.000000000000000036e-10
 minimal distnce in zzd =7.939947e-01 between roots
 

 periodic points are: 
 z = +1.300242590220120419; -0.624810533843826587 exact period = 1 stability = 2.885147480892119399
 z = -0.300242590220120419; +0.624810533843826587 exact period = 1 stability = 1.386410929247594584
 z = -1.000000000000000000; +1.000000000000000000 exact period = 2 stability = 5.656854249492380195
 z = -0.000000000000000000; -1.000000000000000000 exact period = 2 stability = 5.656854249492380195





 attracting limit cycle with exact period : 



 the sum of all roots should be zero by Viete’s formula (this sum should be the negative of the degree d − 1 coefficient)
 Viete sum = 1.045143821782786559e-19 ( it should be zero )```
```


Summary:
* Number of the roots = 4
  * 1 period 2 cycle (= 2 )
  * 2 fixed points 




## Period 3

![3](./images/3_16_15.png)

```
  parameter c of the function fc(z) = z^2+c is c = 0.0000000000000000 ; 1.0000000000000000 

 period = 3  
 degree of polynomial = 8  = 2^period
prime factors of 3  = 3	
 number of roots = number of periodic points = degree of polynomial = 8  
 number of starting points sMax = 16
succes : all 8 distinct points are found !!

 dt = 6.250000e-02
 radius of the circle around all periodic points = 2.000000e+00
 
 maximal allowed number of Newton iterations nMax = 180  =  10*degree + 100, see setup
 maximal used number of Newton iterations maximal_n = 58 
stopping criterion for the Newton iteration is epsilon_stop = 1.000000e-18
m_dist = 1.000000000000000036e-10
 minimal distnce in zzd =3.585710e-01 between roots
 

 periodic points are: 
 z = +1.300242590220120419; -0.624810533843826587 exact period = 1 stability = 2.885147480892119399
 z = +0.100102984663451134; -0.349424850666336804 exact period = 3 stability = 3.168659127802740935
 z = -0.300242590220120419; +0.624810533843826587 exact period = 1 stability = 1.386410929247594584
 z = -0.112077118724660643; +0.930043059065437968 exact period = 3 stability = 3.168659127802740858
 z = -0.852418811174176060; +0.791526907300152676 exact period = 3 stability = 3.168659127802740441
 z = -1.290491233241733357; +0.779281718235989916 exact period = 3 stability = 20.197817884052375560
 z = +0.096796551780185860; -1.140114382717044562 exact period = 3 stability = 20.197817884052375560
 z = +1.058087626696933066; -1.011312451218199194 exact period = 3 stability = 20.197817884052375560





 attracting limit cycle with exact period : 



 the sum of all roots should be zero by Viete’s formula (this sum should be the negative of the degree d − 1 coefficient)
 Viete sum = 1.084202172485504434e-19 ( it should be zero )```
```


Summary:
* Number of the roots = 8
  * 2 period 3 cycles (= 6 )
  * 2 fixed points 
  
  
  
  
Another example:
```
 parameter c of the function fc(z) = z^2+c is c = -0.1225611668766540 ; 0.7448617666197440 

 period = 3  
 degree of polynomial = 8  = 2^period
prime factors of 3  = 3	
 number of roots = number of periodic points = degree of polynomial = 8  
 number of starting points sMax = 16
succes : all 8 distinct points are found !!

 dt = 6.250000e-02
 radius of the circle around all periodic points = 2.000000e+00
 
 maximal allowed number of Newton iterations nMax = 180  =  10*degree + 100, see setup
 maximal used number of Newton iterations maximal_n = 65 
stopping criterion for the Newton iteration is epsilon_stop = 1.000000e-18
m_dist = 1.000000000000000036e-10
 minimal distnce in zzd =3.065014e-01 between roots
 

 periodic points are: 
 z = +1.276337623593117529; -0.479727984309394897 exact period = 1 stability = 2.727032576504271481
 z = +0.000000000000000385; +0.000000000000000798 exact period = 3 stability = 0.000000000000004646
 z = -0.276337623593117529; +0.479727984309394897 exact period = 1 stability = 1.107251409830028695
 z = -0.122561166876653999; +0.744861766619744015 exact period = 3 stability = 0.000000000000004646
 z = -0.662358978622372968; +0.562279512062300512 exact period = 3 stability = 0.000000000000004646
 z = -1.247255127827276373; +0.662949719009367346 exact period = 3 stability = 16.158407096723733116
 z = +0.038593416246120549; -1.061217891258985840 exact period = 3 stability = 16.158407096723733444
 z = +0.993581857080182406; -0.908873106432426830 exact period = 3 stability = 16.158407096723733444





 attracting limit cycle with exact period : 
 z = +0.000000000000000385; +0.000000000000000798 exact period = 3 stability = 0.000000000000004646
 z = -0.122561166876653999; +0.744861766619744015 exact period = 3 stability = 0.000000000000004646
 z = -0.662358978622372968; +0.562279512062300512 exact period = 3 stability = 0.000000000000004646



 the sum of all roots should be zero by Viete’s formula (this sum should be the negative of the degree d − 1 coefficient)
 Viete sum = 2.235140038340936236e-19 ( it should be zero )
```
   


# Files
* [m.c](./src/m.c) - c code ( 1 file program which creates 12.pgm file)
* [n.c](./src/n.c) - c code ( 1 file program which creates pgm files: basins and rays  )
* [p.mac](./src/p.mac) - Maxima CAS batch file ( program) 
* [2.mac](./src/2.mac) - Maxima CAS batch file ( program) for checking peroid 2 case 


# See also
* [LA MÉTHODE DE NEWTON ET SON FRACTAL](http://images.math.cnrs.fr/La-methode-de-Newton-et-son-fractal.html?lang=fr) by Tan Lei
* [Newton Fractal](https://www.mitchr.me/SS/newton/index.html) by Mitch Richling
* [Newton's method for periodic points](http://mathr.co.uk/blog/2018-11-17_newtons_method_for_periodic_points.html) by Claude Heiland-Allen
* [Newton Fractals](https://cindyjs.org/gallery/main/NewtonFractal/) - interactive version ( using CindyJS)
* Newton Fractals by Chris Harshaw: 
  * [blog](http://www.chrisharshaw.com/?p=34)
  * [videos ](https://vimeo.com/user47350684)
  * [python code](https://github.com/crharshaw/blog_materials)


# Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

# License

This project is licensed under the  GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007 - see the [LICENSE](LICENSE) file for details

# technical notes
GitLab uses:
* the Redcarpet Ruby library for [Markdown processing](https://gitlab.com/gitlab-org/gitlab-ce/blob/master/doc/user/markdown.md)
* KaTeX to render [math written with the LaTeX syntax](https://gitlab.com/gitlab-org/gitlab-ce/blob/master/doc/user/markdown.md), but [only subset](https://khan.github.io/KaTeX/function-support.html)




## Git
```
cd existing_folder
git init
git remote add origin git@gitlab.com:adammajewski/periodic-points-of-complex-quadratic-polynomial-using-newton-method.git
git add .
git commit -m "Initial commit"
git push -u origin master
```


```
  git clone git@gitlab.com:adammajewski/periodic-points-of-complex-quadratic-polynomial-using-newton-method.git
```

Subdirectory

```git
mkdir images
git add *.png
git mv  *.png ./images
git commit -m "move"
git push -u origin master
```
then link the images:

```txt
![](./images/n.png "description") 

```

```git
gitm mv -f 
```

to overwrite 






local repo : ~/c/julia/periodic/newton/pgm2

