
Dynamic Data Analysis is a book authored by James Ramsay and Giles Hooker and  published in 2017 by Springer.  It explores systems of linear and nonlinear differential equations as models for functional data.  

Functional data are observations distributed over a continuum (most often time) that can be represented by a smooth curve that can be differentiated up to some required order of derivative.   

Very often in the history of science and in current applications one or more differential equations provide a neat, compact and accurate summary of how the data behave.  A large literature exist when the data to be fit are located only at the start of the observation interval.  These are *initial value* systems.  Slightly more generally, data may also be found at the end of the interval, and such equation are called *boundary value* systems.  But the book treats the most general case where data are available for one or more of the equations over the whole interval, what one might call *distributed value systems*.

It is usual that these equations involve parameters that define the equations but that require estimation from the data in order to completely define the differential equation system.  The book details a procedure for estimating these parameters and the differential equation system in a single analysis.

Code package Data2LD  is a collection of code in Matlab that is specifically for linear differential equations.  

A linear differential equation is a sum of terms, with each term being the product of coefficient function and the value of a derivative of a variable in the equation system.  The product may also include a fixed constant that multiplies the product, such as, in many cases, -1 when the coefficient function is assumed to be non-negative or strictly positive.  The term with the highest order of derivative usually appears on the left side of the equation with no multiplying coefficient.  Terms on the left can involve any lower order derivative of the equation variable, or any derivative of another equation variable.

A critical feature of each term is that the coefficient function may be either a constant or a function of time.  If the later, the function may NOT be a function of any of the variables in the system, but on the other hand may be a function of any set of observations external to the system.  

A linear differential equation may also contain additional terms in the sum that are known functions multiplied by a coefficient function.  These known functions are not one of the differential equations in the system, and are often called “forcing functions” and the terms called “forcing terms.”  

That is, external information can define a linear equation either through its role in defining a coefficient function or as a forcing term, or both.  Again we stress that these external variables are not any of the differential equation variables.

Linear differential equations do a lot of the heavy lifting in science because, first they are often themselves interesting models for the data and, also, can provide useful approximations for more complicated nonlinear equation systems.  Linear equations are also stable and have solutions everywhere in the interval.

##  The package contains five examples of differential equation systems that illustrate various equation structures and data that are to be fit.  Each of the examples carries out one or more analyses, and extensive commentary on how to set up the analysis.  The book *Dynamic Data Analysis *also contains analysis instructions and illustrations.