function XStr = make_Xterm(variable, ncoef, derivative, factor)
%  make_Xterm assembles four arguments into a struct object that
%  defines a single homogeneous term in a linear differentiable 
%  equation.
%  Arguments are:
%  VARIABLE   ... The index of the variable in the system
%  NCOEF      ... The index of the coefficient in cell array coefCell.
%  DERIVATIVE ... The order of the derivative of the variable in the term.
%                 Defaults to 0.
%  FACTOR     ... The fixed known constant multiplying the coefficient
%                 function.  Defaults to 1.

%  Last modified 26 January 2019

if nargin < 4, factor     = 1;  end
if nargin < 3, derivative = 0;  end

if nargin < 2
    error('Less than two arguments supplied.');
end

if floor(variable) ~= variable || variable < 1
    error('Argument VARIABLE is not a positive integer.');
end
if floor(ncoef) ~= ncoef || ncoef < 1
    error('Argument NCOEF is not a positive integer.');
end
if floor(derivative) ~= derivative || derivative < 0
    error('Argument DERIVATIVE is not a positive integer.');
end
XStr.variable   = variable;
XStr.ncoef      = ncoef;
XStr.derivative = derivative;
XStr.factor     = factor;
end
