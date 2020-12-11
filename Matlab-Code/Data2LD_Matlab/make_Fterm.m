function FStr = make_Fterm(ncoef, Ufd, factor)
%  make_Fterm assembles four arguments into a struct object that
%  defines a single forcing term in a linear differentiable 
%  equation.
%  Arguments are:
%  NCOEF    ... The index of the coefficient in cell array coefCell.
%  UFD      ... A functional data object that represents the forcing
%               variable.
%  FACTOR   ... The fixed known constant multiplying the coefficient
%                 function.  Defaults to 1.

%  Last modified 26 March 2020

if nargin < 3, factor = 1;  end
if floor(ncoef) ~= ncoef || ncoef < 1
    error('Argument NCOEF is not a positive integer.');
end
if ~isa_fd(Ufd)
    error('Argument UFD is not a functional data object.');
end
FStr.ncoef    = ncoef;
FStr.Ufd      = Ufd;
FStr.factor   = factor;
end
