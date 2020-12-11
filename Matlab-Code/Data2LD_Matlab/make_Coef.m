function coefStr = make_Coef(funobj, parvec, estimate)
% make_Coef assembles the three arguments into a struct object that is used
% by function make_Model to set up a linear dynamicsystem object.
%  Arguments are:
%  FUNOBJ   ... Either a functional data object that defines the coefficient
%               and its first derivative, or a user-supplied struct object with 
%               fields fd and difdip that compute the value of the function and 
%               its first derivative, respectively, and the field more that
%               is used for any additional information required by the 
%               user-supplied coefficient.
%  PARVEC   ... A vector of values for either the coefficients for the 
%               functional data object or for the user-supplied struct
%               object.
%  ESTIMATE ... If nonzero, estimate the values in PARVEC within function
%               Data2LD_Opt.  Otherwise, keep them fixed.
%               This may also be a logical vector of the same length as
%               PARVEC specifying which coefficients are fixed and which
%               are not.  

%  Last modified 1 February 2019

if nargin < 3, estimate = 1; end

%  check class of funobj

if ~isstruct(funobj) && ~isa_basis(funobj)  && ...
   ~isa_fd(funobj)   && ~isa_fdPar(funobj)
    error(['Argument FUNOBJ is not a struct object, basis object, ', ...
           'fd object or a fdPar object.']);
end

%  check class of parvec

if ~isnumeric(parvec)
    error('Argument PARVEC is not double.');
end

%  check dimension of parvec against number of basis functions

parvecdim = size(parvec);

if ~isstruct(funobj)
    if isa_basis(funobj)
        nfunbasis = getnbasis(funobj);
    end
    if isa_fd(funobj) || isa_fdPar(funobj)
        nfunbasis = getnbasis(getbasis(funobj));
    end
    if parvecdim(1) ~= nfunbasis
        error(['Length of parameter vector not equal to ', ...
               'the number of basis functions.']);
    end
end

if parvecdim(2) ~= 1
    error('Parameter vector is not a column vector.');
end

%  check the length of estimate and, if scalar, convert to vector

if length(estimate) > 1
    if length(estimate) ~= parvecdim(1)
        error(['Length of parameter vector not equal to ', ...
               'the length of estimate vector.']);
    end
else
    if estimate
        estimate = ones( parvecdim(1),1, 'logical');
    else
        estimate = zeros(parvecdim(1),1, 'logical');
    end
end

coefStr.fun      = funobj;
coefStr.parvec   = parvec;
coefStr.estimate = estimate;

end


