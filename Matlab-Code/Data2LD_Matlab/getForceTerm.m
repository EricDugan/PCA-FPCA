
function [coefStr, factor, Ufd, Avec, Ucoef, nUbasis, Aestim, ...
    nAbasis, funtype] = getForceTerm(modelStruct, coefCell)
%  get details for a forcing term
%  last modified 2 December 2018
ncoef   = modelStruct.ncoef;  % index of coefficient object in coefList
factor  = modelStruct.factor; % fixed scale factor
Ufd     = modelStruct.Ufd;    % functional data object for forcing term
coefStr = coefCell{ncoef};    % coefficient object itself
Avec    = coefStr.parvec;     % parameter vector for forcing coefficient
Aestim  = coefStr.estimate;   % parameter to be estimated? (TRUE/FALSE)
Ubasis  = getbasis(Ufd);      % functional basis object for forcing function
Ucoef   = getcoef(Ufd);       % coefficient vector for forcing function
nUbasis = getnbasis(Ubasis);  % number of coefficients for B-spline coeff.
nAbasis = length(Avec);       % number of basis fucntions
funtype = isstruct(coefStr.fun); % is a user-defined coefficient (TRUE/FALSE)
end
