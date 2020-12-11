                
function [coefStr, iv, factor, Bvec, nWbasis, estim, ...
    nderiv, funtype] = getHomoTerm(modelStruct, coefCell)
%  get details for homogeneous term 
%  last modified 2 December 2018
iv      = modelStruct.variable;   % number of variable in this term
ncoef   = modelStruct.ncoef;      % index of coefficient object in coefList
factor  = modelStruct.factor;     % fixed scale factor
nderiv  = modelStruct.derivative; % order of derivative in this term
coefStr = coefCell{ncoef};        % coefficient object itself
Bvec    = coefStr.parvec;         % parameter vector for forcing coefficient
estim   = coefStr.estimate;       % parameter to be estimated? (TRUE/FALSE)
nWbasis = length(Bvec);           % number of basis functions in coefficient
funtype = isstruct(coefStr.fun);  % is a user-defined coefficient (TRUE/FALSE)
