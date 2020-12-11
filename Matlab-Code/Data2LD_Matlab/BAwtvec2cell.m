function coefCellnew = BAwtvec2cell(thetavec, coefCell)
%  Places estimated weight coefficients only in THETAVEC into pre-existing 
%    cell array COEFCELL.

%  Last modified 6 February 2019

if nargin < 2
    error('Number of arguments is less than 2.');
end

if ~iscell(coefCell)
    error('coefCELL is not a cell array.');
end

thetavec = thetavec(:);
ncoef  = length(coefCell);
ntheta = length(thetavec);
coefCellnew = cell(ncoef,1);
m2 = 0;
for icoef=1:ncoef
    coefStructi = coefCell{icoef};
    if coefStructi.estimate
        coefCellnew{icoef} = coefStructi;
        m1 = m2 + 1;
        m2 = m2 + length(coefStructi.parvec);
        coefStructi.parvec = thetavec(m1:m2);
    end
    coefCellnew{icoef} = coefStructi;
end
