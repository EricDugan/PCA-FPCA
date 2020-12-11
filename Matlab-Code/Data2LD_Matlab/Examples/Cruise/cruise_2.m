function DSvec = cruise_2(t, Svec, cruiseCell, coefCell, SetPtfd)
%  Cruise control equation ... A second order equation
%  Last modified 29 March 2020

DSvec = zeros(2,length(t));
SStr  = cruiseCell{1};

SStr1    = SStr.XCell{1};
ncoef1   = SStr1.ncoef;
coefStr1 = coefCell{ncoef1};
coef1    = coefStr1.parvec;
Wbasis1  = coefStr1.fun;
SfdC     = fd(coef1, Wbasis1);

SStr2    = SStr.XCell{2};
ncoef2   = SStr2.ncoef;
coefStr2 = coefCell{ncoef2};
coef2    = coefStr2.parvec;
Wbasis2  = coefStr2.fun;
SfdS     = fd(coef2, Wbasis2);

FfdC     = SfdC;

BTS0  = eval_fd(t,SfdC);
BTS1  = eval_fd(t,SfdS);
APC   = eval_fd(t,FfdC);
U11   = eval_fd(t,SetPtfd);
DSvec(1,:) = Svec(2,:);
DSvec(2,:) = -BTS0'.*Svec(1,:) - BTS1'.*Svec(2,:) + APC'.*U11';


