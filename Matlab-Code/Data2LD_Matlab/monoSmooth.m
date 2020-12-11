function [yfd, bfd, df, gcv, MSE, theta] = ...
    monoSmooth(argvals, y, Xbasis, betafdPar, ...
               gammavec, dbglev, iterlim, conv)

if nargin < 8,  conv     = [1e-4, 1e-3]; end 
if nargin < 7,  iterlim  = 50;           end 
if nargin < 6,  dbglev   =  2;           end 
if nargin < 5,  gammavec = 0:1:8;        end 

[n, ncurve] = size(y);  
nXbasis = getnbasis(Xbasis);
if isa_basis(betafdPar)
    betabasis = betafdPar;
end
if isa_fd(betafdPar)
    betabasis = getbasis(betafdPar);
end
if isa_fdPar(betafdPar)
    betabasis = getbasis(getfd(betafdPar));
end
betanbasis = getnbasis(betabasis);

%  basis for data

XbasisCell = cell(1);
XbasisCell{1} = Xbasis;

%  damping coefficient 
 
coefStr1 = make_Coef(betafdPar, zeros(betanbasis,1), true);

coefCell    = cell(1);
coefCell{1} = coefStr1;

[coefCell, ntheta] = coefcheck(coefCell);

%  The damping coefficient is estimated but the reaction coefficient is 0.

modelStr.order = 2;

XtermStr = make_Xterm(1, 1, 1, 1);
XCell    = cell(1);
XCell{1} = XtermStr;

VariableStr     = make_Variable('f(x)', modelStr.order, XCell, []);
VariableCell    = cell(1);
VariableCell{1} = VariableStr;

ModelCell = make_Model(XbasisCell, VariableCell, coefCell);

DataCell = cell(1);
DataStr.argvals = argvals;

rhoVec   = exp(gammavec)./(1+exp(gammavec));
nrho     = length(rhoVec);
dfesave  = zeros(nrho,ncurve);
gcvsave  = zeros(nrho,ncurve);
MSEsave  = zeros(nrho,ncurve);
thesave  = zeros(nrho,ncurve,ntheta);

ycoef = zeros(nXbasis,ncurve);
bcoef = zeros(betanbasis,ncurve);
for icurve = 1:ncurve
    DataStr.y = y(:,icurve);
    DataCell{1} = DataStr;
    coefCell_opti = coefCell;
    for irho = 1:nrho
        rhoi = rhoVec(irho);
        if dbglev == 2
            disp(['rho = ',num2str(rhoi)])
        end
        theta_opti = Data2LD_opt(DataCell, XbasisCell, ...
                                 ModelCell, coefCell_opti, ...
                                 rhoi, conv, iterlim, dbglev);
        coefCell_opti = BAwtvec2cell(theta_opti, coefCell);
        [MSE, DpMSE, D2ppMSE, ~, df, gcv] = ...
            Data2LD(DataCell, XbasisCell, ModelCell, ...
                    coefCell_opti, rhoi);
        thesave(irho,icurve,:) = theta_opti';
        dfesave(irho,icurve)   = df;
        gcvsave(irho,icurve)   = gcv;
        MSEsave(irho,icurve)   = MSE;
    end
    
    theta    = thesave(nrho,icurve,:);
    coefCell = BAwtvec2cell(theta, coefCell);
    rho      = rhoVec(nrho);
    [MSE, ~, ~, XfdCell, df, gcv] = ...
        Data2LD(DataCell, XbasisCell, ModelCell, coefCell, rho);
    bcoef(:,icurve) = theta;
    ycoef(:,icurve) = getcoef(getfd(XfdCell{1}));
end
yfd = fd(ycoef,Xbasis);
bfd = fd(bcoef,betabasis);