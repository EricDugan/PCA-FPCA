%%%                 Data2LD analyses of weather data
%  Montreal's average temperature, centered on January 1, is fit by
%  a first order linear equation with a time-varying coefficient and
%  forced by a constant and by a cosine function representing solar
%  heating translated by 10 days:
%         DT(t) = - \beta(t) + \alpha_1 + \alpha_2 U(t)
%  where U(t) = -cos((2*pi/365)*(t+192)).
%  Rate function \beta(t) is constrained to be positive by expressing it
%  as the exponential of a B-spline function with seven basis functions.
%  This illustrates the use of a user-defined pair of functions for 
%  \beta and its derivative with respect to the coefficients.
%  The user-defined functions are:
%
% function bval = fun_explinear(t, bvec, Bbasisobj)
% if isa_fdPar(Bbasisobj)
%     Bbasisobj = getbasis(Bbasisobj);
% end
% basismat = eval_basis(t,Bbasisobj);
% bval     = exp(basismat*bvec);
% end
%
% function Dbval = fun_Dexplinear(t, bvec, Bbasisobj)
% nbasis    = length(bvec);
% if isa_fdPar(Bbasisobj)
%     Bbasisobj = getbasis(Bbasisobj);
% end
% basismat  = eval_basis(t, Bbasisobj);
% bval      = exp(basismat*bvec);
% Dbval     = basismat.*repmat(bval,1,nbasis);
% end

addpath('../fdaM')
addpath('Examples/DailyTemp')

load daily

%%  set up five-day block averages with winter centering
%  Here we use 73 five-day block averages as the data with the block
%  centers as the time points in order to speed up computation

%  set up centers of 73 5-day blocks

daytime73 = linspace(2.5,362.5,73)';

%  set up block averages

tempav73 = zeros(73,35);
m2 = 0;
for i=1:73
    m1 = m2 + 1;
    m2 = m2 + 5;
    tempavi = mean(tempav(m1:m2,1:35));
    tempav73(i,:) = tempavi;
end

%  winter-center the data

winterind73  = [ 37:73,1: 36];
tempav73 = tempav73(winterind73,:);

station = 12;  %  Montreal

% plot the data

figure(1)
plot(daytime73, tempav73(:,station), 'bo', ...
     [0,365], [0,0], 'b:', 'LineWidth', 2)
axis([0,365,-15,25])
xlabel('\fontsize{16} Time (days)')
ylabel('\fontsize{16} Block average temperature (deg C)')

%%  Define the two forcing functions 

%  constant U for constant forcing

Uconbasis = create_constant_basis(dayrange);
Uconfd    = fd(1, Uconbasis);
Uconvec   = ones(73,1);

%  cosine U for radiative forcing

uvec = -cos((2*pi/365)*(daytime+10+182));
Ucosbasis = create_fourier_basis(dayrange, 3);
Ucosfd    = smooth_basis(daytime, uvec, Ucosbasis);

%%  plot the cosine forcing function

figure(2)
plot(daytime, uvec, 'b-', 'LineWidth', 2)
axis([0,365,-1,1])
xlabel('\fontsize{16} Time (days)')
ylabel('\fontsize{16} Solar forcing')

%% set up basis objects for homogeneous and forcing terms

%  homogeneous or dynamics term

nWbasis   = 7;
Wbasisobj = create_fourier_basis(dayrange, nWbasis);

%  constant forcing coefficient

nAbasisC   = 1;
AbasisobjC = create_constant_basis(dayrange);

%  cosine forcing coefficient

nAbasisF   = 1;
AbasisobjF = create_constant_basis(dayrange);

%%  set up coefCell

%  struct object for homogeneous term

linfun.fd   = @fun_explinear;  % user-defined function for \beta(t)
linfun.Dfd  = @fun_Dexplinear; % user-defined function for D\beta(t)
linfun.more = Wbasisobj;

% linfun  = Wbasisobj;
StrW.fun      = linfun;
StrW.parvec   = zeros(nWbasis,1);
StrW.estimate = 1;

%  constant forcing coefficient

StrA1.fun      = AbasisobjC;
StrA1.parvec   = 1;
StrA1.estimate = 1;

%  cosine forcing coefficient

StrA2.fun      = AbasisobjF;
StrA2.parvec   = 1;
StrA2.estimate = 1;

%  coefCell constructed

coefCell = cell(3,1);
coefCell{1} = StrW;
coefCell{2} = StrA1;
coefCell{3} = StrA2;

[coefCell, ntheta] = coefCheck(coefCell);

disp(['ntheta = ',num2str(ntheta)])

%%  basis for 5-day block averages

%  saturated order 5 bspline basis

norder    = 5;
daybreaks = 0:5:365;
nbreaks   = length(daybreaks);
nbasis    = norder + nbreaks/2 - 2;
daybasis  = create_bspline_basis(dayrange, nbasis);

XbasisCell = cell(1);
XbasisCell{1} = daybasis;

%%  define model Cell

%  struct objects for each term

tempXStr0.variable   =  1;
tempXStr0.derivative =  0;
tempXStr0.ncoef      =  1;
tempXStr0.factor     = -1;  % term required to be negative

tempFStr1.Ufd    = Uconfd;
tempFStr1.ncoef  = 2;

tempFStr2.Ufd    = Ucosfd;
tempFStr2.ncoef  = 3;

%  model struct object

tempStr.name = 'temperature';
tempStr.order = 1;
tempStr.XCell = {tempXStr0};
tempStr.Fcell = cell(2,1);
tempStr.FCell{1} = tempFStr1;
tempStr.FCell{2} = tempFStr2;

% set up dailyCell

dailyCell = {tempStr};

dailyCell = make_Model(XbasisCell, dailyCell, coefCell);

%%  set up yCell

yCell = cell(1);

yStruct.argvals = daytime73;
yStruct.y       = tempav73(:,station);
yCell{1}        = yStruct;

%%  preliminary analysis to set up four-way tensors

clear functions

[MSE, DpMSE, D2ppMSE] = ...
                    Data2LD(yCell, XbasisCell, dailyCell, coefCell, 0.5);

%%  set up analysis 

%  set constants for estimation algorithm Data2LD_Opt

dbglev  =  1;    
iterlim = 50;    
% convrg  = [1e-5, 1e-4];  
convrg  = 1e-5;  

%  define rhovec using the logit function

gammavec = 0:1:7;
rhoVec = exp(gammavec)./(1+exp(gammavec));
nrho   = length(rhoVec);

%  arrays to save the data

dfesave = zeros(nrho,1);
gcvsave = zeros(nrho,1);
SSEsave = zeros(nrho,1);
thesave = zeros(nrho,ntheta);

%  define initial state of coefCell_opt

coefCell_opt = coefCheck(coefCell);

%%  Analyis of Montreal's temperature

for irho = 1:nrho
    rhoVeci = rhoVec(irho);
    disp(['Rho = ',num2str(rhoVeci)])
    theta_opti = Data2LD_opt(yCell, XbasisCell, ...
                            dailyCell, coefCell_opt, rhoVeci, ...
                            convrg, iterlim, dbglev);
    coefCell_opt = BAwtvec2cell(theta_opti, coefCell);    
    [SSE, ~, ~, ~, df, gcv] = ...
        Data2LD(yCell, XbasisCell, dailyCell, coefCell_opt, rhoVeci);
    thesave(irho,:) = theta_opti;    
    SSEsave(irho)   = SSE;
    dfesave(irho)   = df;
    gcvsave(irho)   = gcv;
end

%  display the stored values

disp('    rho      df         gcv')
disp([rhoVec', dfesave, gcvsave])

%% Evaluate the fit for parameter values at highest rho value

irho = 8;

theta = thesave(irho,:);
coefCell = BAwtvec2cell(theta, coefCell);
rho = rhoVec(irho);
[MSE, ~, ~, XfdCell, df, gcv, ISE, Var_theta] = ...
                     Data2LD(yCell, XbasisCell, dailyCell, coefCell, rho);
disp(['MSE = ', num2str(MSE)])
disp(['df  = ', num2str(df)])
disp(['gcv = ', num2str(gcv)])

tempfd   = getfd(XfdCell{1});
tfine    = linspace(0,365,101)';
tempfine = eval_fd(tfine, tempfd);

figure(3)
plot(tfine, tempfine, 'b-', 'LineWidth', 2)
hold on
plot(daytime73, tempav73(:,station), 'bo', ...
     daytime, uvec, 'b--', ...
     [0,365], [0,0], 'b:', 'LineWidth', 2)
hold off
axis([0,365,-15,25])
xlabel('\fontsize{16} Time (days)')
ylabel('\fontsize{16} Temperature (deg C)')
legend('\fontsize{16} Data fit', '\fontsize{16} Data', ...
       '\fontsize{16} Solar forcing', 'Location', 'North')

%% plot the optimal parameter values as functions of rho

indW = 1:nWbasis;
indA1 = nWbasis+1:nWbasis+nAbasisC;
indA2 = nWbasis+nAbasisC+1:ntheta;

%  plot flow of coefficients for log \beta

figure(4)
subplot(1,1,1)
plot(rhoVec, thesave(:,indW), 'bo-', 'LineWidth', 2)
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} log \beta coefficients')

%  plot flow of alpha 1 and alpha 2

figure(5)
subplot(2,1,1)
phdl = plot(rhoVec, thesave(:,indA1), 'bo-');
set(phdl, 'LineWidth', 2)
ylabel('\fontsize{16} \alpha_1')
subplot(2,1,2)
phdl = plot(rhoVec, thesave(:,indA2), 'bo-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} \alpha2')

%%  plot \beta with confidence intervals

betafd     = fd(thesave(irho,indW)', Wbasisobj);
betavec    = exp(eval_fd(daytime73,betafd));
Var_thetaW = Var_theta(indW,indW);
basismatW  = diag(betavec)*eval_basis(daytime73, Wbasisobj);
CI_beta    = 2*sqrt(diag(basismatW*Var_thetaW*basismatW'));

figure(6)
subplot(1,1,1)
plot(daytime73, betavec , 'b-', ...
     daytime73, betavec  + CI_beta , 'b--', ...
     daytime73, betavec  - CI_beta , 'b--', 'LineWidth', 2)
axis([0,365,0,0.05])
xlabel('\fontsize{16} Time (days)')
ylabel('\fontsize{16} \beta(t)')

%  display 95% confidence limits for \alpha1 and \alpha2

alpha1 = thesave(irho,indA1);
alpha2 = thesave(irho,indA2);

alpha1_stderr = sqrt(Var_theta(indA1));
alpha2_stderr = sqrt(Var_theta(indA2));

disp('    alpha_1   lower CI upper CI')
disp([alpha1, alpha1-2*alpha1_stderr, alpha1+2*alpha1_stderr])
disp('    alpha_2   lower CI upper CI')
disp([alpha2, alpha2-2*alpha2_stderr, alpha2+2*alpha2_stderr])

%  plot fit to derivative, lower panel in Figure 12.10 in the book

dayvec  = eval_fd(daytime73, tempfd);
Ddayvec = eval_fd(daytime73, tempfd, 1);
Ucosvec = eval_fd(daytime73, Ucosfd);

Ddayfit = -betavec .* dayvec + alpha1 + alpha2.*Ucosvec;

figure(8)
phdl = plot(daytime73, Ddayfit, 'b-', daytime73, Ddayvec, 'b--');
set(phdl, 'LineWidth', 2)
V = axis();
axis([0,365,V(3),V(4)])
xlabel('\fontsize{16} Time (days)')
ylabel('\fontsize{16} DT(t)')
legend('\fontsize{16} DT fit', '\fontsize{16} actual DT', ...
       'Location', 'NorthWest')


