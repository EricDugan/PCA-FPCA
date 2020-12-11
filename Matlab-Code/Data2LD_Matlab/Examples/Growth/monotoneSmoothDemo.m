%%%  Monotone smoothing of growth data
%  D2x(t) = beta Dx(t)

%  add path to fda functions

addpath('../fdaM')
addpath('Examples/Growth')

%  load growth data

load growth

%% plot the data

agerng = [1,18];

figure(1)
plot(age, hgtfmat, 'bo', 'LineWidth', 2)
xlabel('\fontsize{16} Age (years)')
ylabel('\fontsize{16} Height (cm)')

%% Set up the cell array coefCell for a second order equation 

%  damping coefficient is unconstrained

betanbasis = 7;
betabasis  = create_bspline_basis(agerng, betanbasis);
betafdPar  = fdPar(betabasis);

figure(2)
plot(betabasis)

coefStr1.fun      = betafdPar;
coefStr1.parvec   = zeros(betanbasis,1);
coefStr1.estimate = 1;
coefStr1.coeftype = 'beta';

coefCell    = cell(1);
coefCell{1} = coefStr1;

[coefCell, ntheta] = coefCheck(coefCell);

disp(['ntheta = ',num2str(ntheta)])

%% construct basis for output x(t), multiple knots at times 14 and 15
%  Order 6 spline basis, with three knots at the impact point and 
%  three knots at impact + delta to permit discontinuous first 
%  derivatives at these points

growthbreaks = age;
growthnbasis = 33;
growthbasis  = create_bspline_basis(agerng, growthnbasis);

growthtfine = linspace(1,18,501)';
growthyfine = eval_basis(growthtfine, growthbasis);

XbasisCell = cell(1);
XbasisCell{1} = growthbasis;

%% Set up the cell array modelCell as a second order linear equation

%  The damping coefficient is estimated but the reaction coefficient is 0.

modelCell = cell(1);

modelStr.order = 2;
modelStr.XCell = cell(1);

Xterm1Str.variable   = 1;
Xterm1Str.derivative = 1;
Xterm1Str.ncoef      = 1;
Xterm1Str.factor     = 1;
modelStr.XCell{1}    = Xterm1Str;

modelCell{1} = modelStr;

modelCell = make_Model(XbasisCell, modelCell, coefCell);

%%  Define the cell array containing the data

% igirl = 1;
% igirl = 2;
% igirl = 3;
igirl = 4;

yCell = cell(1);

yStr.argvals = age;
yStr.y       = hgtfmat(:,igirl);

yCell{1} = yStr;

%% An evaluation of the criterion at the initial values

rho = 0.5;  %  light smoothing;

%  this command causes Data2LD to set up and save the tensors

clear functions

[MSE, DMSE] = Data2LD(yCell, XbasisCell, modelCell, coefCell, rho);

%% Optimization of the criterion

%  algorithm constants

dbglev   =  1;    %  debugging level
iterlim  = 50;    %  maximum number of iterations
conv     = [1e-4, 1e-3];  %  convergence criterion

gammavec = 0:8;
rhoVec   = exp(gammavec)./(1+exp(gammavec));
nrho   = length(rhoVec);
dfesave = zeros(nrho,1);
gcvsave = zeros(nrho,1);
MSEsave = zeros(nrho,1);
thesave = zeros(nrho,ntheta);

coefCell_opt = coefCell;

for irho = 1:nrho
    rhoi = rhoVec(irho);
    theta_opti = Data2LD_opt(yCell, XbasisCell, modelCell, coefCell_opt, ...
                             rhoi, conv, iterlim, dbglev);
    coefCell_opti = BAwtvec2cell(theta_opti, coefCell);
    [MSE, DpMSE, D2ppMSE, XfdCell, df, gcv] = ...
              Data2LD(yCell, XbasisCell, modelCell, coefCell_opti, rhoi);
    thesave(irho,:) = theta_opti';
    dfesave(irho)   = df;
    gcvsave(irho)   = gcv;
    MSEsave(irho)   = MSE;
end

% display the optimal parameter values

disp(['Damping   = ',num2str(thesave(9,:))])

% display degrees of freedom and gcv values

disp('    rho      df         gcv')
disp([rhoVec', dfesave, gcvsave])

%% Evaluate the fit for parameter values at highest rho value

irho = 9;

disp(['Girl ',num2str(igirl)])

theta    = thesave(irho,:);
theta(betanbasis) = 0;
coefCell = BAwtvec2cell(theta, coefCell);
rho      = rhoVec(irho);
[MSE, ~, ~, XfdCell, df, gcv, ISE, Var_theta] = ...
                     Data2LD(yCell, XbasisCell, modelCell, coefCell, rho);
disp(['MSE = ', num2str(MSE)])
disp(['df  = ', num2str(df)])
disp(['gcv = ', num2str(gcv)])

betafd   = fd(theta',betabasis);
betafine = eval_fd(growthtfine, betafd);

growthfd     = getfd(XfdCell{1});
growthfine   = eval_fd(growthtfine, growthfd);
Dgrowthfine  = eval_fd(growthtfine, growthfd, 1);
D2growthfine = eval_fd(growthtfine, growthfd, 2);
D2fitfine    = Dgrowthfine.*betafine;

figure(3)
plot(growthtfine, growthfine, 'b-', age, hgtfmat(:,igirl), 'bo', ...
     'LineWidth', 2)
axis([1,18,70,170])
xlabel('\fontsize{16} Time (years)')
ylabel('\fontsize{16} height (cm)')
title(['\fontsize{16} Girl ',num2str(igirl), ...
       ',  Std. Error = ',num2str(round(sqrt(gcv),2)),' cm'])

figure(4)
plot(1:ntheta, theta', 'bo-')
xlabel('\fontsize{16} Parameter number')
ylabel('\fontsize{16} Parameter')

stderr = sqrt(diag(Var_theta));
figure(5)
plot(1:ntheta, theta', 'bo-', ...
     1:ntheta, theta' + 2*stderr, 'r--', ...
     1:ntheta, theta' - 2*stderr, 'r--')
xlabel('\fontsize{16} Parameter number')
ylabel('\fontsize{16} Parameter')
 
figure(6)
plot(growthtfine, D2growthfine, 'b-', growthtfine, D2fitfine, 'b--', ...
     [1,18], [0,0], 'b:', 'LineWidth', 2)
axis([1,18,-7,2])
xlabel('\fontsize{16} Time (years)')
ylabel('\fontsize{16} Height acceleration (cm/yr^2)')


%%  plot the evolution of the parameters over the values of rho

figure(7)
plot(rhoVec, thesave, 'bo-', 'LineWidth', 2)
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} parameter \beta_0')

