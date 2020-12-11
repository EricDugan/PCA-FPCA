%%%  Analyses of headimpact impact data
%  D2x(t) = -beta0 x(t) - beta1 Dx(t) + alpha u(t)

%  add path to fda functions

addpath('../fdaM')
addpath('Examples/HeadImpact')

%% Set up the full data

load motorcycledata.txt

motot = motorcycledata(:,2);  %  time in milliseconds
motoy = motorcycledata(:,3);  %  deformation in 0.1 millimeters

%  adjust the data for baseline and plot

impact  = 14.0;  %  impact time
baseind = motot < impact;
basey   = mean(motoy(baseind));

%  remove baseline, change time, and convert to centimeters

motoy   = (basey - motoy);  

%% plot the data

figure(1)
plot(motot, motoy, 'bo', [0,60], [0,0], 'b:', 'LineWidth', 2)
% axis([0,60,-100,150])
xlabel('\fontsize{16} Time (milliseconds)')
ylabel('\fontsize{16} Acceleration (cm/msec^2)')

%% set up a pulse function located at times 14-15

motorng = [0,60];
Ubasis = create_bspline_basis(motorng, 3, 1, [0,14,15,60]);
Ufd = fd([0;1;0],Ubasis);

%% Set up the cell array coefCell

conbasis = create_constant_basis(motorng);
confdPar = fdPar(conbasis, 0, 0, 1);

fourierbasis = create_fourier_basis(motorng, 3, motorng(2));

% coefStr1.fun      = confdPar;
% coefStr1.parvec   = 0.073;
coefStr1.fun      = fourierbasis;
coefStr1.parvec   = [0.073; 0; 0];
coefStr1.estimate = 1;
coefStr1.coeftype = 'beta';

% coefStr2.fun      = confdPar;
% coefStr2.parvec   = 0.01;
coefStr2.fun      = fourierbasis;
coefStr2.parvec   = [0.01; 0; 0];
coefStr2.estimate = 1;
coefStr2.coeftype = 'beta';

coefStr3.fun      = confdPar;
coefStr3.parvec   = 0.25;
coefStr3.estimate = 1;
coefStr3.coeftype = 'alpha';

coefCell  = cell(3,1);
coefCell{1} = coefStr1;
coefCell{2} = coefStr2;
coefCell{3} = coefStr3;

[coefCell, ntheta] = coefCheck(coefCell);

disp(['ntheta = ',num2str(ntheta)])

%% construct basis for output x(t), multiple knots at times 14 and 15
%  Order 6 spline basis, with three knots at the impact point and 
%  three knots at impact + delta to permit discontinuous first 
%  derivatives at these points

knots     = [0,14,14,14,15,15,linspace(15,60,11)];
norder    = 6;
delta     = 1;
nbasis    = 21;
motobasis = create_bspline_basis(motorng,nbasis,norder,knots);

mototfine = linspace(0,60,501)';
motoyfine = eval_basis(mototfine, motobasis);

figure(2)
subplot(1,1,1)
plot(mototfine, motoyfine, 'b-', ...
     [14,14], [0,1], 'b--', [15,15], [0,1], 'b--', 'LineWidth', 2)
axis([0,60,0,1])
xlabel('\fontsize{16} Time t')
ylabel('\fontsize{16} Basis functions \phi_k(t)')

XbasisCell = cell(1);
XbasisCell{1} = motobasis;

%% Set up the cell array modelCell

modelCell = cell(1);

modelStr.order = 2;
modelStr.XCell = cell(2,1);
modelStr.FCell = cell(1);

Xterm0Str.variable   = 1;
Xterm0Str.derivative = 0;
Xterm0Str.ncoef      = 1;
Xterm0Str.factor     = -1;
modelStr.XCell{1}    = Xterm0Str;

Xterm1Str.variable   = 1;
Xterm1Str.derivative = 1;
Xterm1Str.ncoef      = 2;
Xterm1Str.factor     = -1;
modelStr.XCell{2}    = Xterm1Str;

FtermStr.ncoef = 3;
FtermStr.Ufd   = Ufd;
modelStr.FCell{1} = FtermStr;

modelCell{1} = modelStr;

modelCell = make_Model(XbasisCell, modelCell, coefCell);

%%  Define the cell array containing the data

yCell = cell(1);

yStr.argvals = motot;
yStr.y       = motoy;

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
conv     = [1e-8, 1e-4];  %  convergence criterion

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

disp(['Stiffness = ',num2str(thesave(9,1:3))]) 
disp(['Damping   = ',num2str(thesave(9,4:6))])
disp(['Forcing   = ',num2str(thesave(9,7))])

% display degrees of freedom and gcv values

disp('    rho      df         gcv')
disp([rhoVec', dfesave, gcvsave])

%% Evaluate the fit for parameter values at highest rho value

irho = 9;

theta = thesave(irho,:);
coefCell = BAwtvec2cell(theta, coefCell);
rho = rhoVec(irho);
[MSE, ~, ~, XfdCell, df, gcv, ISE, Var_theta] = ...
                     Data2LD(yCell, XbasisCell, modelCell, coefCell, rho);
disp(['MSE = ', num2str(MSE)])
disp(['df  = ', num2str(df)])
disp(['gcv = ', num2str(gcv)])

motofd   = getfd(XfdCell{1});
motofine = eval_fd(mototfine, motofd);
figure(3)
plot(mototfine, motofine, 'b-', 'LineWidth', 2)
hold on
plot(motot, motoy, 'bo', ...
          [impact,      impact],       [0,1], 'b--', ...
          [impact+delta,impact+delta], [0,1], 'b--', ...
          [impact,      impact+delta], [1,1], 'b--', ...
          [0,60], [0,0], 'b:', 'LineWidth', 2)
hold off
axis([0,60,-100,150])
xlabel('\fontsize{16} Time (msec)')
ylabel('\fontsize{16} Acceleration (cm/msec^2)')

% compute standard error and confidence limits for forcing

stderr = sqrt(diag(Var_theta));
theta = thesave(irho,:)';
thetaUP = theta + 2*stderr;
thetaDN = theta - 2*stderr;

disp(['rate constant       beta0:  ', ...
      num2str([theta(1), thetaDN(1), thetaUP(1), stderr(1)])])
disp(['rate constant       beta1:  ', ...
      num2str([theta(2), thetaDN(2), thetaUP(2), stderr(2)])])
disp(['forcing coefficient alpha: ', ...
      num2str([theta(3), thetaDN(3), thetaUP(3), stderr(3)])])
 
%%  plot the evolution of the parameters over the values of rho

figure(4)
subplot(3,1,1)
plot(rhoVec, thesave(:,1), 'bo-', 'LineWidth', 2)
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} parameter \beta_0')
subplot(3,1,2)
plot(rhoVec, thesave(:,2), 'bo-', 'LineWidth', 2)
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} parameter \beta_1')
subplot(3,1,3)
plot(rhoVec, thesave(:,3), 'bo-', 'LineWidth', 2)
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} parameter \alpha')

