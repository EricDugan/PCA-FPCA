%% The cruise control demonstration.

% Last modified 31 March 2020

% Add paths to the functional data analuysis functions and to the 
% demonstration code, which is assumed to be two folders inside the
% Data2LD code.

addpath('../../fdaM')
addpath('Examples/Cruise')

%% The cruise control problem
%
% A driver starts his car and (1) accelerates to the 60 km/h speed limit
% holding on most main streets of Ottawa, (2) reduces speed to 40 km/h in a
% school zone, (3) enters a controlled-access boulevard with an 80 km/h
% limits and finally (4) against enters a typical principal street.
%
% Operating on a snowy winter morning, the driver takes 8 seconds to reach
% 60 km/h, or 4/30 seconds per unit change in speed.
%
% The speed of the car under acceleration is modelled by a first order
% constant coefficient equation:
%
%% The feedback equations: 
%
% There are two variables in this version of the problem:
%
%  S(t) is the speed of the vehicle at time t in kilometres per hour
%  C(t) is a process that controls the speed of the vehicle, and can
%  be thought of as the amount of depression of the accelerator or 
%  the set speed in a cruise control mechanism.  It's units are 
%  arbitrary.
%  Also, we have a variable that is an input to the second controller
%  equation:  S_0(t).  This the driver-controlled desired speed. 
%  In process control jargon, it is the set point variable.
%  This input is changed by stepping up and down its level, with the 
%  level S_0(t) being the desired speed in km/hr at time t.
%
% The two differental equations specifying the feedback model are:
%
% In standard format:
%
% $$DS(t) = -\beta_{11} S(t) + \beta_{12} C(t)$$
% $$DC(t) =  \beta_{21} S(t) - \beta_{22} C(t) + \alpha S_0(t)$$
%
% This standard input format hides how the equations work, so here's
% an alternative formulation that reveals more clearly how a feedback
% system works:
%
% In feedback format:
%
% $$DS(t) = -\beta_{11} S(t) + \beta_{12} C(t)$$
% $$DC(t) = -\beta_{22} C(t) + \alpha (S_0(t) - S(t))$$
%
% In this format it becomes clear that \beta_{21} = -\alpha, so that
% we really only have four parameters instead of five.
% In this form, we see that the controller's input is the difference
% between the set point and the actual speed.  
% When actual speed exceeds the set point, the controller receives
% negative input, and changes downwards.
% When actual speed is below the set point, the controller receives
% positive input, and changes downwards.
% When actual and set speed are equal, the controller receives
% no input, and does not change.
%
% There are five parameters in this model, all of them constants,
% with the constant values that we will use to generate some data:
% $\beta_{11}$ =  1:   The rate of increase of speed when it responds to 
%                      an input from the controller
% $\beta_{12}$ =  1/4: Determines the impact of the controller input on
%                      speed.
% $\beta_{21}$ =  -1:  The impact of speed on the controller.
% $\beta_{22}$ =  0:   The rate of increase of the controller level.
%                      When this is 0, any input to the equation
%                      directly controls the rate of increase in the
%                      controller level, that is, instantaneously
%                      changes the controller level.
% $\alpha$     =  1    The impact of the set point on, in this case,
%                      the rate of change of the controller level.
%
% Inserting the proposed coefficient values, the standard equations are:
%
% $$DS(t) = -1*S(t) + (1/4)*C(t)$$
% $$DC(t) = -1*S(t) +     0*C(t) + 1*S_0$$
%
% Notice that in both equations there are constants installed in front
% of some of the terms.  These are considered standard in the field of
% process controll, and therefore are not viewed as being estimated from 
% the data.  Coefficients \beta_{11} and \beta_{22} are 
% multiplied by the fixed values -1, and this is usual in equations
% that link change with level.  It is assumed that 
% \beta_{11} and \beta_{22} are in fact non-negative. 
% The feedback form of the equation shows that 
% We will use a special case of these equations to generate some 
% simulated data and then estimate the parameters defining the 
% equation using these simulated data.
%
% These feedback version of the equations tell us that:
%
% when the speed $S(t)$ is less than the set point value S_0(t), 
%      the control level increases  and forces the speed to increase, and
% when the speed $S(t)$ is greater than the set point value S_0(t), 
%      the control level decreases and forces the speed to decrease.
% The value $\beta_{22}$ = 0 implies that the controller responds
% instantly to a change in the difference $S_0(t) - S(t)$.  That may
% an exaggeration, but control process usually react faster than what
% they control.
%
%% Where we go in these experiments:
%
% The cruise control is not so much about process control as it is about
% working with sets of linear differential equations.
% For that reason, we will work with the standard five-coefficient
% model, and will do various things to make it behave like a process
% control system.  
% First we will estimate four parameters, holding beta_{22} fixed at 0.
% Then we will set the restriction beta_{22} + alpha = 0, keep beta_{22} 
% fixed at 0, and thus estimate only three coefficient values.
% Finally, we will combine the two first order equation together into a 
% single second order equation, and solve that.
%% ------------------------------------------------------------------------
%  Experiment 1:  Four coefficients to be estimated, one held fixed at 0.
%  ------------------------------------------------------------------------
%
% In the first we estimate the four parameters 
%  $\beta_{11}$, $\beta_{12}$, $\beta_{21}$ and $\alpha_2$,
% but fix the value of \beta_22 to be 0 rather than estimating it.
% The steps to setting up the analysis are:
% 1. Define the time range, the number of observation poinmts, and the
%    times of observation.
% 2. Define the set point function $S_0(t)$, which is a known input to
%    the controller variable.
% 3. Define the five coefficients.  Each is expressed as a functional
%    parameter object with a set multiplying factor and a specificatiion
%    of whether or not it will be estimated.
% 4. Define the functional basis system that will be used to define the
%    variation in the S and C variables.  
% 5. Define the model in a form that the code can understand.  Each of
%    the two variables is defined in a struct object, which contains
%    the definition of the internal structure of the linear differential
%    equation.
% 6. We simulate data, which is ideal for demos because the user can 
%    change the data in various ways.  To do this, we have to solve the
%    equations using the five parameter values set above, and then
%    define the data by adding some random noise to the solutions.
% 7. Set up a cell array containing the data to be analyzed.
% 8. Evaluate the model at the initial values of the parameters.  This
%    is important as a check, but also sets up some big arrays that only
%    need to be calculated once.
% 9. Set up a sequence of values of rho, the value between 0 and 1 
%    that controls the emphasis on fitting the equation (1) and 
%    fitting the data (0).
% 10. Carrying out the analysis over these rho values, pausing at 
%    each to display results.

%%  1. Defining the problem
% Set up the time span, a set of observation times, and a fine time grid 
% for plotting:

T     = 80;     %  maximum time in seconds
rng   = [0,T];  %  range of times to be used
n     = 41;     %  number of observations
tobs  = linspace(0,T,n)';  %  times of the observations
nfine = 501;  % number of points in a fine mesh for plotting purposes
tfine = linspace(0,T,nfine)';  %  mesh values

%% 2.  Set up the set-point forcing function $S_0(t)$.
%
% The set-point function uses an order 1 step function B-spline basis.  The 
% knots are placed at the points where the set-point changes values.

steporder  = 1;  %  step function basis
stepnbasis = 4;  %  four basis functions
stepbreaks = [0,20,40,60,80];
stepbasis  = create_bspline_basis(rng, stepnbasis, steporder, stepbreaks);
stepcoef   = [60;40;80;60];  % target speed for each step
SetPtfd    = fd(stepcoef, stepbasis);  % define the set point function
Ufine      = eval_fd(tfine, SetPtfd);    

%% 3.  Set up the coefficient dictionary in cell array coefCell
%
%  The total number of coefficients defining the estimated coefficient 
%  functions is three because each coefficient function has only one 
%  coefficient, and only three of them are estimated.

%  Set up a constant basis over [0,T], 
%  all of our coefficients will constants

conbasis  = create_constant_basis(rng);
confdPar  = fdPar(conbasis);

%  Function make_Coef defines for each coefficient a struct object
%  whose fields are: (1) the coefficient function object, 
%  (2)the initial value of its parameter vector, and 
%  (3) whether the function is to be estimated or held fixed.
%  You can make as many as you want, but each coefficient that you 
%  need should be in this list.

% Arguments:           fun       parvec  estimate  
coefStrS_S = make_Coef(confdPar,      1,     true); % speed   beta_ss
coefStrS_C = make_Coef(confdPar,    1/4,     true); % speed   beta_sc
coefStrC_S = make_Coef(confdPar,      1,     true); % control beta_cs
coefStrC_C = make_Coef(confdPar,      0,    false); % control beta_cc
coefStrC_F = make_Coef(confdPar,      1,     true); % control alpha_c

%  Set up the cell array object coefCell for five coefficients 
%  and assign the coefficient definitions to the cells.
%  Although there are five in the full equation, we won't use them all.
%  beta_cc will be fixed at 0 and beta_cs will not be independently 
%  estimated, but instead will be set up as the negative of the 
%  forcing coefficient alpha multiplying the set-point function S_0
%  in the contgrol equation.

coefCell    = cell(5,1);
coefCell{1} = coefStrS_S;
coefCell{2} = coefStrS_C;
coefCell{3} = coefStrC_S;
coefCell{4} = coefStrC_C;
coefCell{5} = coefStrC_F;

%  Check the coefficients in coefCell,
%  and compute the total number of parameters defining them.
%  Since each coefficient is a constant and there are five of them,
%  this total will be simply 5.  It would be higher if there more than
%  one parameter defining some coefficient.

%  This is a mandatory step.  Function coefCheck checks the set up and
%  also adds some additional information.

[coefCell, ntheta, nparam] = coefCheck(coefCell);

%  Display the number of parameters to be estimated, here 4 because
%  one of the five coefficients is fixed at 0.

disp(['The number of parameters to estimate is ',num2str(nparam)])

%% 4.  Define XbasisCell containing the basis system for each varible.
%
% We also have to provide a basis system for each variable that is large
% enough to allow for any required sharp curvature in the solution to the
% differential equation system.  
%
% First we set the order of the B-spline basis to 5 so that the first 
% derivative will be smooth when working with a third order derivative in
% the penalty term.  Then we position knots at 41 positions where we willl
% simulate noisy observations.  
% We load these bases into a cell array of length 2.

XbasisCell = cell(2,1);

%  Speed basis, use order 5 splines because speed behaves like
%  the solution of a second order linear equation.

nSorder = 5;
% nSbasis = 7;  % this is for quick demos, but won't be accurate.
% Sbasis  = create_bspline_basis(rng, nSbasis, nSorder);
knotrep = ones(1,nSorder-3);
Sknots   = ...
    [0:2:20,20*knotrep,20:2:40,40*knotrep,40:2:60,60*knotrep,60:2:80];
nSbasis  = length(Sknots) + nSorder - 2;
Sbasis  = create_bspline_basis(rng, nSbasis, nSorder, Sknots);

%  Control basis:  use only order 4 splines because control behaves like
%  the solution of a first order linear equation.

nCorder = 4;
% nCbasis = 6;  % this is for quick demos, but won't be accurate.
% Cbasis  = create_bspline_basis(rng, nCbasis, nCorder);
Cknots   = ...
    [0:2:20,20*knotrep,20:2:40,40*knotrep,40:2:60,60*knotrep,60:2:80];
nCbasis  = length(Cknots) + nCorder - 2;
Cbasis  = create_bspline_basis(rng, nCbasis, nCorder, Cknots);

%  load the two bases into the cell array XbasisCell

XbasisCell{1} = Sbasis;
XbasisCell{2} = Cbasis;

%% 5.  Set up the model structure in cruiseCell
%  The coefficient values used here are for the true system
%  In this version of the model, there are four estimated parameters
%  and one fixed parameter.  The fixed parameter is the coefficent
%  for the controller in the control equation, set to zero to imply
%  a feed back reaction that is proportional speed.

%  Now set up the struct objects for each of the two cells in 
%  cell array cruiseCell.
%  Function make_Variable sets up a struct object whose fields are: 
%  (1) the variable name, 
%  (2) the index of the variable in cruiseCell,
%  (3) the order or derivative for each equation,
%  (4) a cell array object defining the homogeneous terms for the right
%      side of the equation which involve the variables and their
%      derivatives, and
%  (5) a cell array containing struct objects defining the forcing
%      functions for the variable.

% Struct objects for the homogeneous terms in the speed equation
SStr.XCell = cell(2,1);
%  Fields:                 variable ncoef derivative factor
SStr.XCell{1} = make_Xterm(1,       1,    0,         -1); % Speed term
SStr.XCell{2} = make_Xterm(2,       2,    0,          1); % Control term
SStr.FCell    = {};  %  No forcing, this is empty
SStr = make_Variable('Speed', 1, SStr.XCell, SStr.FCell);

% Struct objects for the homogeneous terms in the control equation
CStr.XCell = cell(2,1);
%  Fields:                 variable ncoef derivative factor
CStr.XCell{1} = make_Xterm(1,       3,    0,         -1); % Speed term
CStr.XCell{2} = make_Xterm(2,       4,    0,          0); % Control term

% Struct objects for the forcing term in the control equation
%  Fields:                      ncoef      force factor
CStr.FCell{1} = make_Fterm(         5,    SetPtfd,    1); % Control forcing
CStr = make_Variable('Control', 1, CStr.XCell, CStr.FCell);

%  Cell array for the whole system

cruiseCell = cell(2,1);
cruiseCell{1} = SStr;
cruiseCell{2} = CStr;

%  Check the system specification for consistency
%  This is a mandatory step

cruiseCell = make_Model(XbasisCell, cruiseCell, coefCell);

%% 6. Solving the equations for known parameter values 
%  and generating simulated data from the solutions
%
% In order to simulate data, we need to know the true values of _S(t)_
% and _C(t)_ at the time points where the process is observed.  We also 
% need tospecify the initial state of the system at time 0, which we define 
% to be zero for both variables.  We get the solution by using an initial 
% value approximation algorithm, which is here the Runge-Kutta fourth order
% method coded in Matlab's function |ode45|.
% The function |cruuise_0| evaluates the right side of the equations at a 
% set of time values.  

% We first have a look at the solution at a fine mesh of values by solving
% the equation for the points in |tfine| and plotting them

%  This function evaluates the right side the specialized cruise control 
%  model at time t with variable values in vector Svec and set point
%  values defined by functional data object SetPtfd

odeoptions = odeset;

%  The function cruise_1 defining the right sides of the equation is
%  at the end of this file.

Svec0 = [0,0]';  %  Set initial values of the right sides

%  Use Matlab function ode23 to solve the differential equation
%  given the initial values for a fine mesh of values

[tfine, y0fine] = ode23(@cruise_1, tfine, Svec0, odeoptions, SetPtfd);

%  Use Matlab function ode23 to solve the differential equation
%  given the initial values for times in tobs of values

[tobs, y0obs] = ode23(@cruise_1, tobs, Svec0, odeoptions, SetPtfd);

%  Add some Gaussian random error to the solution values

sigerr = 2;
yobs = zeros(n,2);
yobs(:,1) = y0obs(:,1) + randn(n,1).*sigerr;
yobs(:,2) = y0obs(:,2) + randn(n,1).*sigerr*4;

%  Plot the solution

figure(1)
subplot(2,1,1)
plot(tfine,  y0fine(:,1), 'b-',  [ 0,20],[60,60],'b--', ...
     [20,20],[60,40],     'b--', [20,40],[40,40],'b--', ...
     [40,40],[40,80],     'b--', [40,60],[80,80],'b--', ...
     [60,60],[80,60],     'b--', [60,80],[60,60],'b--', ...
     tobs, yobs(:,1), 'bo', 'LineWidth', 2)
ylabel('\fontsize{16} Speed S (mph)')
subplot(2,1,2)
plot(tfine, y0fine(:,2), 'b-', tobs, yobs(:,2), 'bo', 'LineWidth', 2)
xlabel('\fontsize{16} Time (mins)')
ylabel('\fontsize{16} Control level C')

%%  7. Define cell array yCell, inserting data structs into its two cells.

%  load the data into two struct objects

yStr1.argvals = tobs;
yStr2.argvals = tobs;

% load yobs.txt

yStr1.y       = yobs(:,1);
yStr2.y       = yobs(:,2);

%  load the struct objects into the cell array

yCell = cell(2,1);
yCell{1} = yStr1;
yCell{2} = yStr2;

%%  8. A preliminary evaluation of the function and its derivatives
%  This step is strongly recommended, if only to trap problems in the 
%  set up of the analysis.  

%  Set the initial values of the smoothing parameter $\rho$

rhoVec = 0.9*ones(1,2);  

%  Evaluate at the corrent parameter values, which in this case are
%  the right values since the data are without error

[MSE0, DpMSE0, D2ppMSE0, XfdCell, ...
          df, gcv, ISE, Var_theta, Rmat0, Smat0, fitMap, ...
          DRarray0, DSarray0] = ...
          Data2LD(yCell, XbasisCell, cruiseCell, coefCell, rhoVec);
      
%  Display the total mean squared error and its gradient and hessian

disp(MSE0)
disp(DpMSE0)
disp(D2ppMSE0)

%%  9. Set up a loop through a series of values of rho
%  We know, because the signal is smooth and the data are rough, that the 
%  optimal value of rho will be rather close to one, here we set up a 
%  range of rho values using the logistic transform of equally spaced
%  values between 0 and 5.
%  For each value of rho we save the degrees of freedom, the gcv value,
%  the error sum of squares for each equation, the mean squared errors for 
%  the parameters, and the parameter values.

Gvec    = 0:1:7;
rhoMat  = ones(2,1)*(exp(Gvec)./(1+exp(Gvec)));

%  Matrices for saving values:

nrho    = size(rhoMat,2);  %  smoothing parameter values
dfesave = zeros(nrho,1);   %  degrees of freedom measure
gcvsave = zeros(nrho,1);   %  generalized cross-validation value
MSEsave = zeros(nrho,2);   % total mean squared error
thesave = zeros(nrho,nparam);  %  parameter values

%  Set up values defining optimization process

convrg  = 1e-4;  % convergence criterion
iterlim = 20;    % maximum number of iterations
dbglev  =  1;    % display level (2 gives results for line searches)

%%  10. Loop through rho values

%  Initialize the cell array used for each $\rho$-value set by the
%  the optimization function Data2LD_Opt.  This is reset after each
%  optimization.

coefCell_opti = coefCell;

figure(1)
for irho = 1:nrho   
    rhoVeci = rhoMat(:,irho);
    disp([' ------------------  rhoVeci = ', num2str(rhoVeci'), ...
          ' ------------------'])
    theta_opti = Data2LD_opt(yCell, XbasisCell, cruiseCell, ...
                             coefCell_opti, rhoVeci, ...
                             convrg, iterlim, dbglev);
    coefCell_opti = BAwtvec2cell(theta_opti, coefCell);
    [MSE, DMSE, D2MSE, XfdParCell, df, gcv] = ...
            Data2LD(yCell, XbasisCell, cruiseCell, coefCell_opti, rhoVeci);
    thesave(irho,:) = theta_opti';
    dfesave(irho) = df;
    gcvsave(irho) = gcv;
    x1fd   = getfd(XfdParCell{1});
    x1vec  = eval_fd(tobs, x1fd);
    msex1  = mean((x1vec - y0obs(:,1)).^2);
    x2fd   = getfd(XfdParCell{2});
    x2vec  = eval_fd(tobs, x2fd);
    msex2  = mean((x2vec - y0obs(:,2)).^2);
    MSEsave(irho,1) = msex1;
    MSEsave(irho,2) = msex2;
    Xfd1 = getfd(XfdParCell{1});
    Xfd2 = getfd(XfdParCell{2});    
    Xvec1 = eval_fd(tfine, Xfd1);
    Xvec2 = eval_fd(tfine, Xfd2); 
    disp(['theta = ',num2str(theta_opti')])
    subplot(2,1,1)
    yStr1 = yCell{1};
    plot(tfine,          Xvec1,        'b-',  ...
         tfine,          y0fine(:,1),  'r--', ...
         yStr1.argvals,  yStr1.y,      'bo',  ...
         tfine,          Ufine,        'b:', 'LineWidth', 2)
    axis([0,80,0,100])
    ylabel('\fontsize{16} Speed')
    title(['\fontsize{16} rho = ',num2str(rhoMat(1,irho)), ...
           ',  RMSE = ',num2str(sqrt(msex1))])
    subplot(2,1,2)
    yStr2 = yCell{2};
    plot(tfine,        Xvec2,       'b-',  ...
        tfine,         y0fine(:,2), 'r--', ...
        yStr2.argvals, yStr2.y,     'bo', 'LineWidth', 2)
    axis([0,80,0,400])
    xlabel('\fontsize{16} Time (sec)')
    ylabel('\fontsize{16} Control')
    title(['\fontsize{16} RMSE = ',num2str(sqrt(msex2))])
    pause
end

%%  11. Display some final resultls

%  The parameter values for each set of $\rho$-values

disp(thesave)

%  The values of df, gcv and MSE for each set of $\rho$-values

disp('    rho      df        gcv       MSE:')
disp([rhoMat(1,:)', dfesave(:), gcvsave(:), MSEsave])

%  plot parameters as a function of \rho

figure(2)
ind = 1:nrho;
subplot(1,1,1)
phdl = plot(rhoMat(1,ind)', thesave(ind,:), 'o-');
set(phdl, 'LineWidth', 2)
% axis([0.5,1.0,-2,0.8])
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} \beta(\rho)')
legend('\fontsize{16} \beta_{SS}', '\fontsize{16} \beta_{SC}', ...
       '\fontsize{16} \beta_{CS}', '\fontsize{16} \alpha_{C}', ...
       'location', 'NorthWest')

%  plot root sum of errors for fit as a function of \rho

figure(3)
phdl = plot(rhoMat(1,:)', sqrt(MSEsave), 'o-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} RMSE(\rho)')
legend('\fontsize{16} Speed', '\fontsize{16} Control', ...
       'location', 'West')
   
%  Evaluate the fit at the final valuer of rho:

[MSE, DMSE, D2MSE, XfdParCell, df, gcv, ISE, Var_theta] = ...
     Data2LD(yCell, XbasisCell, cruiseCell, coefCell_opti, rhoMat(:,nrho));

%  data RMSE's 
disp(['Fitting function RMSEs: ',num2str(sqrt(MSE)')])
%  estimated degrees of freedom
disp(['Degrees of freedom = ',num2str(df)])
%  GCV values
disp(['gcv criterion = ',num2str(gcv)])

%  compute estimates of the standard error for each parameter

stddev_opt = sqrt(diag(Var_theta));

%  True parameter values

theta_tru = [1, 0.25, 1, 1]';
theta_opt = thesave(nrho,:);

%  display the estimated and true parameter values 
%  with 95% confidence limits

disp('    True      Est.      Std. Err. Low CI    Upr CI:')
for i=1:nparam
    disp([theta_tru(i), theta_opt(i), stddev_opt(i), ...
          theta_opt(i)-2*stddev_opt(i), ...
          theta_opt(i)+2*stddev_opt(i)])
end

%% ------------------------------------------------------------------------
%  Fitting the constrained model $\beta_22 - \alpha = 0$
%  ------------------------------------------------------------------------

%  This experiment requires that the first experiment be run beforehand.

%  In the first analysis we ignored the fact that that the speed 
%  coefficents $\beta_{CS}$ and $\alpha_C$ are identical since they
%  define a difference.  Thus the four-parameter model over-fits the
%  data by one degree of freedom.  We saw that the two estimated values
%  were close to each other.  Here we re-analyze the data by
%  incorporating the constraint into the analysis so that they are really
%  equal
%
%  Suppose that we have n parameters and these are subject to k
%  constraints.   We can write this in matrix terms as A' p = 0 where
%  A is an n by k matrix, each column of which defines a constraint on
%  the n parameter values in column vector p.  In our example, matrix
%  has a single column and its values are 0, 0, 1, -1, respectively.
%
%  We deal with this in the optimization by defining an n by n-k matrix
%  P in such a way that q = P'*p is a vector of length n-k whose values
%  are unconstrained, and that p = P*q is a vector length n that satisfies
%  the constraints.
%
%  There are several ways of defining matrix P with this property, but
%  perhaps the easiest way is to use the qr decomposition.  If we
%  evaluate statement [Q,R] = qr(A), the resulting matrix Q is n by n and
%  has orthonormal columns, so that Q'*Q = I.  Matrix R, not used here,
%  is upper triangular.
%  Matrix P is defined as the last n-k columns of Q.
%  
%  The optimization function Data2LD_Opt now requires an additional
%  argument:  ParMap is a 4 by 3 matrix that transforms four
%  parameters into a vector space where the last two sum to 0.  We
%  set up ParMap by using the full qr decomposition of the 
%  single constraint vector that expresses the fact that the last two 
%  parameters are equal.

A = [0,0,1,-1]';
[Q,R] = qr(A);
ParMap = Q(:,2:4);

%% Fitting the constrained model

coefCell_opti = coefCell;

%  Loop through rho values

figure(4)
for irho = 1:nrho   
    rhoVeci = rhoMat(:,irho);
    disp([' ------------------  rhoVeci = ', num2str(rhoVeci'), ...
          ' ------------------'])
    theta_opti = Data2LD_opt(yCell, XbasisCell, cruiseCell, ...
                             coefCell_opti, rhoVeci, ...
                             convrg, iterlim, dbglev, ParMap);
    coefCell_opti = BAwtvec2cell(theta_opti, coefCell);
    [MSE, DMSE, D2MSE, XfdParCell, df, gcv] = ...
            Data2LD(yCell, XbasisCell, cruiseCell, coefCell_opti, rhoVeci);
    thesave(irho,:) = theta_opti';
    dfesave(irho) = df;
    gcvsave(irho) = gcv;
    x1fd   = getfd(XfdParCell{1});
    x1vec  = eval_fd(tobs, x1fd);
    msex1  = mean((x1vec - y0obs(:,1)).^2);
    x2fd   = getfd(XfdParCell{2});
    x2vec  = eval_fd(tobs, x2fd);
    msex2  = mean((x2vec - y0obs(:,2)).^2);
    MSEsave(irho,1) = msex1;
    MSEsave(irho,2) = msex2;
    disp(['theta = ',num2str(theta_opti')])
    yStr1 = yCell{1};
    subplot(2,1,1)
    plot(tfine,          Xvec1,        'b-',  ...
         tfine,          y0fine(:,1),  'r--', ...
         yStr1.argvals,  yStr1.y,      'bo',  ...
         tfine,          Ufine,        'b:', 'LineWidth', 2)
    axis([0,80,0,100])
    ylabel('\fontsize{16} Speed')
    title(['\fontsize{16} rho = ',num2str(rhoMat(1,irho)), ...
           ',  RMSE = ',num2str(sqrt(msex1))])
    subplot(2,1,2)
    yStr2 = yCell{2};
    plot(tfine,        Xvec2,       'b-',  ...
        tfine,         y0fine(:,2), 'r--', ...
        yStr2.argvals, yStr2.y,     'bo', 'LineWidth', 2)
    axis([0,80,0,400])
    xlabel('\fontsize{16} Time (sec)')
    ylabel('\fontsize{16} Control')
    title(['\fontsize{16} RMSE = ',num2str(sqrt(msex2))])
    pause
end

%%  The results may now be displayed using the code above.  
%  The last two parameters will now be equal.  The MSE values will
%  inevitably be slightly higher.  In our run, the final unconstrained RMSE
%  was 8.1058 and the constrained version was 8.1074.  This exceedingly
%  small change is consistent with the unconstrained values being 
%  estimated as very close together.

%  The parameter values for each set of $\rho$-values

disp(thesave)

%  The values of df, gcv and MSE for each set of $\rho$-values

disp('    rho      df        gcv       MSE:')
disp([rhoMat(1,:)', dfesave(:), gcvsave(:), MSEsave])

ind = 1:nrho;

%  plot parameters as a function of \rho

figure(5)
subplot(1,1,1)
phdl = plot(rhoMat(1,ind)', thesave(ind,:), 'o-');
set(phdl, 'LineWidth', 2)
% axis([0.5,1.0,-2,0.8])
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} \beta(\rho)')
legend('\fontsize{16} \beta_{SS}', '\fontsize{16} \beta_{SC}', ...
       '\fontsize{16} \beta_{CS}', '\fontsize{16} \alpha_{C}', ...
       'location', 'NorthWest')

%  plot root sum of errors for fit as a function of \rho

figure(6)
phdl = plot(rhoMat(1,:)', sqrt(MSEsave), 'o-');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} \rho')
ylabel('\fontsize{16} RMSE(\rho)')
legend('\fontsize{16} Speed', '\fontsize{16} Control', ...
       'location', 'West')

%% Evaluate aspects of the optimal solution

%  We see that the gcv criterion favors the 8th rho value, 0.9991.

rho_opt   = rhoMat(:,nrho);
theta_opt = thesave(nrho,:)';

%  convert the optimal parameter values to optimal coefCell

coefCell_opt = BAwtvec2cell(theta_opt, coefCell);

%  Evaluate the solution at the optimal solution

[MSE, DpMSE, D2ppMSE, XfdCell, df, gcv, ISE, Var_theta] = ...
      Data2LD(yCell, XbasisCell, cruiseCell, coefCell_opt, rho_opt);
  
%  display results

%  data RMSE's 
disp(['Fitting function RMSEs: ',num2str(sqrt(MSE)')])
%  estimated degrees of freedom
disp(['Degrees of freedom = ',num2str(df)])
%  GCV values
disp(['gcv criterion = ',num2str(gcv)])

%  compute estimates of the standard error for each parameter

stddev_opt = sqrt(diag(Var_theta));

%  True parameter values

theta_tru = [1, 0.25, 1, 1]';
theta_opt = thesave(nrho,:);

%  display the estimated and true parameter values 
%  with 95% confidence limits

disp('    True      Est.      Std. Err. Low CI    Upr CI:')
for i=1:nparam
    disp([theta_tru(i), theta_opt(i), stddev_opt(i), ...
          theta_opt(i)-2*stddev_opt(i), ...
          theta_opt(i)+2*stddev_opt(i)])
end

%%  Things that you might like to try:

% 1.  Play around with sigerr in step 6. How large does it have to be
%     to have the estimation break down?  
% 2.  If it breaks down for 41 observations, what if there are 81?
% 3.  What if you only had data for speed?  Set yCell{2} to {} in
%     step 7.  How small would sigerr have to be to get good recovery
%     of the controller level?
% 4.  Do the same thing if there are only the controller observations.
% 5.  What if beta_{22} is a positive value instead of 0?
% 6.  What some of the coefficients vary over time?  As, for example,
%     driving conditions change.  Try defining one or more of them as
%     functional data objects defined by, for example, order 2 splines
%     with no interior knots. 

%% A single order two Feedback equation
%
% We can, without losing much modelling power, simplify the equations by 
% assuming that the gain $K = \alpha/\beta$ is fixed at 1.0 and the 
% controller reaction speed is so large that the controller responds 
% virtually immediately to a change in speed/set--point discrepancy.  
% That is, $DC = S_0 - S$, or 
%%
% $$C(t) = \int_0^t (S_0(u) - S(u)) \, du + C_0$$
%
% where $C_0$ is an arbitrary constant of integration. Inserting the right 
% side of this equation into the _DS_ equation and differentiating both 
% sides of the result, we have the single second order equation
%%
% $$D^2 S(t) = -\beta_0 S(t) - \beta_1 DS(t) + \alpha S_0(t) + C_0$$
%
% where we have the constraints that $-\beta_0 = \alpha$ and 
% $-\beta_1 > 0$.  Now we assume that only speed is observed.
%
%% Generate some replicated sample data for analysis
% Now we will improve the power of the estimation by assuming that the
% experiment has been repeated ten times.  The forcing functions will 
% vary from replication to replication, but their coefficient will not.
% That, SetPtfd now contains 10 functions instead of one, varying randomly.

%%  1. Defining the problem
% Set up the time span, a set of observation times, and a fine time grid 
% for plotting:

T     = 80;     %  maximum time in seconds
rng   = [0,T];  %  range of times to be used
n     = 41;     %  number of observations
tobs  = linspace(0,T,n)';  %  times of the observations
nfine = 501;  % number of points in a fine mesh for plotting purposes
tfine = linspace(0,T,nfine)';  %  mesh values

%% 2.  Set up the set-point forcing function $S_0(t)$.
%
% The set-point function uses an order 1 step function B-spline basis.  The 
% knots are placed at the points where the set-point changes values.
%
% Here we define 10 randomly varying set functions

nrep = 10;  %  the number of replications
% nrep = 1;  %  the number of replications

stepbasis = create_bspline_basis(rng, 4, 1, [0,20,40,60,80]);
stepcoef = [60;40;80;60];

sigstep = 0.1;
stepcoef   = stepcoef*exp(randn(1,nrep)*sigstep);
SetPtfd    = fd(stepcoef, stepbasis);

%% 3.  Set up the coefficient dictionary in cell array coefCell
%
%  The total number of coefficients defining the estimated coefficient 
%  functions is three because each coefficient function has only one 
%  coefficient, and only three of them are estimated.

%  Set up a constant basis over [0,T], 
%  all of our coefficients will constants

conbasis  = create_constant_basis(rng);
confd     = fd(1,conbasis);
confdnrep = fd(ones(1,nrep),conbasis);
confdPar  = fdPar(conbasis);

coefStr1 = make_Coef(conbasis, 1/4, true);
coefStr2 = make_Coef(conbasis,   1, true);
coefStr3 = make_Coef(conbasis,   1, true);
coefStr4 = make_Coef(conbasis,   1, true);

coefCell = cell(3,1);
coefCell{1} = coefStr1;
coefCell{2} = coefStr2;
coefCell{3} = coefStr3;

coefCell = cell(4,1);
coefCell{1} = coefStr1;
coefCell{2} = coefStr2;
coefCell{3} = coefStr3;
coefCell{4} = coefStr4;

[coefCell, ntheta, nparam] = coefCheck(coefCell);

disp(['ntheta = ',num2str(ntheta)])
disp(['nparam = ',num2str(nparam)])

%% 4.  Define XbasisCell containing the basis system for each varible.
%
% We also have to provide a basis system for each variable that is large
% enough to allow for any required sharp curvature in the solution to the
% differential equation system.  
%
% First we set the order of the B-spline basis to 5 so that the first 
% derivative will be smooth when working with a third order derivative in
% the penalty term.  Then we position knots at 41 positions where we willl
% simulate noisy observations.  
% We load these bases into a cell array of length 2.

XbasisCell = cell(1);

nSorder = 5;
% nSbasis = 7;  % this is for quick demos, but won't be accurate.
% Sbasis  = create_bspline_basis(rng, nSbasis, nSorder);
Sknots = tobs;
nSbasis  = length(Sknots) + nSorder - 2;
Sbasis  = create_bspline_basis(rng, nSbasis, nSorder, Sknots);

%  load the basis into the cell array XbasisCell

XbasisCell{1} = Sbasis;

%% 5.  Set up the model structure in cruiseCell

%  First term in S equation
Xterm1 = make_Xterm(1, 1, 0, -1);
%  Second in term in S equation
Xterm2 = make_Xterm(1, 2, 1, -1);

% Struct object for set point forcing term in C equation
Fterm1 = make_Fterm(3, SetPtfd,   1);
% Struct object for constant forcing term C_0
Fterm2 = make_Fterm(4, confdnrep, 1);

% Struct object for the speed equation

SStr.XCell = cell(2,1);
SStr.XCell{1} = Xterm1;
SStr.XCell{2} = Xterm2;

SStr.FCell = cell(1);
SStr.FCell{1} = Fterm1;

SStr.FCell = cell(2,1);
SStr.FCell{1} = Fterm1;
SStr.FCell{2} = Fterm2;

SpeedStr = make_Variable('Speed', 2, SStr.XCell, SStr.FCell);

%  Cell array for the whole system

cruiseCell = cell(1);
cruiseCell{1} = SpeedStr;

cruiseCell = make_Model(XbasisCell, cruiseCell, coefCell);

%% 6. Solving the equations for known parameter values 
%  and generating simulated data from the solutions

y0fine = zeros(nfine,nrep);
y0obs  = zeros(41, nrep);

Svec0 = [0,0]';

odeoptions = odeset;

for irep=1:nrep
    [tfine, y0finei] = ...
        ode23(@cruise_2, tfine, Svec0, odeoptions, ...
              cruiseCell, coefCell, SetPtfd(irep));
    [tobs, y0obsi]   = ...
        ode23(@cruise_2, tobs,  Svec0, odeoptions, ...
              cruiseCell, coefCell, SetPtfd(irep));
    y0fine(:,irep)   = y0finei(:,1);
    y0obs(:,irep)    = y0obsi(:,1);
end

sigerr = 2;

yobs = zeros(41,nrep);
for irep=1:nrep
    yobs(:,irep) = y0obs(:,irep) + randn(41,1)*sigerr;
end

yStruct.argvals = tobs;
yStruct.y       = yobs;

%  Define yCell

yCell = cell(1);
yCell{1} = yStruct;

%% 7.  Compute starting values
%
% We didn't do this for the two variable case because we used the
% values that defined the differential equation.  This is cheating 
% since in practice we wouldn't have these.
%
% Here we now generate initial parameter values directly from the data.
% We use regression analysis for this.  First, smooth the data.  The
% evaluate the smooth and its first and second derivatives at the 
% observation points, as well as evaluating the forcing function.  
% Assemble these elements into a vectorized version of the second 
% derivative and a covariate matrix with three columns, the first 
% containing the vectorized smooth values, the second the vectorized first
% derivative values, and the final the vectorized forcing function values.

%  smooth the data

lambda = 1e1;
XfdPar = fdPar(Sbasis,2,lambda);

% evaluate the smooth, its first and second derivatives and the forcing
% functions

yfd    = smooth_basis(tobs, yobs, XfdPar);
ymat   = eval_fd(tobs, yfd);
Dymat  = eval_fd(tobs, yfd, 1);
D2ymat = eval_fd(tobs, yfd, 2);
Umat   = eval_fd(tobs, SetPtfd);

% assemble this into the linear model

D2Yvec = zeros(nrep*n,1);
Zmat   = zeros(nrep*n,4);
m2 = 0;
for irep = 1:nrep
    m1 = m2 + 1;
    m2 = m2 + n;
    D2Yvec(m1:m2) = D2ymat(:,irep);
    Zmat(m1:m2,1) =  -ymat(:,irep);
    Zmat(m1:m2,2) =  Dymat(:,irep);    
    Zmat(m1:m2,3) =   Umat(:,irep);
    Zmat(m1:m2,4) =   ones(n,1);
end

% evaluate the regression coefficients

beta = Zmat\D2Yvec;

% set up an initial estimates of parameters

coefStr1.parvec = mean([beta(1),beta(3)]);
coefStr2.parvec = beta(2);
coefStr3.parvec = mean([beta(1),beta(3)]);
coefStr4.parvec = beta(4);

coefCell{1} = coefStr1;
coefCell{2} = coefStr2;
coefCell{3} = coefStr3;
coefCell{4} = coefStr4;
[coefCell,ntheta] = coefCheck(coefCell);
disp(['ntheta = ',num2str(ntheta)])

%%  8. A preliminary evaluation of the function and its derivatives
%  This step is strongly recommended, if only to trap problems in the 
%  set up of the analysis.  

rho = 0.5;

clear functions

[MSE, DpMSE, D2ppMSE] = Data2LD(yCell, XbasisCell, cruiseCell, ...
                                coefCell, rho);
disp(MSE)
disp(DpMSE)
disp(D2ppMSE)

%%  9. Set up a loop through a series of values of rho
% Now the first and the third parameters must be equal and of opposite 
% sign.
%

conv    = 1e-6;
iterlim = 50;
dbglev  = 1;

%  set up constraint that coefficients 1 and 3 must add to 0

A = [-1,0,1]';
[Q,R] = qr(A);
ParMap = Q(:,2:3);

%  set up a sequence of analyses with varying value of rho

Gvec = 0:7;
rhovec = exp(Gvec)./(1+exp(Gvec));
nrho   = length(rhovec);

%%  10. Loop through rho values

%  Initialize the cell array used for each $\rho$-value set by the
%  the optimization function Data2LD_Opt.  This is reset after each
%  optimization.

coefCell_opt = coefCell;

thetasave = zeros(nrho,4);
dfsave    = zeros(nrho,1);
gcvsave   = zeros(nrho,1);
figure(2)
for irho = 1:nrho
    rhoi = rhovec(irho);
    disp(['rho = ',num2str(rhovec(irho))])
    theta_opt = Data2LD_opt(yCell, XbasisCell, ...
        cruiseCell, coefCell, rhoi, conv, iterlim, dbglev);
    coefCell_opt = BAwtvec2cell(theta_opt, coefCell_opt);
    [SSE, DSSE, D2SSE, XfdParCell, df, gcv, ISE] = ...
        Data2LD(yCell, XbasisCell, cruiseCell, coefCell, rhoi);
    thetasave(irho,:) = theta_opt';
    dfsave(irho)      = df;
    gcvsave(irho)     = gcv;
    disp(theta_opt')
    disp(df)
    disp(gcv)
    subplot(1,1,1)
    yStruct1 = yCell{1};
    Xfd  = getfd(XfdParCell{1});
    Xmat = eval_fd(tfine, Xfd);
    plot(tfine, Xmat, 'b-', yStruct1.argvals, yStruct1.y, 'bo', ...
         'LineWidth', 2);
    xlabel('\fontsize{16} Time (sec)')
    ylabel('\fontsize{16} Speed (km/h)')
    title(['\fontsize{16} rho = ',num2str(rhoi)])
    pause
end

disp('rho-values, degrees of freedom and gcv values:')
ind = 1:nrho;
disp([rhovec(ind)',dfsave(ind),gcvsave(ind)])

%% Evaluation the solution for best gcv value

irho = 2;

disp(['theta: ',num2str(thetasave(irho,:))])

[SSE, DSSE, D2SSE, XfdParCell, df, gcv] = ...
        Data2LD(yCell, XbasisCell, cruiseCell, coefCell, rhovec(irho));

%  display the estimated solutions

Xfd  = getfd(XfdParCell{1});
Xmat = eval_fd(tfine, Xfd);

figure(2)
subplot(1,1,1)
yStruct1 = yCell{1};
phdl = plot(tfine, Xmat, 'b-', yStruct1.argvals, yStruct1.y, 'bo');
set(phdl, 'LineWidth', 2)
xlabel('\fontsize{16} Time (sec)')
ylabel('\fontsize{16} Speed (km/h)')


