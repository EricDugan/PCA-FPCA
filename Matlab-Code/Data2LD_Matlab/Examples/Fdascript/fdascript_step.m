%%  -----------------------------------------------------------------------
%                 Registered Handwriting Data
%  -----------------------------------------------------------------------

% These data are the X-Y coordinates of 20 replications of writing
% the script "fda".  The subject was Jim Ramsay.  Each replication
% is represented by 1401 coordinate values.  The scripts have been 
% extensively pre-processed.  They have been adjusted to a common
% length that corresponds to 2.3 seconds or 2300 milliseconds, and
% they have already been registered so that important features in
% each script are aligned.
% 
% This analysis is designed to illustrate techniques for working
% with functional data having rather high frequency variation and
% represented by thousands of data points per record.  Comments
% along the way explain the choices of analysis that were made.
% 
% The final result of the analysis is a third order linear 
% differential equation for each coordinate forced by a 
% constant and by time.  The equations are able to reconstruct
% the scripts to a fairly high level of accuracy, and are also
% able to accommodate a substantial amount of the variation in
% the observed scripts across replications.  by contrast, a 
% second order equation was found to be completely inadequate.
% 
% An interesting suprise in the results is the role placed by
% a 120 millisecond cycle such that sharp features such as cusps
% correspond closely to this period.  This 110-120 msec cycle
% seems is usually seen in human movement data involving rapid
% movements, such as speech, juggling and so on.

%  Last modified 20 November 2017

%% Add paths to functional data functions and to the data

%  Input the data.  These 20 records have already been
%  normalized to a common time interval of 2300 milliseconds
%  and have beeen also registered so that prominent features
%  occur at the same times across replications.
%  Time will be measured in milliseconds and space in meters.
%  The data will require a small amount of smoothing, since
%  an error of 0.5 mm is characteristic of the OPTOTRAK 3D
%  measurement system used to collect the data.

%  paths to necessary folders

addpath('../fdaM')
addpath('Examples/Fdascript')

%  load the data

load fda

%  set up data array, time values and range 

fdaarray = fda.fdaarray;
fdatime  = fda.fdatime;  % size(fdatime) = 1401   1

n        = 1401; % number of observations (coordinate values)
sec      = linspace(0,2.3,n)';
centisec = sec*100;  %  at better choice for using splines
fdarange = [0,230];

%  set which of the 20 records are to be analyzed and displayed

nrecord = 20;
recordindex = 1:20;
recordindex = 1;

%  create an array object containing observations in centimetres

XY = reshape(fdaarray, [n,nrecord,2]).*1000;
X = XY(:,:,1);
Y = XY(:,:,2);

%  Here we position a number at each of 20 equally spaced times
%  located within 20 equal-sized intervals

nevent    = 20;
eventtime = linspace(fdarange(1)/10,fdarange(2)/10,2*nevent+1);
eventind  = 2:2:2*nevent;

%%  The basis functions will be B-splines

%  32 basis functions per second, 4 per 8 herz cycle
%  2.3 seconds 2.3*32=73.6, 
%  This implies norder + no. of interior knots = 74 - 1 + 6 = 79 
%  basis functions.  
nbasis   =   79; 
norder   =    6; %  order 6 for a smooth 2nd deriv.
fdabasis = create_bspline_basis(fdarange, nbasis, norder);

%  Don't need penalizing function
%  fitting error is smaller than 0.05 mm, the known error level in
%  the OPTOTRACK equipment.

% fdafd  = fd(zeros(nbasis,20,2), fdabasis);
lambda = 0;
fdaPar = fdPar(fdabasis, 5, lambda);

%  Add suitable names for the dimensions of the data.

% fdafd_fdnames{1} = 'Milliseconds';
% fdafd_fdnames{2} = 'Replications';
% fdafd_fdnames{3} = 'Metres';
% fdafd = putnames(fdafd, fdafd_fdnames);

%  set up the functional data structure 

[fdafd, df, gcv] = smooth_basis(centisec, XY(:,recordindex,:), fdaPar);

%  display degrees of freedom and total GCV criterion

disp(['degrees of freedom = ',num2str(df)])  %  79
totalgcv = sum(sum(gcv)); % 1.2257e-09.  0.12257
disp(['total GCV = ',num2str(totalgcv)])  
RMSgcv = sqrt(totalgcv)*1000; % 0.035009  350.0934
disp(['RMS GCV = ',num2str(RMSgcv)])  

%%  plot the fit to the data.  It's good to within plotting
%  accuracy, but the title displays RMSE values

figure(1)
plotfit_fd(XY(:,recordindex,:), centisec, fdafd);

%%  -----------------------------------------------------------------------
%             Analyze the two equations seperately
%   -----------------------------------------------------------------------


%% set up abasis for X

XbasisCell = cell(1,1);
XbasisCell{1} = fdabasis;

%% Set up the cell array coefCell

%Constant basis

rng = [0,230];

basisobjC = create_constant_basis(rng);
confd     = fd(ones(getnbasis(basisobjC),1),basisobjC);
confdPar  = fdPar(confd);
Ufd       = fd(1,basisobjC);

%  values controlling optimization

dbglev   = 1;    %  debugging level
iterlim  = 50;    %  maximum number of iterations
conv     = [1e-8, 1e-4];  %  convergence criterion

%  set up sequence of rho values
rhoVec1 = 0.5;

gamvec = 0:7;
rhoVec = exp(gamvec)./(1+exp(gamvec));
nrho   = length(rhoVec);

%  Forcing function basis

nAorder   = 1;
% nratio    = 5;
% ratiovec = linspace(0.01,0.1,nratio);
% ratio  = ratiovec(iratio);
nAbasis   = nevent;
nAknots   = nAbasis + 1;
Aknots    = linspace(rng(1),rng(2),nAknots);
Abasisobj = create_bspline_basis(rng,nAbasis,nAorder,Aknots);

%% Beta coefficients

%  arrays to hold results over rho

%% Order 2 model with step input

% Set up coefficient struct objects for X term and forcing function

% Fields:           fun         parvec            estimate 
coefStr1 = make_Coef(confdPar,  0.04,             1);
coefStr2 = make_Coef(Abasisobj, zeros(nAbasis,1), 1);

% cell array containing the coefficient structs

coefCell    = cell(2,1);
coefCell{1} = coefStr1;
coefCell{2} = coefStr2;

% Check the coefficients

[coefCell, ntheta] = coefCheck(coefCell);
disp(['ntheta = ',num2str(ntheta)]) 

%%  Define each of the two terms in the right sides of equations

%  Set up the single term in the equation D^2X(t) = -beta X(t)

% Fields:          variable derivative ncoef factor
XStr  = make_Xterm(1,       1,         1,    -1);
XCell = {XStr};

% Struct object for constant forcing term in X equation

%  Set up forcing term in the equation D^2X(t) = -beta X(t) + alpha F(t)

% Fields:          ncoef Ufd  factor
FStr  = make_Fterm(2,    Ufd, 1);
FCell = {FStr};

%% Set up the cell array modelCell

% Struct object for the single variable in the equation

% Fields:                name order XCell  FCell
coordStr = make_Variable('X', 2,    XCell, FCell);

%  Cell array for the whole system

modelCell = {coordStr};

%  check the system specification for consistency

modelCell = make_Model(XbasisCell, modelCell, coefCell);

%% -----------------------------
%          X direction
% ------------------------------

%%  Single evaluation in order to set up the 4-way tensors

yStrX.argvals = centisec;
yStrX.y       = XY(:,recordindex,1);
yCellX = {yStrX};

rhoVec = 0.5;

clear functions;

tic;
[fvec, grad] = Data2LD(yCellX, XbasisCell, modelCell, coefCell, rhoVec);
toc

%  set up sequence of rho values

gamvec = 0:7;
rhoVec = exp(gamvec)./(1+exp(gamvec));
nrho   = length(rhoVec);

dbglev  = 0;
iterlim = 20;
convrg  = 1e-6;

%%  Set up arrays to hold results over rho

nrec = length(recordindex);

%  Dimensions:      parmeters  rho   record   coordinate
thetastore = zeros(ntheta,    nrho, nrec, 2);
dfstore    = zeros(           nrho, nrec, 2);
gcvstore   = zeros(           nrho, nrec, 2);
MSEstore   = zeros(           nrho, nrec, 2);
ISEstore   = zeros(           nrho, nrec, 2);

%%  loop through the 20 replications analyzing the X coordinate

yStrX.argvals = centisec;
yStrX.y       = XY(:,recordindex,1);
yCellX = {yStrX};

tic;
for record = recordindex
    
    disp(['------------------------- X Record ', num2str(record), ...
          '  -------------------------'])
    
    yStrX.y = XY(:,record,1);
    
    yCell = cell(1,1);
    yCell{1} = yStrX;
    
    %%  step through rho values, optimizing at each step
    
    coefCell_opti = coefCheck(coefCell);
    
    ind = 1:n;
    
    for irho = 1:nrho
        rhoi = rhoVec(irho);
        disp(['rho = ',num2str(rhoi')])
        theta_opti = Data2LD_opt(yCell, XbasisCell, modelCell, ...
            coefCell_opti, rhoi, convrg, iterlim, dbglev);
        coefCell_opti = BAwtvec2cell(theta_opti, coefCell_opti);
        [MSE, DSSE, D2SSE, XfdParCell, df, gcv, ISE] = ...
            Data2LD(yCell, XbasisCell, modelCell, coefCell_opti, rhoi);
        thetastore(:,irho,record,1) = theta_opti;
        MSEstore(    irho,record,1) = MSE;
        ISEstore(    irho,record,1) = ISE;
        dfstore (    irho,record,1) = df;
        gcvstore(    irho,record,1) = gcv;
    end
    
end
toc

%%  loop through the 20 replications analyzing the Y coordinate

yStrY.argvals = centisec;
yStrY.y       = XY(:,recordindex,2);
yCellY = {yStrY};

tic;
for record = recordindex
    
    disp(['------------------------- Y Record ', num2str(record), ...
        '  -------------------------'])
    
    %%  step through rho values, optimizing at each step
    
    coefCell_opti = coefCheck(coefCell);
    
    ind = 1:n;
    
    for irho = 1:nrho
        rhoi = rhoVec(irho);
        disp(['rho = ',num2str(rhoi')])
        theta_opti = Data2LD_opt(yCellY, XbasisCell, modelCell, ...
            coefCell_opti, rhoi, conv, iterlim, dbglev);
        coefCell_opti = BAwtvec2cell(theta_opti, coefCell_opti);
        [SSE, DSSE, D2SSE, XfdParCell, df, gcv, ISE, Var_theta] = ...
            Data2LD(yCellY, XbasisCell, modelCell, coefCell_opti, rhoi);
        thetastore(:,irho,record,2) = theta_opti;
        MSEstore(    irho,record,2) = MSE;
        ISEstore(    irho,record,2) = ISE;
        dfstore(     irho,record,2) = df;
        gcvstore(    irho,record,2) = gcv;
    end
    
end
toc

%% Plots for the highest rho value

nbeta = 1;
for record = recordindex
    
    disp(['-------------------------  Record ', num2str(record), ...
        '  -------------------------'])
    
    irho = nrho;
    rhoi = rhoVec(irho);
    
    Xtheta_opti = thetastore(:,irho,record,1);
    XcoefCell_opti = BAwtvec2cell(Xtheta_opti, coefCell);
    [SSE, DSSE, D2SSE, XfdParCell_X, df, gcv, ISE, Var_thetaX] = ...
        Data2LD(yCellX, XbasisCell, modelCell, XcoefCell_opti, rhoi);
    
    Ytheta_opti = thetastore(:,irho,record,2);
    YcoefCell_opti = BAwtvec2cell(Ytheta_opti, coefCell);
    [SSE, DSSE, D2SSE, XfdParCell_Y, df, gcv, ISE, Var_thetaY] = ...
        Data2LD(yCellY, XbasisCell, modelCell, YcoefCell_opti, rhoi);
    
    %  -------  display values of beta  --------
    
    disp(['betaX = ',num2str(Xtheta_opti(1))])
    disp(['betaY = ',num2str(Ytheta_opti(1))])
      
    %  -------  plot fits to coordinate data  --------
    
    Xfd   = getfd(XfdParCell_X{1});
    Xfit  = eval_fd(centisec,Xfd);
    XRMSE = sqrt(mean((yStrX.y - Xfit).^2));
    
    Yfd   = getfd(XfdParCell_Y{1});
    Yfit  = eval_fd(centisec,Yfd);
    YRMSE = sqrt(mean((yStrY.y - Yfit).^2));
    
    %  -----------------  plot fit to script  ----------------------
    
    Xvec = eval_fd(centisec, Xfd);
    Yvec = eval_fd(centisec, Yfd);
    ind50 = 1:5:n;  % plot every 5th point
    
    figure(1)
    subplot(1,1,1)
    phdl = plot(Xvec, Yvec, 'b-', X(ind50,record), Y(ind50,record), 'bo');
    set(phdl, 'LineWidth', 1)
    grid on
    %axis([-20,20,-12,10,-1,4])
    xlabel('\fontsize{16} X(t) (cm)')
    ylabel('\fontsize{16} Y(t) (cm)')
    title(['\fontsize{16} record = ',num2str(record)])
    
    %  -----------------  plot fit to coordinates  ----------------------

    figure(2)
    subplot(2,1,1)
    phdl = plot(centisec, yStrX.y, 'b--', ...
                centisec, Xfit,    'b-');
    set(phdl, 'LineWidth', 2)
    ylabel('\fontsize{16} X(t)')
    title(['\fontsize{16} XRMSE = ',num2str(XRMSE)])
    
    subplot(2,1,2)
    phdl = plot(centisec, yStrY.y, 'b--', ...
        centisec, Yfit, 'b-');
    set(phdl, 'LineWidth', 2)
    xlabel('\fontsize{16} t (msec)')
    ylabel('\fontsize{16} Y(t)')
    title(['\fontsize{16} YRMSE = ',num2str(YRMSE)])
    
    %  -------  plot forcing functions  -------
    
    AfdX = fd(Xtheta_opti(nbeta+1 : nAbasis+nbeta), Abasisobj);
    AcoefX = getcoef(AfdX);
    XAvec = eval_fd(centisec, AfdX);
    
    AfdY = fd(Ytheta_opti(nbeta+1 : nAbasis+nbeta), Abasisobj);
    AcoefY = getcoef(AfdY);
    YAvec = eval_fd(centisec, AfdY);
    
    figure(3)
    subplot(2,1,1)
    fdascript_force_plot(Aknots, AcoefX , ['\alpha_','X(t)'], [-1.1,1.1])    
    subplot(2,1,2)
    fdascript_force_plot(Aknots, AcoefY , ['\alpha_','Y(t)'], [-1.1,1.1])    
    
    %  -----------  plot fit to second derivative  -----------------
    
    ind50 = 1:50:n;
    ind   = 1:n;
    
    D2Xvec = eval_fd(centisec, Xfd, 2);
    D2Xfit = Xtheta_opti(1).*Xvec + XAvec;
    D2Xfit(ind50) = D2Xvec(ind50);
    DXRMSE = sqrt(mean((D2Xvec - D2Xfit).^2));
    
    D2Yvec = eval_fd(centisec, Yfd, 2);
    D2Yfit = Ytheta_opti(1).*Yvec + YAvec;
    D2Yfit(ind50) = D2Yvec(ind50);
    DYRMSE = sqrt(mean((D2Yvec - D2Yfit).^2));
    
    figure(4)
    subplot(2,1,1)
    phdl = plot(centisec(ind), D2Xvec(ind), 'b--', ...
        centisec(ind), D2Xfit(ind), 'b-', ...
        rng, [0,0], 'b--');
    set(phdl, 'LineWidth', 2)
    %axis([rng,-0.2,0.2])
    ylabel('\fontsize{16} D^2X(t)')
    title(['\fontsize{16} L-XRMSE = ',num2str(DXRMSE)])
    
    subplot(2,1,2)
    phdl = plot(centisec(ind), D2Yvec(ind), 'b--', ...
        centisec(ind), D2Yfit(ind), 'b-', ...
        rng, [0,0], 'b--');
    set(phdl, 'LineWidth', 2)
    %axis([rng,-0.2,0.2])
    xlabel('\fontsize{16} t (msec)')
    ylabel('\fontsize{16} D^2Y(t)')
    title(['\fontsize{16} L-YRMSE = ',num2str(DYRMSE)])
    
    if length(recordindex) > 1
        pause
    end
    
end

