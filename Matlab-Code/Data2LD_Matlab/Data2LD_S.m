function [Smat, DSarray] = Data2LD_S(XbasisCell, modelCell, coefCell, ...
                                     rhoVec, nthetaH, nthetaF, nrep, nforce)
%  Data2LD ... stands for "Data to Linear Dynamics"
%  Data2LD_S computes the penalty matrix S associated with the forcing
%  portion of a linear differential operator L, as well as its partial
%  derivatives with  respect to the parameter vector.
%  For a single variable whose approximated in terms of an exansion in
%  terms of a vector \phi of basis functions, S is
%                 S = \int [L \phi(t)] U' dt.
%  S has dimensions K and NREP, where K is the number of basis
%  functions in the expansion of the variable, NFORCE is the number of
%  forcing functions, and NREP is the number of replications.  The
%  forcing functions are assumed to vary from one replication to another.
%
%  If multiple variables are involved, then S is a composite matrix
%  constructed from inner products and cross-products of the basis 
%  function vectors associate with each variable.  It's dimension will be
%  \sum K_i by NFORCE*NREP.
%
%  Arguments:
%
%  XBASISCELL ... A functional data object or a BASIS object.  If so, the 
%               smoothing parameter LAMBDA is set to 0.
%  MODELCELL  ... A cell aray of length NVAR. Each cell contains a 
%                struct object with members:              
%                XCell ... cell array of length number of homogeneous terms
%                          Each cell contains a struct object with members:
%                          WfdPar     ... fdPar object for the coefficient
%                          variable   ... the index of the variable
%                          derivative ... the order of its derivative
%                          ncoef      ... if coefficient estimated, its 
%                                         location in the composite vector 
%                          factor     ... a scalar multiplier (def. 1)
%                FCell ... cell arrau of length number of forcing terms
%                          Each cell contains a struct object with members:
%                          AfdPar ... an fdPar object for the coefficient
%                          Ufd    ... an fd object for the forcing function
%                          ncoef  ... if coefficient estimated, its 
%                                     location in the composite vector 
%                          factor ... a scalar multiplier (def. 1)
%                order     ... the highest order of derivative
%                name      ... a  tag for the variable
%                nallXterm ... the number of homogeneous terms
%                nallFterm ... the number of forcing functions
%  COEFCELL  ... A cell array of length NCOEF containing struct objects
%                with fields:
%               parvec   ... a vector of parameters
%               estimate ... 0, held fixed, otherwise, estimated 
%               coeftype ... homogeneous or forcing
%               fun      ... functional basis, fd, or fdPar object, 
%                            or a struct object for a general function 
%                            with fields:
%                 fd      ... function handle for evaluating function
%                 Dfd     ... function handle for evaluating 
%                             partial derivative with respect to parameter
%                 more    ... object providing additional information for 
%                             evaluating coefficient function
%  RHOVEC    ... A vector of length NVAR containing values in [0,1].  
%                The data sums of squares are weighted by P and 
%                the roughness penalty by 1-P.
%  NREP       ... The number of replications of the system.  

%  Last modified 31 March 2020

%  ------------------------------------------------------------------------
%                         Set up analysis
%  ------------------------------------------------------------------------

%  compute number of variables

nvar = length(modelCell);

%  Set up a vector NCOEFVEC containing number of coefficients used
%  for the expansion of each variable

ncoefvec = zeros(nvar,1);
for ivar=1:nvar
    ncoefvec(ivar) = getnbasis(XbasisCell{ivar});
end
ncoefsum = sum(ncoefvec);
ncoefcum = [0;cumsum(ncoefvec)];

%  get the width of the time domain

Xrange = getbasisrange(XbasisCell{1});
T      = Xrange(2) - Xrange(1);

coefStri.parvec   = 1;
coefStri.estimate = 0;
coefStri.fun      = fd(1,create_constant_basis(Xrange));

%--------------------------------------------------------------------------
%       Compute the penalty vector or matrix Smat(theta)
%--------------------------------------------------------------------------

%  If there are no forcing functions, return Smat and DSarray as empty

if sum(nforce) == 0
    Smat = [];
    DSarray = [];
    return
end

%  Set up Smat to have as many rows as basis functions for representing
%  all variables and as many columns as replications of forcing functions.

Smat = zeros(ncoefsum,nrep);

%  Loop through variables

m2 = 0;
for ivar=1:nvar
    m1  = m2 + 1;
    m2  = m2 + ncoefvec(ivar);
    ind = m1:m2;
    if nforce(ivar) > 0
        % This variable is forced, set up order, weight, basis and 
        % both X and F terms
        modelStructi = modelCell{ivar};
        order    = modelStructi.order;
        weighti  = modelStructi.weight;
        Xbasisi  = XbasisCell{ivar};
        nXbasisi = ncoefvec(ivar);
        nXtermi  = modelStructi.nallXterm;
        nFtermi  = modelStructi.nallFterm;
        for jforce = 1:nFtermi
            % For this forcing term select parameter vectors for 
            % coefficient, factor, known forcing function and whether
            % coefficient has a basis expanstion (funtypej 0) or is
            % user-defined (funtypej 1)
            modelStructij = modelStructi.FCell{jforce};
            %  get details for forcing term
            [coefStrj, factorj, Ufdj, Avecj, Ucoefj, nUbasisj, ~, ...
             nAbasisj, funtypej] = ...
                   getForceTerm(modelStructij, coefCell);           
            %  Crossproducts of homogeneous terms with forcing terms
            for iw=1:nXtermi
                % For this homogeneous term obtain variable basis,
                % coefficient parameter vector and type of coefficient
                modelStructiw = modelStructi.XCell{iw};
                [coefStrw, ivw, factorw, Bvecw, nWbasisw, ~, ...
                 derivw, funtypew] = ...
                       getHomoTerm(modelStructiw, coefCell);
                indw     = (ncoefcum(ivw)+1):ncoefcum(ivw+1);
                Xbasisw  = XbasisCell{ivw};
                nXbasisw = ncoefvec(ivw);
                if funtypej || funtypew
                    % Either the homogeneous coefficient or the 
                    % forcing coefficient or both are user-defined,
                    % use numerical integration to get inner products
                    Smatjw = ...
                        inprod_basis_Data2LD(Xbasisw,  Ufdj,     ...
                                             coefStrw, coefStrj, ...
                                             derivw,   0);
                else
                    % Both coefficients are basis function expansions,
                    % use previously computed inner product values
                    BAtenswj = modelStructi.BAtens{iw,jforce};
%                     disp([nXbasisw, nWbasisw, nUbasisj, nAbasisj, nrep])
%                     disp([size(Bvecw),size(Bvecw),size(Avecj), ...
%                           size(Ucoefj), size(BAtenswj)])
                    Smatjw = inprodijw(nXbasisw, nWbasisw, ...
                                       nUbasisj, nAbasisj, nrep, ...
                                       Bvecw, Avecj, Ucoefj, BAtenswj);
                end
                Smatjw = weighti.*factorj*factorw*rhoVec(ivar).*Smatjw/T;
                Smat(indw,:) = Smat(indw,:) + Smatjw;
            end
            %  Crossproducts of D^m with forcing terms
            if funtypej
                Smatji = inprod_basis_Data2LD(Xbasisi,  Ufdj, ...
                                                      1, coefStrj, ...
                                                      order,    0);
            else
                BAtenswj = modelStructi.BAtens{nXtermi+1,jforce};
                Smatji = inprodijw(nXbasisi, 1, ...
                                   nUbasisj, nAbasisj, nrep, ...
                                   1, Avecj, Ucoefj, BAtenswj);
            end
            Smatji = weighti.*factorj.*rhoVec(ivar).*Smatji/T;
            Smat(ind,:) = Smat(ind,:) - Smatji;
        end  %  end of forcing term loop for Smat
    end  %  end of if nforce(ivar) > 0 loop for Smat
end  %  end of variable loop for Smat

%--------------------------------------------------------------------------
%    Compute the gradient DSarray(theta) of the penalty vector or matrix 
%--------------------------------------------------------------------------

if nargout > 1
    if nrep == 1
        % no replications of forcing terms, DSarray is a matrix
        % with columns corresponding to parameters theta
        DSarray = zeros(ncoefsum,nthetaH+nthetaF);
    else
        % DSarray is a 3-D array with second dimension containing theta
        DSarray = zeros(ncoefsum,nthetaH+nthetaF,nrep);
    end
    %  ---------------  Computation of DASarray -------------------
    %  partial derivatives of product of homogeneous terms
    %  and forcing terms with respect to all non-fixed forcing coefficients
    m2 = 0;
    for ivar=1:nvar
        m1   = m2 + 1;
        m2   = m2 + ncoefvec(ivar);
        modelStructi = modelCell{ivar};
        indi = m1:m2;
        weighti  = modelStructi.weight;
        nXbasisi = ncoefvec(ivar);
        nXtermi  = modelStructi.nallXterm;
        nFtermi  = modelStructi.nallFterm;
        Xbasisi  = XbasisCell{ivar};
        order    = modelStructi.order;
        n2 = nthetaH;
        for jforce = 1:nFtermi
            modelStructij = modelStructi.FCell{jforce};
            %  get details for forcing term
            [coefStrj, factorj, Ufdj, Avecj, Ucoefj, nUbasisj, Aestimj, ...
                nAbasisj, funtypej] = ...
                getForceTerm(modelStructij, coefCell);
            if Aestimj
                %  The index set for the parameters for this forcing
                %  coefficient
                n1 = n2 + 1;
                n2 = n2 + length(Avecj);
                indtha  = n1:n2;
                %  partial derivatives of products of homogeneous terms
                %  with forcing terms wrt forcing coefficients
                for iw=1:nXtermi
                    modelStructiw = modelStructi.XCell{iw};
                    [coefStrw, ivw, factorw, Bvecw, nWbasisw, ...
                        ~, nderivw, funtypew] = ...
                        getHomoTerm(modelStructiw, coefCell);
                    indw = (ncoefcum(ivw)+1):ncoefcum(ivw+1);
                    nXbasisw = ncoefvec(ivw);
                    %  compute the matrix or array of products of
                    %  partial derivatives of these forcing
                    %  coefficients and homogeneous basis functions
                    if funtypew || funtypej
                        % one or more user-defined coefficients
                        DASarrayjw = ...
                            inprod_Dbasis_Data2LD(Ufdj, Xbasisw,  ...
                            coefStrj, coefStrw, 0, nderivw);
                        DASarrayjw = permute(DASarrayjw,[2,3,1]);
                    else
                        BAtenswj = modelStructi.BAtens{iw,jforce};
                        DASarrayjw  = inprodjw(nXbasisw, nWbasisw, ...
                            nUbasisj, nAbasisj, nrep, Bvecw, ...
                            Ucoefj, BAtenswj);
                    end
                    DASarrayjw = weighti.*factorj.*factorw* ...
                        rhoVec(ivar).*DASarrayjw(:,Aestimj,:)./T;
                    if nrep == 1
                        %  disp([indtha,size(DASarrayjw)])
                        DSarray(indw,indtha) = ...
                            DSarray(indw,indtha) + DASarrayjw;
                    else
                        for irep=1:nrep
                            DSarray(indw,indtha,irep) = ...
                                DSarray(indw,indtha,irep) + ...
                                DASarrayjw(:,:,irep);
                        end
                    end
                end  %  end of homogeneous term loop for DASarray
                %  partial derivatives of cross-products of D^m
                %  with forcing terms wrt forcing coefficients
                if funtypej
                    DASarrayji = ...
                        inprod_Dbasis_Data2LD(Ufdj,  Xbasisi,  ...
                        coefStrj, coefStri, 0, order);
                    DASarrayji = permute(DASarrayji,[2,3,1]);
                else
                    BAtensij = modelStructi.BAtens{nXtermi+1,jforce};
                    DASarrayji  = inprodjw(nXbasisi, 1, ...
                        nUbasisj, nAbasisj, ...
                        nrep, 1, Ucoefj, BAtensij);
                end
                DASarrayji = weighti.*factorj.*rhoVec(ivar).* ...
                    DASarrayji/T;
                if nrep == 1
                    DSarray(indi,indtha) = ...
                        DSarray(indi,indtha) - DASarrayji;
                else
                    DSarray(indi,indtha,:) = ...
                        DSarray(indi,indtha,:) - DASarrayji;
                end
            end  %  end of if Aestimj for DASarray
        end  %  end of forcing term loop for DASarray row ivar
    end  %  end of first variable loop for DSarray
            
    %  ---------------  Computation of DBSarray -------------------
    %  partial derivatives of product of homogeneous terms and
    %  forcing terms with respect to all non-fixed homogeneous terms.
    l2 = 0;
    for ivar=1:nvar
        modelStructi = modelCell{ivar};
        weighti  = modelStructi.weight;
        nXtermi  = modelStructi.nallXterm;
        nFtermi  = modelStructi.nallFterm;
        %  Crossproducts of homogeneous terms with forcing terms,
        %  derivative with respect to homogeneous coefficient
        for iw=1:nXtermi
            modelStructiw = modelStructi.XCell{iw};
            %  get details for homogeneous term iw
            [coefStrw, ivw, factorw, Bvecw, nWbasisw, ...
                Bestimw, nderivw, funtypew] = ...
                getHomoTerm(modelStructiw, coefCell);
            nXbasisw = ncoefvec(ivw);
            Xbasisw  = XbasisCell{ivw};
            if Bestimw
                indw   = (ncoefcum(ivw)+1):ncoefcum(ivw+1);
                l1 = l2 + 1;
                l2 = l2 + length(Bvecw);
                indthw = l1:l2;
                for jforce = 1:nFtermi
                    modelStructij = modelStructi.FCell{jforce};
                    %  get details for forcing term jforce
                    [coefStrj, factorj, Ufdj, Avecj, Ucoefj, ...
                        nUbasisj, ~, nAbasisj, funtypej] = ...
                        getForceTerm(modelStructij, coefCell);
                    if funtypej || funtypew
                        % one or more user-defined coefficients
                        DBSarrayjw = ...
                            inprod_Dbasis_Data2LD(Xbasisw,  Ufdj, ...
                            coefStrw, coefStrj, nderivw,  0);
                        DBSarrayjw = permute(DBSarrayjw,[1,3,2]);
                    else
                        BAtenswj   = modelStructi.BAtens{iw,jforce};
                        DBSarrayjw = inprodwj(nXbasisw, nWbasisw, ...
                            nUbasisj, nAbasisj, ...
                            nrep, Avecj, Ucoefj, ...
                            BAtenswj);
                    end
                    DBSarrayjw = weighti.*factorw.*factorj.* ...
                        rhoVec(ivar).*DBSarrayjw./T;
                    if nrep == 1
                        DSarray(indw,indthw) = ...
                            DSarray(indw,indthw) + DBSarrayjw;
                    else
                        DSarray(indw,indthw,:) = ...
                            DSarray(indw,indthw,:) + DBSarrayjw;
                    end
                end  %  end of forcing term loop for DBSarray
            end  %  end of if Bestimj for DBSarray
        end  %  end of homogeneous term loop for DSarray row ivar
    end  %  end of second variable loop for DSarray
end  

%  ------------------------------------------------------------------------

function Smatjw = inprodijw(nXbasisw, nWbasisw, nUbasisj, nAbasisj, ...
                            nrep, Bvecw, Avecj, Ucoefj, BAtenswj)
ncum = cumprod([nAbasisj, nUbasisj, nWbasisw, nXbasisw]);
Smatjw = zeros(nXbasisw,nrep);
for irep=1:nrep
    for i=1:nXbasisw
        for j=1:nWbasisw
            for k=1:nUbasisj
                for l=1:nAbasisj
                    ijkl = (i-1)*ncum(3) + ...
                           (j-1)*ncum(2) + ...
                           (k-1)*ncum(1) + l;
                    Smatjw(i,irep) = Smatjw(i,irep) + ...
                        Bvecw(j)*Avecj(l)*BAtenswj(ijkl)*Ucoefj(k,irep);
                end
            end
        end
    end
end

%  ------------------------------------------------------------------------

function DASarray = inprodjw(nXbasisw, nWbasisw, nUbasisj, nAbasisj, ...
                             nrep, Bvecw, Ucoefj, BAtenswj)
ncum = cumprod([nAbasisj, nUbasisj, nWbasisw, nXbasisw]);
if nrep == 1
    DASarray = zeros(nXbasisw,nAbasisj);
else
    DASarray = zeros(nXbasisw,nAbasisj,nrep);
end
for irep=1:nrep
    for i=1:nXbasisw
        for j=1:nWbasisw
            for k=1:nUbasisj
                for l=1:nAbasisj
                    ijkl = (i-1)*ncum(3) + ...
                           (j-1)*ncum(2) + ...
                           (k-1)*ncum(1) + l;
                    if nrep == 1
                        DASarray(i,l) = DASarray(i,l) + ...
                            Bvecw(j).*BAtenswj(ijkl).*Ucoefj(k);
                    else
                        DASarray(i,l,irep) = DASarray(i,l,irep) + ...
                            Bvecw(j).*BAtenswj(ijkl).*Ucoefj(k,irep);
                    end
                end
            end
        end
    end
end

%  ------------------------------------------------------------------------

function DBSarray = inprodwj(nXbasisw, nWbasisw, nUbasisj, nAbasisj, ...
                            nrep, Avecj, Ucoefj, BAtenswj)
ncum = cumprod([nAbasisj, nUbasisj, nWbasisw, nXbasisw]);
if nrep == 1
    DBSarray = zeros(nXbasisw,nWbasisw);
else
    DBSarray = zeros(nXbasisw,nWbasisw,nrep);
end
for irep=1:nrep
    for i=1:nXbasisw
        for j=1:nWbasisw
            for k=1:nUbasisj
                for l=1:nAbasisj
                    ijkl = (i-1)*ncum(3) + ...
                           (j-1)*ncum(2) + ...
                           (k-1)*ncum(1) + l;
                    if nrep == 1
                        DBSarray(i,j) = DBSarray(i,j) + ...
                            Avecj(l).*BAtenswj(ijkl).*Ucoefj(k);
                    else
                        DBSarray(i,j,irep) = DBSarray(i,j,irep) + ...
                            Avecj(l).*BAtenswj(ijkl).*Ucoefj(k,irep);
                    end
                end
            end
        end
    end
end


