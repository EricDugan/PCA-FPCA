function [Rmat, DRarray] = Data2LD_R(XbasisCell, modelCell, coefCell, ...
                                      rhoVec, ntheta)
%  Data2LD ... stands for "Data to Linear Dynamics"
%  Data2LD_R computes the penalty matrix R associated with the homogeneous
%  portion of a linear differential operator L as well as its partial
%  derivative with respect to parameters defining the homogeneous portion.
%  For a single variable whose approximated in terms of an exansion in
%  terms of a vector \phi of basis functions, R is
%                 R = \int [L \phi(t)] [L \phi(t)]' dt.
%  R is of order K, the number of basis functions, symmetric, and of rank
%  K - m where m is the order of the largest derivative in the operator.
%  The eigenanalysis of R can be used to construct an alternative basis
%  expansion defined in terms an increasing order of complexity of shape.
%
%  If multiple variables are involved, then R is a composite matrix
%  constructed from inner products and cross-products of the basis 
%  function vectors associate with each variable.  It's order will be
%  \sum K_i.
%
%  This version approximates the integrals in the penalty terms by using 
%  inprod_basis to compute the cross-product matrices for the  
%  \beta-coefficient basis functions and the corresponding derivative of 
%  the x-basis functions,and the cross-product matrices for the 
%  \alpha-coefficients and the corresponding U functions.  
%  These are computed upon the first call to Data2LD4, and then retained 
%  for subsequent calls by using the persistent command.  See lines about 
%  560 to 590 for this code.
%
%  This version disassociates coefficient functions from equation 
%  definitions to allow some coefficients to be used repeatedly and for
%  both homogeneous and forcing terms.  It requires an extra argument
%  COEFCELL that contains the coefficients and the position of their
%  coefficient vectors in vector THETA.
%
%  Arguments:
%
%  XBASISCELL ... A functional data object or a BASIS object.  If so, the 
%               smoothing parameter LAMBDA is set to 0.
%
%  MODELCELL...  A cell aray of length NVAR. Each cell contains a 
%                struct object with members:              
%                XCell ... cell array of length number of homogeneous terms
%                          Each cell contains a struct object with members:
%                          WfdPar     ... fdPar object for the coefficient
%                          variable   ... the index of the variable
%                          derivative ... the order of its derivative
%                          ncoef      ... if coefficient estimated, its 
%                                         location in the composite vector 
%                          factor     ... a scalar multiplier (def. 1)
%                FCell ... cell array of length number of forcing terms
%                          Each cell contains a struct object with members:
%                          AfdPar ... an fdPar object for the coefficient
%                          Ufd    ... an fd object for the forcing function
%                          ncoef... if coefficient estimated, its location
%                                   in the composite vector 
%                order     ... the highest order of derivative
%                name      ... a  tag for the variable
%                nallXterm ... the number of homogeneous terms
%                nallFterm ... the number of forcing functions
%  COEFCELL  ... A cell array of length NCOEF containing struct objects
%                with fields:
%               parvec   ... a vector of parameters
%               coeftype ... homogeneous or forcing
%               estimate ... 0, held fixed, otherwise, estimated 
%               fun      ... functional basis, fd, or fdPar object, 
%                            or a struct object for a general function 
%                            with fields:
%                 fd      ... function handle for evaluating function
%                 Dfd     ... function handle for evaluating 
%                             partial derivative with respect to parameter
%                 more    ... object providing additional information for 
%                             evaluating coefficient function
%  RHOVEC    ... A vector of length NVAR containing values in [0,1].  
%                The data sums of squares are weighted by rho and 
%                the roughness penalty by 1-rho.
%
%  BTENSORCELL 

%  Last modified 27 March 2020

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

%  ------------------------------------------------------------------------
%                  Compute penalty matrix Rmat(theta)
%  ------------------------------------------------------------------------

%  The order of Rmat is the sum of the number of numbers of spline
%  coefficients for each variable or trajectory

Rmat = zeros(ncoefsum,ncoefsum);

%  loop through variables or equations

m2 = 0;
for ivar=1:nvar
    %  select struct object for this variable that defines this variable's
    %  structure.
    modelStructi = modelCell{ivar};
    %  weight for this variable for computing fitting criterion
    weighti = modelStructi.weight;
    %  indices in Rmat for this order nXbasisi submatrix
    nXbasisi = ncoefvec(ivar);
    m1  = m2 + 1;
    m2  = m2 + nXbasisi;
    indi = m1:m2;
    %  index of homogeneous terms
    nXtermi  = modelStructi.nallXterm;
    %  order of the equation
    order    = modelStructi.order;
    %  the functional basis object
    Xbasisi  = XbasisCell{ivar};
    %  First we compute the term in Rmat for the left side of the equation
    Btensii = modelStructi.Btens{nXtermi+1,nXtermi+1};
    Rmatii  = inprodwx(nXbasisi, 1, nXbasisi, 1, 1, 1, Btensii);
    %  multiply by weight, value of smoothing parameter rho, and divide by
    %  duration of interval
    Rmatii   = weighti*rhoVec(ivar)*Rmatii/T;
    %  initialize the left side submatrix withinin supermatrix Rmat
    Rmat(indi,indi) = Rmat(indi,indi) + Rmatii;
    %  now we compute the contribution to R mat for each pair of 
    %  homogeneous terms in the equation
    for iw=1:nXtermi
        %  select homogeneous term iw
        modelStructiw = modelStructi.XCell{iw};
        %  get the details for homogeneous coefficient iw
        [coefStrw, ivw, factorw, ...
                 Bvecw, nWbasisw, ~, derivw, funtypew] = ...
                          getHomoTerm(modelStructiw, coefCell);
        indw     = (ncoefcum(ivw)+1):ncoefcum(ivw+1);
        nXbasisw = ncoefvec(ivw);
        Xbasisw  = XbasisCell{ivw};
        for ix=1:nXtermi
            %  select homogeneous term ix
            modelStructix = modelStructi.XCell{ix};
            %  get the details for homogeneous coefficient ix
            [coefStrx, ivx, factorx, ...
                 Bvecx, nWbasisx, ~, derivx, funtypex] = ...
                          getHomoTerm(modelStructix, coefCell);
            %  define coefficient of active variable and
            %  it's derivative index
            indx     = (ncoefcum(ivx)+1):ncoefcum(ivx+1);
            nXbasisx = ncoefvec(ivx);
            Xbasisx  = XbasisCell{ivx};
            %  set up the inner products for the two basis systems
            if funtypew || funtypex
                %  the user-supplied coefficient function case ...
                %  must be computed anew for each pair
                Rmatwx = inprod_basis_Data2LD(Xbasisw,  Xbasisx, ...
                                              coefStrw, coefStrx, ...
                                              derivw,   derivx);
            else
                %  the much faster case where both homogeneous coefficients
                %  are B-spline functional data objects
                Btenswx  = modelStructi.Btens{iw,ix};
                %  reformat the vector of values
                Rmatwx   = inprodwx(nXbasisw, nWbasisw, ...
                                    nXbasisx, nWbasisx, ...
                                    Bvecw, Bvecx, Btenswx);
            end
            %  apply the scale factors
            Rmatwx   = weighti*factorw*factorx*rhoVec(ivar)*Rmatwx/T;
            %  increment the submatrix
            Rmat(indw,indx) = Rmat(indw,indx) + Rmatwx;
        end
        %  now we need to compute the submatrices for the product of
        %  the left side of the equation and each homogeneous term in
        %  the equation
        if funtypew
            %  user-supplied coefficient case
            Rmatiw = inprod_basis_Data2LD(Xbasisi, Xbasisw, ...
                                          coefStri, coefStrw, ...
                                          order,    derivw);
        else
            %  reformat the vector of previous computed values
            Btensiw = modelStructi.Btens{nXtermi+1,iw};
            Rmatiw  = inprodwx(nXbasisi, 1, nXbasisw, nWbasisw, ...
                1, Bvecw, Btensiw);
        end
        %  apply the scale factors
        Rmatiw  = weighti*factorw*rhoVec(ivar)*Rmatiw/T;
        %  subtract the increment for both off-diagonal submatrices
        Rmat(indi,indw) = Rmat(indi,indw) - Rmatiw;
        Rmat(indw,indi) = Rmat(indw,indi) - Rmatiw';
    end
end

%  ------------------------------------------------------------------------
%  Compute partial derivatives with respect to the parameters in vector
%  theta that are to be estimated if they are required
%  This is a three-dimensional array DRarray with the final dimension being
%  The verbose commentary above will not be continued.
%  ------------------------------------------------------------------------

if nargout > 1
    
%  set up DRarray

DRarray = zeros(ncoefsum,ncoefsum,ntheta);

m2 = 0;
for ivar=1:nvar
    modelStructi = modelCell{ivar};
    weighti = modelStructi.weight;
    m1   = m2 + 1;
    m2   = m2 + ncoefvec(ivar);
    indi = m1:m2;
    nXtermi  = modelStructi.nallXterm;
    nXbasisi = ncoefvec(ivar);
    Xbasisi  = XbasisCell{ivar};
    nderivi  = modelStructi.order;
    %  loop through homogeneous terms using their derivatives
    for iw=1:nXtermi
        modelStructiw = modelStructi.XCell{iw};
        %  get the details for homogeneous coefficient iw
        [coefStrw, ivw, factorw, ...
                 ~, nWbasisw, Bestimw, derivw, funtypew] = ...
                          getHomoTerm(modelStructiw, coefCell);
        %  disp([ivar,iw])
        %  disp(coefStrw.estimate)
        %  disp(Bestimw)
        if any(Bestimw)
            %  define coefficient of estimated variable and
            %  it's derivative index
            indw     = (ncoefcum(ivw)+1):ncoefcum(ivw+1);
            indthw   = coefStrw.index;
            nXbasisw = ncoefvec(ivw); 
            Xbasisw  = XbasisCell{ivw};
            %  loop through all terms using their values rather than their
            %  derivatives
            for ix=1:nXtermi
                modelStructix = modelStructi.XCell{ix};
                %  get the details for homogeneous coefficient ix
                [coefStrx, ivx, factorx, ...
                    Bvecx, nWbasisx, ~,  derivx, funtypex] = ...
                    getHomoTerm(modelStructix, coefCell);
                
                indx     = (ncoefcum(ivx)+1):ncoefcum(ivx+1);
                nXbasisx = ncoefvec(ivx);
                Xbasisx  = XbasisCell{ivx};
                %  get the tensor vector for this pair of coefficients
                %  and derivatives
                if funtypew || funtypex
                    %  user-supplied coefficient functions
                    DRarraywx = ...
                        inprod_Dbasis_Data2LD(Xbasisw,  Xbasisx,  ...
                                              coefStrw, coefStrx, ...
                                              derivw,  derivx);
                else
                    %  B-spline expansions
                    Btenswx   = modelStructi.Btens{iw,ix};
                    DRarraywx = ...
                        inprodDwx(nXbasisw, nWbasisw, ...
                                  nXbasisx, nWbasisx, ...
                                  Bvecx,    Btenswx);
                end
                %  rescale the inner product
                DRarraywx = weighti*factorw*factorx*rhoVec(ivar)* ...
                            DRarraywx(:,:,Bestimw)/T;
                %  increment the inner product and its transpose for
                %  the appropriate location in DRarray
                %  diagonal case
                if iw == ix 
                    DRarray(indw,indw,indthw(Bestimw)) = ...
                    DRarray(indw,indw,indthw(Bestimw)) + 2*DRarraywx;
                else
                    %  lower left off-diagonal
                    DRarray(indw,indx,indthw(Bestimw)) = ...
                    DRarray(indw,indx,indthw(Bestimw)) + ...
                    DRarraywx;
                    %  upper right off-diagonal, transpose required
                    DRarray(indx,indw,indthw(Bestimw)) = ...
                    DRarray(indx,indw,indthw(Bestimw)) + ...
                    permute(DRarraywx,[2,1,3]);
                end
            end
            %  partial derivatives wrt Wcoef for cross-products with D^m
            %  here x = ivar, Wbasisx is the constant basis, and
            %  Bvecx = 1;
            %  get the tensor vector for this pair of coefficients
            %  and derivatives
            %  compute inner product with respect to the ix dimension
            if funtypew
                % user code
                DRarraywi = inprod_Dbasis_Data2LD(Xbasisw,  Xbasisi,  ...
                                                  coefStrw, coefStri, ...
                                                  derivw,  nderivi);
            else
                %  B-spline expansion code
                Btenswi   = modelStructi.Btens{iw,nXtermi+1};                
                DRarraywi = inprodDwx(nXbasisw, nWbasisw, nXbasisi, 1, ...
                                      1, Btenswi);               
            end
            %  rescale the inner product
            DRarraywi = weighti*factorw*rhoVec(ivar)* ...
                        DRarraywi(:,:,Bestimw)/T;
            %  decrement the inner product and its transpose for
            % the appropriate location in DRarray
            if nXbasisw == 1
                DRarray(indw,indi,indthw(Bestimw)) = ...
                DRarray(indw,indi,indthw(Bestimw)) - ...
                    DRarraywi;
                DRarray(indi,indw,indthw(Bestimw)) = ...
                DRarray(indi,indw,indthw(Bestimw)) - ...
                    DRarraywi';
            else
                DRarray(indw,indi,indthw(Bestimw)) = ...
                DRarray(indw,indi,indthw(Bestimw)) - ...
                    DRarraywi;
                DRarray(indi,indw,indthw(Bestimw)) = ...
                DRarray(indi,indw,indthw(Bestimw)) - ...
                    permute(DRarraywi,[2,1,3]);
            end
        end
    end
end

else
    
    DRarray = [];
end

%  ------------------------------------------------------------------------
%  These are the two reformating functions that are required.
%  For large scale problems it would be better to supply these in 
%  C,  C++, or etc. because Matlab does not handle orders of looping
%  beyond 2 efficiently.

function Rmatwx = inprodwx(nXbasisw, nWbasisw, nXbasisx, nWbasisx, ...
                           Bvecw, Bvecx, Btenswx)
Rmatwx = zeros(nXbasisw,nXbasisx);
ncum   = cumprod([nWbasisx, nXbasisx, nWbasisw, nXbasisw]);
for i=1:nXbasisw
    for k=1:nXbasisx
        for j=1:nWbasisw
            for l=1:nWbasisx
                ijkl = (i-1)*ncum(3) + ...
                       (j-1)*ncum(2) + ...
                       (k-1)*ncum(1) + l;
                Rmatwx(i,k) = Rmatwx(i,k) + ...
                    Btenswx(ijkl)*Bvecw(j)*Bvecx(l);
            end
        end
    end
end

%  ------------------------------------------------------------------------

function  DRarraywx = inprodDwx(nXbasisw, nWbasisw, nXbasisx, nWbasisx, ...
                               Bvecx, Btenswx)
DRarraywx = zeros(nXbasisw, nXbasisx, nWbasisw);
ncum = cumprod([nWbasisx, nXbasisx, nWbasisw, nXbasisw]);
for i=1:nXbasisw
    for j=1:nWbasisw
        for k=1:nXbasisx
            for l=1:nWbasisx
                ijkl = (i-1)*ncum(3) + ...
                       (j-1)*ncum(2) + ...
                       (k-1)*ncum(1) + l;
                if abs(ijkl) > eps
                    DRarraywx(i,k,j) = DRarraywx(i,k,j) + ...
                        Bvecx(l)*Btenswx(ijkl);
                else
                    error('abs(ijkl) <= eps in DRarraywx.');
                end
            end
        end
    end
end


