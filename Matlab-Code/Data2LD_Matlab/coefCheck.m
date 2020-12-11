function [coefCellnew, ntheta, nparam] = coefCheck(coefCell)
%  COEFCHECK checks the cell array of coefficient functions defining a
%  possibly forced linear system of llinear differential equations.
%
%  COEFCHECK goes through the following steps:
%
%  1. Each cell is checked for being a struct object.  
%  If this test is not passed, the checking terminates.
%  
%  2. The presence and class of required fields are checked.  
%  These are:
%  Name       Class
%  parvec     real vector
%  estimate   numeric, and converted to logical
%  fun        basis, fd, fdPar, struct
%  If this test is not passed, the checking terminates.
%
%  3. for each coefStruct, the index field indicating estimated parameters 
%  is set up, with the indices of the estimated homogeneous or beta 
%  coefficients coming first.
%
%  4. for any general coefficient functions, the fields and their objects
%  are checked.
%  The required fields are:
%  Name       Class
%  fd         function_handle
%  difdip     function_handle
%  more       any object
%  

%  Last modified 27 March 2020

%  -------------------------------  Step 1  -------------------------------

%  Check that coefCell is a cell object and obtain number of coefficient
%  functions

if ~iscell(coefCell)
    error('COEFCELL is not a cell object.');
end

ncoef = length(coefCell);

%  check that each cell contains a struct object

errwrd = 0;
for icoef=1:ncoef
    coefStri = coefCell{icoef};
    if ~isstruct(coefStri)
        disp(['Cell ', num2str(icoef), ...
              ' does not contain a struct object.']);
        errwrd = 1;
    end
end

if errwrd
    error('Errors found, cannot continue.');
end

%  -------------------------------  Step 2  -------------------------------

%  check that required fields are present and that their contents
%  are in the right class

errwrd = 0;
for icoef=1:ncoef
    coefStri = coefCell{icoef};
    
    %  check that all struct objects contain a field named 'parvec'
    
    if ~isfield(coefStri, 'parvec')
        disp(['Struct object for cell ', num2str(icoef), ...
            ' does not have a parvec field.']);
        errwrd = 1;
    else
        parvec = coefStri.parvec;
        npar   = length(parvec);
        if ~isnumeric(parvec)
            disp(['Field parvec for cell ', num2str(icoef), ...
                  ' is not numeric.']);
            errwrd = 1;
        else
            parvecsize = size(parvec);
            if length(parvecsize) ~= 2
                disp(['Field parvec for cell ', num2str(icoef), ...
                      ' is not a matrix.']);
                errwrd = 1;
            elseif ~(parvecsize(1) == 1 || parvecsize(2) == 1)
                disp(['Field parvec for cell ', num2str(icoef), ...
                    ' is not a vector.']);
                errwrd = 1;
            else
                coefStri.parvec = parvec(:);
            end
        end
    end
    
    %  check that all struct objects contain a field named 'estimate',
    %  and if they do, check their lengths.
    
    if ~isfield(coefStri, 'estimate')
        disp(['Struct object for cell ', num2str(icoef), ...
            ' does not have a estimate field.']);
        errwrd = 1;
    else
        estimatei  = coefStri.estimate;
        nestimate = length(estimatei);
        if ~isnumeric(estimatei) && ~islogical(estimatei)
            disp(['Field estimate for cell ', num2str(icoef), ...
                  ' is neither numeric or logical.']);
            errwrd = 1;
        elseif nestimate ~= 1 && nestimate ~= parvecsize(1)
            disp(['Field estimate for cell ', num2str(icoef), ...
                  ' is not of length 1 or of length # parameters.']);
            errwrd = 1;
        else
            if nestimate == 1 && npar > 1
                if estimatei
                    estimatei = ones(npar,1,'logical');
                else
                    estimatei = zeros(npar,1,'logical');
                end
            else
                estimatei = logical(estimatei);
            end
            coefStri.estimate = estimatei;
        end
    end
    
    %  check that all struct objects contain a field named 'fun'
    
    if ~isfield(coefStri, 'fun')
        disp(['Struct object for cell ', num2str(icoef), ...
            ' does not have a fun field.']);
        errwrd = 1;
    else
        fun = coefStri.fun;
        if ~(isa_basis(fun) || isa_fd(fun) || isa_fdPar(fun) || ...
             isstruct(fun))
            disp(['Field fun for cell ', num2str(icoef), ...
                  ' is not a class basis, fd, fdPar, or struct.']);
            errwrd = 1;
        end
    end
    
    coefCell{icoef} = coefStri;
    
end

if errwrd
    error('Errors found, cannot continue.');
end

%  -------------------------------  Step 3  -------------------------------

%  generate index fields sequentially

ntheta = 0;
nparam = 0;
m2 = 0;
for icoef=1:ncoef
    coefStri = coefCell{icoef};
    nparami = length(coefStri.parvec);
    m1 = m2 + 1;
    m2 = m2 + nparami;
    coefStri.index = m1:m2;
    ntheta = ntheta + nparami;
    if coefStri.estimate
        nparam = nparam + nparami;
    end
    coefCell{icoef} = coefStri;
end

%  -------------------------------  Step 4  -------------------------------

%  if the fun field is a struct object, check the object

errwrd = 0;
for icoef=1:ncoef
    coefStri = coefCell{icoef};
    if isstruct(coefStri.fun)
        fdStruct = coefStri.fun;
        if isfield(fdStruct, 'fd')
            coeffd = fdStruct.fd;
            if ~strcmp(class(coeffd), 'function_handle')
                disp('Field fun for cell ', num2str(icoef), ...
                      ' does not contain a function_handle object.');
                errwrd = 1;
            end
        else
            disp('Field fd for cell ', num2str(icoef), ...
                      ' does not have a field fd.');
            errwrd = 1;
        end
        if isfield(fdStruct, 'Dfd')
            coefDfd = fdStruct.Dfd;
            if ~strcmp(class(coefDfd), 'function_handle')
                disp(['Field DfD for cell ', num2str(icoef), ...
                      ' does not contain a function_handle object.']);
                errwrd = 1;
            end
        end
        if ~isfield(fdStruct, 'more')
            disp(['Field fun for cell ', num2str(icoef), ...
                      ' does not have a field more.']);
            errwrd = 1;
        end
    end
end

if errwrd
    error('Errors found, cannot continue.');
end

coefCellnew = coefCell;



