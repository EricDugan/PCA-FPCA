function inprodvec = ...
    Inner_Loop_Matlab(basismat1, basismat2, basismat3, basismat4, wvec)
%InnerLoop Summary of this function goes here
%   Detailed explanation goes here

% output_args = Inner_Loop.mexmaci64(basismat1, basismat2, basismat3, basismat4, wvec)';
% output_args = loop_mex(basismat1', basismat2', basismat3', basismat4');

ncol1 = size(basismat1,2);
ncol2 = size(basismat2,2);
ncol3 = size(basismat3,2);
ncol4 = size(basismat4,2);

nrow = length(wvec);
if size(basismat1,1) ~= nrow || size(basismat2,1) ~= nrow || ...
   size(basismat2,1) ~= nrow || size(basismat4,1) ~= nrow 
    error('number of rows do not match')
end

inprodvec = zeros(1,ncol1*ncol2*ncol3*ncol4);

for i=1:ncol1
    for j=1:ncol2
        for k=1:ncol3
            for l=1:ncol4
                for m=1:nrow
                    index = (i-1)*ncol4*ncol3*ncol2 ...
                          + (j-1)*ncol4*ncol3 ...
                          + (k-1)*ncol4 ...
                          + l;
                    inprodvec(index) = inprodvec(index) + ...
                        basismat1(m,i)*basismat2(m,j)* ...
                        basismat3(m,k)*basismat4(m,l)*wvec(m);
                end
            end
        end
    end
end
                      
