function inprodvec = Inner_Loop(in1, in2, in3, in4, wvec)
%InnerLoop Summary of this function goes here
%   Detailed explanation goes here

% output_args = Inner_Loop.mexmaci64(in1, in2, in3, in4, wvec)';
% output_args = loop_mex(in1', in2', in3', in4');

ncol1 = size(in1,2);
ncol2 = size(in2,2);
ncol3 = size(in3,2);
ncol4 = size(in4,2);

nrow = length(wvec);
if size(in1,1) ~= nrow || size(in2,1) ~= nrow || ...
   size(in2,1) ~= nrow || size(in4,1) ~= nrow 
    error('number of rows do not match')
end

inprodvec = zeros(1,ncol1*ncol2*ncol3*ncol4);
for i=1:nrow
    m1 = 1;
    m2 = ncol1*ncol2;
    prod = wvec(i).*(in1(i,:)'*in2(i,:));
    prod = prod(:)';
    inprodvec(m1:m2) = inprodvec(m1:m2) + prod;
    m1 = m2 + 1;
    m2 = m2 + ncol2*ncol3;
    prod = wvec(i).*(in2(i,:)'*in3(i,:));
    prod = prod(:)';
    inprodvec(m1:m2) = inprodvec(m1:m2) + prod;
    m1 = m2 + 1;
    m2 = m2 + ncol3*ncol4;
    prod = wvec(i).*(in3(i,:)'*in4(i,:));
    prod = prod(:)';
    inprodvec(m1:m2) = inprodvec(m1:m2) + prod;
end