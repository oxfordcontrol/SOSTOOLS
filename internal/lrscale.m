function X = lrscale(X,L,R)
%LRSCALE  Applies left and right scaling matrices.
%
%   Y = LRSCALE(X,L,R) forms Y = diag(L) * X * diag(R) in 
%   2mn flops if X is m-by-n.  L=[] or R=[] is interpreted
%   as the identity matrix.

%   Copyright 1986-2011 The MathWorks, Inc.
cL = size(L,2);
if cL==0
   L = ones(size(X,1),1);
elseif cL>1
   % Make L a column vector
   L = L.';
end
rR = size(R,1);
if rR==0
   R = ones(1,size(X,2));
elseif rR>1
   % Make R a row vector
   R = R.';
end
% Note (9/11): This implementation is on par with C++ implementation
% and faster than BSXFUN for small to midsize arrays.
X = X .* (L * R);