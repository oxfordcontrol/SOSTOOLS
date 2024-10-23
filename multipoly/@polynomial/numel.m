function N = numel(P)
% N = NUMEL(P) determines the number of elements of the matrix-valued 
% 'polynomial' object P 
% 
% INPUTS:
% P:    mxn object of type 'polynomial';
% 
% OUTPUTS:
% N:    'double' computed as m*n, providing total number of elements in the
%       object P;
% 
%
% Initial coding DJ - 10/23/2024

N = prod(P.matdim);

end