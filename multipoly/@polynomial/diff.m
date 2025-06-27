function b = diff(a,x,d)
% function B=diff(A,X,D)
%
% DESCRIPTION
%   Element-by-element differentiation of a polynomial with respect
%   to a single variable up to order d.
%
% INPUTS
%   A: polynomial
%   X: Differentiate with respect to the (single) variable X.
%   D: Order of the derivative to be taken.
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B = diff(A,X,D);
%     Differentiate the polynomial, A, with respect to X, up to order D.
%     A should be a polynomial and X should be a polynomial variable or a 
%     string. Differentiation is done element-by-element if A is a matrix.
%
% EXAMPLE
%   pvar x y z;
%   f = 2*x^3+5*y*z-2*x*z^2;
%   df = diff(f,x,2)
%
% See also: jacobian, int

% 11/5/2002 PJS  Initial Coding
% 10/21/2010 PJS diff calls jacobian (reduces redundant code)
% 06/27/2025, DJ: Add support for higher-order derivatives.

% Error Checking
if ~isa(a,'polynomial')
    error("First input argument must be of type 'polynomial'.")
end
if nargin<=1
    if isscalar(a.varname)
        x = a.varname{1};
    else
        error("Not enough input arguments.")
    end
    d = 1;
elseif nargin<=2
    if isnumeric(x)
        if isscalar(a.varname) 
            d = x;
            x = a.varname{1};
        else
            error("Second argument must be a single polynomial variable or a string.")
        end
    else
        d = 1;
    end
elseif nargin>3
    error("Too many input arguments.")
end

if ispvar(x) && length(x)==1
    x = char(x);
    x = x{1};
elseif ~ischar(x)
    error('X must be a single polynomial variable or a string');
end

if ~isnumeric(d) || ~isscalar(d) || d<0 || round(d)~=d
    error("Order of derivative must be specified as nonnegative integer.")
end

% Call jacobian to perform differentiation
% jacobian works on a column vector and returns a column vector
sza = size(a);
b = a(:);
for iter=1:d                                                                % 06/27/2025, DJ
    b = jacobian(b,x);
end

% Reshape derivative back to size of a
b = reshape(b,sza);

