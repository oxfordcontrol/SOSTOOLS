function varargout=p2s(varargin)
% function s=p2s(p)
%
% DESCRIPTION 
%   Converts from a multipoly polynomial to a symbolic toolbox polynomial.
%   
% INPUTS 
%   p: Polynomial created using the multipoly toolbox
%
% OUTPUTS  
%   s: Polynomial created using the symbolic math toolbox  
%  
% SYNTAX 
%   s = p2s(p)
%
% See also s2p
  
% 1/30/2003  PJS  Initial Coding  
% 11/23/2010 PJS  Updated for matrix polynomials


for i=1:nargin

% Get polynomial info
[nr nc] = size(varargin{i});
[nt,nv] = size(varargin{i}.degmat);
coef = varargin{i}.coefficient;
deg = varargin{i}.degmat;
var = varargin{i}.varname;

% Create symbolic toolbox variables
for i1 = 1:nv
    eval(['syms ' var{i1} ' real']);
end

% Construct monomial by monomial
varargout{i} = zeros(nr,nc);
for i1 = 1:nt
    termexp = deg(i1,:);
    idx = find(termexp);
    term = 1;
    for i2 = 1:length(idx)
        eval(['term = term*(' var{idx(i2)} ')^termexp(idx(i2));']);
    end
    varargout{i} = varargout{i} + reshape(coef(i1,:),[nr,nc])*term;
end  
end