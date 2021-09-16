function varargout=s2p(varargin)
% function p=s2p(s)
%
% DESCRIPTION
%   Converts from a symbolic toolbox polynomial to a multipoly polynomial.
%
% INPUTS
%   s: Polynomial created using the symbolic toolbox
%
% OUTPUTS
%   p: Polynomial created using the multipoly toolbox
%
% SYNTAX
%   p = s2p(s)
%
% See also p2s

% copied and modified from s2p() by PJS to work with newer matlab versions,
% 9/2/2021 SS

for i=1:nargin
if ~isa(varargin{i},'sym')
    error('Input must be a symbolic toolbox object');
end

% Convert to variable precision arithmetic
varargin{i} = vpa(varargin{i});

% Create pvars for each variable in s
vars = symvar(varargin{i});

if isempty(vars)
    nv = 0;
else
    nv = length(vars);
end
for i1 =1:nv
    varname = char(vars(i1));
    pvar(varname);
end

% Expand the polynomial 
se = expand(varargin{i});

% Evaluate in symbolic expression to convert to multipoly
varargout{i} = eval(se);
end
