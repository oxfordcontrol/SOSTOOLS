function b = int_p(a,x,x1,x2)
% b = int_p(a,x,x1,x2) integrates the polynomial a in independent variable x
% from lower limit x1 to upper limit x2
% 
% INPUTS:
% a: polynomial to integrate
% x: polynomial variable or string specifying var to integrate with respect to
% x1: lower limit of definite integral
% x2: upper limit of definite integral
% 
% OUTPUTS:
% b: integrated polynomial object
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% 6/20/2010: MMP  Initial Coding   -   based on diff.m by PJS
% 6/16/2013: MMP  Modified to eliminate for-loops


% Error Checking
if nargin~=4
    error('Incorrect number of input arguments. Format int_p(f, intVar, lowerlimit, upperlimit)');
end

if isa(x,'polynomial')
    x = x.varname{1};
elseif ~ischar(x)
    error('X must be a polynomial variable or a string');
end

% Get polynomial info about a
a = polynomial(a);
adeg = a.degmat;
avar = a.varname;
nta = size(adeg,1);
[nra,nca]=size(a);
acoef = reshape(a.coefficient,nta,nra*nca);

% Get format x1 and x2
if ~isa(x1,'double')
    x1 = polynomial(x1);
    x1 = combine(x1);
end
if ~isa(x2,'double')
    x2 = polynomial(x2);
    x2 = combine(x2);
end

% Find variable we are integrating with respect to
varnumb=find(strcmp(x,avar));  % MMP - 6.16.2013

if isempty(varnumb)
    % Constant function case
    b=a*(x2-x1);
else    
    % Integrate
    acoef = diag(1./(adeg(:,varnumb)+1))*acoef;
    adeg(:,varnumb) = adeg(:,varnumb)+1;
    a_int = polynomial(acoef,adeg,avar,[nra nca]);
    b = combine(subs_p(a_int,x,x2)-subs_p(a_int,x,x1));
end