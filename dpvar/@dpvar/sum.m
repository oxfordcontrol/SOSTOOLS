function Dout = sum(Din,dir)
% Dout = sum(Din,dir) computes the trace sum of all elements of a
% dpvar object Din along direction dir.
% 
% INPUTS:
% Din: dpvar object of size mxn
% dir: positive integer scalar or string optional input (defaults to 1)
%       dir = 1: compute sum over all rows of Din
%       dir = 2: compute sum over all columns of Din
%       dir = [1,2], ':', 'all': compute sum over all rows of Din
% 
% OUTPUTS:
% Dout: dpvar object of size 1xn if dir=1, mx1 if dir=2, or 1x1 if 
%       dir='all', corresponding to the sum of elements of Din along dir.
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
% Initial coding DJ - 02/15/2022

% Check that a direction is properly specified
if nargin==1
    dir = 1;
elseif strcmp(dir,':') || strcmp(dir,'all')
    dir = [2,1];    % maybe slightly cheaper than [1,2]
elseif ~isnumeric(dir) || any(dir~=1 & dir~=2)
    error('Dimension argument must be either 1, 2, or ''all''.')
end

% Avoid redundant sums, but don't change order of dir, just in case
dir = unique(dir,'stable');

% Exclude cases where summation is not necessary
[m,n] = size(Din);
if isempty(Din) || all([m,n]==[1,1]) || (m==1 && all(dir==1)) || (n==1 && all(dir==2))
    Dout = compress(Din);
    return
end

% Take the sum along each direction
Dout = combine(Din);
for j=1:length(dir)
   
    [m,n] = size(Dout);
    nd = length(Dout.dvarname)+1;       % Number of decision vars + constant term
    nz = max(size(Dout.degmat,1),1);    % Number of monomials (=1 if no monomials)
    
    if dir(j)==1       % Sum along the rows
        P = repmat(speye(nd),1,m);
        Cn = P*Dout.C;              % Sum coefficients associated to same dvar, column, and monomial
    
        Dout = combine(dpvar(Cn,Dout.degmat,Dout.varname,Dout.dvarname,[1,n]),'extended');
        
    elseif dir(j)==2   % Sum along the colums
        
        C = reshape(Dout.C,m*nd*nz,n);  % Place coefficients associated to each column side-by-side
        Cn = sum(C,2);                  % Sum over the columns
        Cn = reshape(Cn,m*nd,nz);       % Reshape to appropriate size
        
        Dout = compress(dpvar(Cn,Dout.degmat,Dout.varname,Dout.dvarname,[m,1]));
        
    end  
end
end