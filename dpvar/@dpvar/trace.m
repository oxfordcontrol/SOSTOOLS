function Dout = trace(Din)
% Dout = trace(Din) computes the trace (sum of diagonal elements) of a
% dpvar object Din.
% 
% INPUTS:
% Din: dpvar object of size nxn
% 
% OUTPUTS:
% Dout: dpvar object of size 1x1 corresponding to the trace of Din
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
% Initial coding DJ - 02/08/2022

sz = size(Din);
if sz(1)~=sz(2)
    error('dpvar matrix must be square');
elseif isempty(Din)
    Dout = dpvar(0);
    return;
elseif all(sz==[1,1])
    Dout = compress(Din);
    return
else
    % Grab diagonal elements
    n = sz(1);
    nd = length(Din.dvarname)+1;      % Number of decision vars + constant term
    nz = max(size(Din.degmat,1),1);   % Number of monomials (=1 if no monomials)
    
    indx = (1:nd)';                 % Indices of coeffs associated to first element and firt monomial
    indx = indx + (0:nz-1)*nd*n;    % Indices of coeffs associated to first element
    indx = reshape(indx,nz*nd,1);
    indx = indx + (0:n-1)*(nz*nd*n + nd);   % Indices of coeffs associated to diagonal elements
    indx = reshape(indx,nz*nd*n,1);
    
    C = reshape(Din.C(indx),nz*nd,n); % Collect coefficients associated to diagonal elements
    C = sum(C,2);                   % Add coefficients corresponding to same dvar and monomial
    C = reshape(C,nd,nz);           % Reshape to the appropriate size
    
    Dout = compress(dpvar(C,Din.degmat,Din.varname,Din.dvarname,[1,1]));
end
end