function Eint = int(E,var,L,U)
% Eint = int(E,var,L,U) integrates the dpvar E in independent variable var
% from lower limit L to upper limit U
% 
% INPUTS:
% E: A dpvar 
% var: A scalar pvar variable (1x1 polynomial object)
% L: Either a scalar pvar variable or scalar value
%       - If a pvar variable, it cannot match input "var"
%       - This input provides the lower limit of the definite integral
% U: Either a scalar pvar variable or scalar value
%       - If a pvar variable, it cannot match input "var"
%       - This input provides the upper limit of the definite integral
% 
% OUTPUTS:
% Eint: integrated dpvar object: Eint = int_L^U (E) d var
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, MP, SS - 07/12/2021

% correction for zero dimension dpvar case
if numel(E.C)==0
    Eint = dpvar(zeros(size(E)));
    return
end

% Integrate the monomial terms specified by E.degmat
nz = size(E.degmat,1);
Z = polynomial(speye(nz),E.degmat,E.varname,[nz,1]);
temp = int_p(Z,var,L,U);

% Adjust the coefficients to account for the integration
% multiply each row/col block of dpvar.C with the integration matrix temp
Cn = E.C*kron(speye(E.matdim(2)),temp.coeff');

% Set numerically insignificant coefficients to zero
[i,j,vals] = find(Cn);
idx = (abs(vals)>=1e-15);
Cn = sparse(i(idx),j(idx),vals(idx),size(Cn,1),size(Cn,2));

Eint = dpvar(Cn, temp.degmat, temp.varname, E.dvarname, E.matdim);

end