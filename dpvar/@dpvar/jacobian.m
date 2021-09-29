function J = jacobian(obj,vars)
% J = jacobian(obj,vars) computes the Jacobian matrix for a dpvar object
% in the independent variables vars
% 
% INPUTS:
% obj: a dpvar class object
% vars: a cell or polynomial object specifying the varnames bases on which 
% to compute the Jacobian. If not specified, it will default to obj.varname
% 
% OUTPUTS:
% J: dpvar object corresponding to the Jacobian matrix of obj
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
% Initial coding, DJ, MP, SS - 09/27/2021
% reusing some code from multipoly toolbox polynomial jacobian - start

if nargin == 1
    vars = obj.varname;
end

% Convert x to a cell array of strings
if ispvar(vars)
    vars = char(vars);
elseif ischar(vars)
    vars = cellstr(vars);
elseif ~iscellstr(vars)
    error('vars must be an array of polynomial variables or strings');
end
vars = vars(:);
Nx = size(vars,1);

% Get polynomial info about F
coef = obj.C;
dmat = obj.degmat;
vname = obj.varname;
dvname = obj.dvarname;
mdim = obj.matdim;
nr = mdim(1);
nc = mdim(2);
nz = size(dmat,1);
np = length(vname);


% First build the block matrix
Jdeg = spalloc(nz*Nx,np,Nx*nnz(dmat));
coef = mat2cell(coef, (length(dvname)+1)*ones(nr,1), size(dmat,1)*ones(nc,1));
Cnew = cell(nr,nc);
for i=1:Nx
    xi = vars(i);
    
    % Find variable we are differentiating with respect to
    varnumb = find(strcmp(vname,xi));
    
    % Differentiate
    if isempty(varnumb)
        Jcoef = zeros(1,size(dmat,1));
    else
        Jcoef = dmat(:,varnumb)';
        
        tmpdeg = dmat;
        tmpdeg(:,varnumb) = max( dmat(:,varnumb)-1 , 0);
        Jdeg((i-1)*nz+1:i*nz,:) = tmpdeg;
    end
    zleft = sparse([],[],[],(length(dvname)+1),size(dmat,1)*(i-1));
    zright = sparse([],[],[],(length(dvname)+1),size(dmat,1)*(Nx-i));
    Cnew = cellfun(@(x,y) [x,[zleft,y.*Jcoef,zright]], Cnew, coef,'un',0);
    
end
% Combine blocks to build coefficient matrix
Cnew = cell2mat(Cnew);

% Build the dpvar Jacobian object, and get rid of non-unique monomial terms
J = dpvar(Cnew,Jdeg,vname,dvname,[nr nc*Nx]);
J = combine(J);
end