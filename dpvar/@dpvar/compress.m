function outdpvar = compress(indpvar)
% outdpvar = DPcompress(indpvar) removes decision variables and monomial 
% terms from D that do not appear in the actual polynomial object
% 
% INPUTS:
% indpvar: a dpvar variable
% 
% OUTPUTS:
% outdpvar: a dpvar object with no zero rows or columns in outdpvar.C
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
% Initial coding DJ, MP, SS - 07/01/2021

if ~isa(indpvar,'dpvar')
    error("DPcompress() needs a dpvar variable as input");
end

% Extract components of the dpvar
C = indpvar.C;
degmat = indpvar.degmat;
pvarname = indpvar.varname;
dvarname = indpvar.dvarname;
mdim = size(indpvar,1); 
ndim = size(indpvar,2);
nd = length(dvarname) + 1;
ntt = size(degmat,1);

% Establish which dvars appear in the actual dpvar
Crows = any(C,2);                       % Logical array indicating which rows contain a nonzero value
dvars = any(reshape(Crows,nd,mdim),2);  % Logical array inidcating which dvars {1;dvar1;...} actually appear in the object indpvar
dvars(1) = 1;                           %  - make sure not to exclude constant term
Crows_new = repmat(dvars,mdim,1);       % Logical array indicating which rows of C should be retained

% Establish which monomial terms appear in the actual dpvar
Ccols = vec(any(C,1));                      % Logical array indicating which columns contain a nonzero value
degrows = any(reshape(Ccols,ntt,ndim),2);   % Logical array inidcating which monomial terms {1;...;ntt} actually appear in the object indpvar
Ccols_new = repmat(degrows,ndim,1);         % Logical array indicating which columns of C should be retained

% Compress C and degmat to get rid of dvars and monomials that do
% not contribute
C_new = C(Crows_new,Ccols_new);
if isempty(C_new)
    outdpvar = dpvar(zeros(mdim,ndim),zeros(1,0),{},{},[mdim,ndim]);
    return
end
degmat_new = degmat(degrows,:);

% Establish which pvars appear in the actual polynomial
degcols = any(degmat_new,1);

% Retain only dvars and pvars that contribute to the polynomial
dvarname_new = dvarname(dvars(2:end));   
pvarname_new = pvarname(degcols);

% Compress degmat to get rid of pvars that do not contribute
degmat_new = degmat_new(:,degcols);

% Build the compressed dpvar
outdpvar = dpvar(C_new,degmat_new,pvarname_new,dvarname_new,[mdim,ndim]);

end