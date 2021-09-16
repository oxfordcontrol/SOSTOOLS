function outdpvar = combine(indpvar,opts)
% outdpvar = DPcombine(indpvar,opts) performs combine on degmat of dpvar 
% object
% 
% INPUTS:
% indpvar: dpvar class object
% opts: optional char input
%       - If opts='extended', an extended version of DPcombine will be
%         used, removing non-unique varnames and dvarnames of indpvar, in
%         addition to the non-unique rows of indpvar.degmat.
% 
% OUTPUTS:
% outdpvar: A dpvar object with degmat that has unique rows (in other words
%           unique monomials in pvars only). Note that outdpvar.varname and 
%           outdpvar.dvarnames may still contain repeats, unless
%           opts='extended' is used.
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
% Initial coding DJ,MP, SS  - 06/27/2021

if ~isa(indpvar,'dpvar')
    error("DPcombine() needs dpvar variables as inputs");
end

% When extended DPcombine is requested, use alternative function, also
% combining non-unique varnames and dvarnames
if nargin==2 && strcmpi(opts,'extended')
    outdpvar = DPcombine_extended(indpvar);
    return
end

% First extract monomials and find unique monomials
dmat = indpvar.degmat;
[d, ~, idxd] = unique(dmat,'rows');
A = sparse(1:length(idxd),idxd,ones(length(idxd),1),size(dmat,1), size(d,1));

% Add columns corresponding to repeated monomials
C = indpvar.C;
Cnew = C*kron(eye(size(indpvar,2)),A);
outdpvar = dpvar(Cnew, d, indpvar.varname,indpvar.dvarname,indpvar.matdim);
end



function outdpvar = DPcombine_extended(indpvar)
% outdpvar = DPcombine_extended(indpvar) combines identical varnames, 
% dvarnamaes, and rows in the degmat of a dpvar object indpvar
% 
% INPUTS:
% indpvar: dpvar class object
% 
% OUTPUTS:
% outdpvar: A dpvar object with unique rows in outdpvar.degmat, and no
%           repeats in either outdpvar.dvarname or outdpvar.varname.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ,MP,SS  - 06/27/2021

% Extract necessary values
dvarname = indpvar.dvarname;
varname = indpvar.varname;
dmat = indpvar.degmat;
C = indpvar.C;

% Combine dvars that appear multiple times
nd = length(dvarname);
if nd>1
    [newC,newdvar] = DPVuniquedvar(C,dvarname);
else
    newC = C;
    newdvar = dvarname;
end

% Combine pvars that appear multiple times
np = length(varname);
if np>1
    [newdeg,newpvar] = DPVuniquepvar(dmat,varname);
else
    newdeg = dmat;
    newpvar = varname;
end

% Combine monomial terms that appear multiple times.
nt = size(dmat,1);
if nt>1
    [newnewC,newnewdeg] = DPVuniqueterm(newC,newdeg);
else
    newnewC = newC;
    newnewdeg = newdeg;
end

% Build the combined object
outdpvar = dpvar(newnewC,newnewdeg,newpvar,newdvar,indpvar.matdim);

end