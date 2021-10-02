function Esub = subs(E, pvars, val)
% Esub = subs(E, pvars, val) substitutes the pvar values by val and returns
% the resulting dpvar Esub
% 
% INPUTS:
% E: a dpvar 
% pvars: polynomial or cell indicating the pvar variable(s) to substitute
% val: double, polynomial, or cell of same size as pvars, indicating the
% values or new variable name assigned to variables pvars
% 
% OUTPUTS:
% Esub: a dpvar object, Esub = E(pvars=val)
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
% Initial coding DJ, MP, SS - 07/19/2021
% 10/01/2021, DJ: Small correction in case where substitution returns a 
% double.

% Account for some trivial cases
if numel(E.C)==0
    Esub = dpvar(zeros(size(E)));
    return
end
if isempty(E.varname)
    Esub = E;
    return
end

% In case pvars is a cell, construct a corresponding polynomial array
if isa(pvars,'cell')
    tmp = [];
    for indx=1:length(pvars)
        if isa(pvars{indx},'char')
            tmp = [tmp; polynomial(pvars(indx))];
        else
            tmp = [tmp; polynomial(pvars{indx})];
        end
    end
    pvars = tmp;
end
% In case val is a cell, construct a corresponding polynomial array
% NOTE that this only works if val is a cell of varnames, not actual values
if isa(val,'cell')
    tmp = [];
    for indx=1:length(val)
        if isa(val{indx},'char')
            tmp = [tmp; polynomial(val(indx))];
        else
            tmp = [tmp; polynomial(val{indx})];
        end
    end
    val = tmp;
end

% Perform polynomial substitution on the monomial terms specified by E.degmat
nz=size(E.degmat,1);
Z=polynomial(speye(nz),E.degmat,E.varname,[nz,1]);
temp=subs(Z,pvars,val);

% Build the new dpvar based on these substituted monomial terms
if isa(temp,'double') % DJ 10/01/21
    Esub = dpvar(E.C*kron(speye(E.matdim(2)),temp), zeros(1,0), {}, E.dvarname, E.matdim);
else % multiply only coefficients if temp is a polynomial 
    Esub = dpvar(E.C*kron(speye(E.matdim(2)),temp.coeff'), temp.degmat, temp.varname, E.dvarname, E.matdim);
end
end
