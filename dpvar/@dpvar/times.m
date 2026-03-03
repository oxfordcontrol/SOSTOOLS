function C = times(E,F)
% C = times(E,F) takes the elementwise product of a dvpar object and a pvar 
% polynomial object.
%
% INPUTS:
% E: A dpvar or pvar variable (only one object can be a dpvar), dimension
%       [m n]
% F: A dpvar or pvar variable (only one object can be a dpvar), dimension 
%       [m n]
%
% OUTPUTS:
% C: a dpvar object, dimension [m n]
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  M. Peet, S. Shivakumar, D. Jagt
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
% DJ, 03/03/2026: Initial coding;

% If either object is double, convert to polynomial
if ~isa(E,'dpvar') && ~isa(E,'polynomial')
    E = polynomial(E);
end
if ~isa(F,'dpvar') && ~isa(F,'polynomial')
    F = polynomial(F);
end

% Start multiplication process, find which object is dpvar
if isa(E,'dpvar') && isa(F,'dpvar') % 2 dpvars should not be multiplied
    error('Cannot multiply two dpvars');

elseif all(size(E)==1) || all(size(F)==1)
    % Use mtimes function in scalar case
    C = E*F;
    return
elseif isa(F,'dpvar') % right object is dpvar --> multiply transposes instead
    C_transpose = F.'.*E.';
    if all(size(C_transpose.C)==0) % correction for zero dimension dpvars
        C_transpose.C = sparse(size(C_transpose,1)*(length(C_transpose.dvarname)+1), size(C_transpose,2)*size(C_transpose.degmat,1));
    end
    C = C_transpose.';
    return
end
    
% At this point, E must be the dpvar object left object is dpvar
Edegmat=E.degmat;       Fdegmat=F.degmat; % we assume there are no repeats in degmat...
Epvarname = E.varname;  Fpvarname = F.varname;
ntE=size(Edegmat,1);    ntF=size(Fdegmat,1);

% synchronize the pvarnames
[pvarnames_all,IE_p,~] = union(Epvarname,Fpvarname,'stable');   % returns indices in E
[~,~,IF_p] = intersect(Fpvarname,pvarnames_all,'stable');       % returns indices in F
np_all=length(pvarnames_all);

% adjust the degmats accordingly
Edegmat_new = spalloc(ntE,np_all,nnz(Edegmat));
Fdegmat_new = spalloc(ntF,np_all,nnz(Fdegmat));
Edegmat_new(:,IE_p) = Edegmat;
Fdegmat_new(:,IF_p) = Fdegmat;

% adjust E and F to use the new varnames and degmats
F = polynomial(F.coef,Fdegmat_new,pvarnames_all,F.matdim);
E.degmat = Edegmat_new;
E.varname = vec(pvarnames_all); % NOTE: if E and F have only 1 varname, need vec to get appropriate dim

% Make sure the dimensions of E and F match
if size(E,1)~=size(F,1)
    if size(E,1)==1
        E.C = repmat(E.C,size(F,1),1);
        E.matdim(1) = size(F,1);
    elseif size(F,1)==1
        mdim = [size(E,1),size(F,2)];
        CF_new = repelem(F.coefficient,1,size(E,1));
        F = polynomial(CF_new,F.degmat,F.varname,mdim);
    else
        error("Row dimensions of objects must do not match.")
    end
end
if size(E,2)~=size(F,2)
    if size(E,2)==1
        E.matdim(2) = size(F,2);
        E.C = repmat(E.C,1,size(F,2));
    elseif size(F,2)==1
        mdim = [size(F,1),size(E,2)];
        CF_new = repmat(F.coefficient,1,size(E,2));
        F = polynomial(CF_new,F.degmat,F.varname,mdim);
    else
        error("Column dimensions of objects do not match.")
    end
end

% Finally, we perform the multiplication
[m,n] = size(E);
Cmat_E = E.C;
ndecvars = numel(E.dvarname)+1;
nZ = size(E.C,2)/E.matdim(2);

% Establish a unique list of monomials appearing in E.*F
degmat_full = repelem(E.degmat,size(F.degmat,1),1) + repmat(F.degmat,size(E.degmat,1),1);
[Pmat, degmat] = uniquerows_integerTable(degmat_full);
nZ_new = size(degmat,1);

% Loop over all elements of the matrix-valued object
ridcs_new = zeros(0,1);
cidcs_new = zeros(0,1);
vals_new = zeros(0,1);
for k=1:numel(E)
    % Establish row and column number in E
    [rnum,cnum] = ind2sub(size(E),k);
    ridcs_E = (rnum-1)*ndecvars+1:rnum*ndecvars;
    cidcs_E = (cnum-1)*nZ+1:cnum*nZ;
    % Take the product of the kth element of E and F
    Ck = kron(Cmat_E(ridcs_E,cidcs_E),F.coefficient(:,k)');
    Ck = Ck*Pmat;
    % Add the coefficients representing this element to the list
    [ridcs_k,cidcs_k,vals_k] = find(Ck);
    ridcs_new = [ridcs_new; ridcs_k+(rnum-1)*ndecvars];
    cidcs_new = [cidcs_new; cidcs_k+(cnum-1)*nZ_new];
    vals_new = [vals_new; vals_k];
end
% Declare the new coefficient matrix
Cnew = sparse(ridcs_new,cidcs_new,vals_new,m*ndecvars,n*nZ_new);

% Declare the elementwise product
C = E;
C.C = Cnew;
C.degmat = degmat;

end