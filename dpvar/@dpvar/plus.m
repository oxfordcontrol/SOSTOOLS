function G = plus(E,F)
% G = plus(E,F) performs the addition operation E+F 
% 
% INPUTS:
% E,F: dpvar or polynomials of same dimension (at least one is a dpvar
%   object)
% 
% OUTPUTS:
% G: a dpvar object, G = E+F
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
% Initial coding DJ,MP,SS - 07/07/2021

% For scalar addition, add the scalar to all terms in the dpvar object
if all(size(E)==1)
    E = E*ones(size(F));
elseif all(size(F)==1)
    F = F*ones(size(E));
end

% Otherwise, dimensions must match
if ~all(size(E)==size(F))
    error('Objects being added do not have the same dimension');
end


% Treat case where one object is a matrix separately
if isa(E,'double')
    % Check if a constant term is already present in F
    [r_indx,~] = find(F.degmat);
    z_indx = setdiff((1:size(F.degmat,1))',unique(r_indx));
    if isempty(z_indx)
        % We have to add a zero row in F.degmat, and add E to F.C accordingly
        Gdegmat = [F.degmat;zeros(1,size(F.degmat,2))];
        nt = size(Gdegmat,1);   nd = length(F.dvarname)+1;
        [Ci,Cj,Cval] = find(F.C);
        Cj = Cj + floor((Cj-1)/(nt-1));     % take into account the additional row in the degmat
        newi = repmat((1: nd : nd*size(F,1))',[size(E,2),1]);       % row indices specifying constant term
        newj = kron((nt : nt : nt*size(F,2))',ones(size(E,1),1));   % column indices specifying constant term
        newC = vec(E);      % constant values to be assigned to row and column indices
        GC = sparse([vec(Ci);newi],[vec(Cj);newj],[vec(Cval);newC],nd*size(F,1),nt*size(F,2));
    else
        % We have to add E to F.C for the already present zero row in F.degmat
        Gdegmat = F.degmat;
        nt = size(Gdegmat,1);   nd = length(F.dvarname)+1;
        [Ci,Cj,Cval] = find(F.C);
        
        % Add the coefficients in E to F.C in the appropriate locations
        newi = repmat((1: nd : nd*size(F,1))',[size(E,2),1]);
        newj = kron((z_indx : nt : nt*size(F,2))',ones(size(E,1),1));
        newC = vec(E);
        GC = sparse([vec(Ci);newi],[vec(Cj);newj],[vec(Cval);newC],nd*size(F,1),nt*size(F,2));
        % Note that matlab will take care of assignments (i,j) appearing
        % twice in sparse by adding the associated coefficients
    end
    G = dpvar(GC,Gdegmat,F.varname,F.dvarname,F.matdim);
    return
    
elseif isa(F,'double')
    % Check if a constant term is already present in E
    [r_indx,~] = find(E.degmat);
    z_indx = setdiff((1:size(E.degmat,1))',unique(r_indx));
    if isempty(z_indx)
        % We have to add a zero row in E.degmat, and add F to E.C accordingly
        Gdegmat = [E.degmat;zeros(1,size(E.degmat,2))];
        nt = size(Gdegmat,1);  nd = length(E.dvarname)+1;
        [Ci,Cj,Cval] = find(E.C);
        Cj = Cj + floor((Cj-1)/(nt-1)); % take into account the additional row in the degmat
        newi = repmat((1: nd : nd*size(E,1))',[size(F,2),1]); % row indices specifying constant term
        newj = kron((nt : nt : nt*size(E,2))',ones(size(F,1),1)); % column indices specifying constant term
        newC = vec(F); % constant values to be assigned to row and column indices
        GC = sparse([vec(Ci);newi],[vec(Cj);newj],[vec(Cval);newC],nd*size(E,1),nt*size(E,2));
    else
        % We have to add F to E.C for the already present zero row in E.degmat
        Gdegmat = E.degmat;
        nt = size(Gdegmat,1);   nd = length(E.dvarname)+1;
        [Ci,Cj,Cval] = find(E.C);
        
        % Add the coefficients in F to E.C in the appropriate locations
        newi = repmat((1: nd : nd*size(E,1))',[size(F,2),1]);
        newj = kron((z_indx : nt : nt*size(E,2))',ones(size(F,1),1));
        newC = vec(F);
        GC = sparse([vec(Ci);newi],[vec(Cj);newj],[vec(Cval);newC],nd*size(E,1),nt*size(E,2));
        % Note that matlab will take care of assignments (i,j) appearing
        % twice in sparse by adding the associated coefficients
    end
    G = dpvar(GC,Gdegmat,E.varname,E.dvarname,E.matdim);
    return
end
    

% Convert polynomials to dpvar with same set of decision variables
E_temp = E; F_temp = F;
if isa(E,'polynomial')
    E_temp = poly2dpvar(E);
elseif isa(F,'polynomial')
    F_temp = poly2dpvar(F);
elseif ~isa(E,'dpvar') || ~isa(F,'dpvar')
    error('Addition is only supported for objects of type double, polynomial, or dpvar');
end

% Check the density of the coefficient matrices
dnstyE = nnz(E_temp.C)/(size(E_temp.C,1)*size(E_temp.C,2));
dnstyF = nnz(F_temp.C)/(size(F_temp.C,1)*size(F_temp.C,2));

% Use sparse implementation depending on density of coefficient matrices 
if dnstyE>=0.5 || dnstyF>=0.5 % (note that this bound is purely heuristic)
    
    % Give E and F the same varname, dvarname, and degmat...
    [E_temp2,F_temp2] = common_basis(E_temp,F_temp);
    % ...then we can add the coefficients...
    C = E_temp2.C + F_temp2.C;
    % ...and then we can build the summed dpvar, getting rid of
    % zero-contributions using DPcompress
    G = compress(dpvar(C,E_temp2.degmat,E_temp2.varname,E_temp2.dvarname,E_temp2.matdim));
    
else
    % In sparse case, use sparse_plus (see below)
    G = compress(DP_sparse_plus(E_temp,F_temp));
end

end



function G = DP_sparse_plus(E,F)
% G = DP_sparse_plus(E,F) performs the addition operation E+F 
% 
% INPUTS:
% E,F: dpvar or polynomials of same dimension (at least one is a dpvar
%   object)
% 
% OUTPUTS:
% G: a dpvar object, G = E+F
% 
% NOTES:
% The function works by introducing a new dpvar Gfull, with
% Gfull.varname = [E.varname;F.varname];
% Gfull.dvarname = [E.dvarname;F.varname];
% Gfull.degmat = blkdiag(E.degmat,F.degmat);
% and with Gfull.C an appropriate combination of E.C and F.C;
% Then, all non-unique varnames, dvarnames, and rows in degmat are combined
% using DPcombine_extended, yielding the sum G = E+F;
%
% The function was found to work best when E.C and F.C are sparse
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ,MP,SS - 07/07/21

nd1 = length(E.dvarname);   nd2 = length(F.dvarname);
nt1 = size(E.degmat,1);     nt2 = size(F.degmat,1);
mdim = size(E,1);           ndim = size(E,2);

% Concatenate the dvars and degmat rows of the two objects
dvarname = [E.dvarname; F.dvarname];
pvarname = [E.varname; F.varname];
degmat = blkdiag(E.degmat,F.degmat);
ndn = nd1 + nd2;
ntn = nt1 + nt2;

% % Adjust the coefficents accordingly:

% Extract row and column indices of nonzeros in E.C, and the associated
% coefficients
[CEi,CEj,CEval] = find(E.C);
CEi = vec(CEi);     CEj = vec(CEj);     CEval = vec(CEval);

% Determine which dvar and row in degmat each i and j index corresponds to
CEi_uniq = mod(CEi-1,nd1+1) + 1;
CEj_uniq = mod(CEj-1,nt1) + 1;

% Shift the row and column indices to account for the coefficients from F 
% that will be added
CEi_shift = floor((CEi-1)/(nd1+1)) * (ndn+1);
CEi_new = CEi_uniq + CEi_shift;

CEj_shift = floor((CEj-1)/(nt1)) * ntn;
CEj_new = CEj_uniq + CEj_shift;

% Extract row and column indices of nonzeros in F.C, and the associated
% coefficients
[CFi,CFj,CFval] = find(F.C);
CFi = vec(CFi);     CFj = vec(CFj);     CFval = vec(CFval);

% Determine which dvar and row in degmat each i and j index corresponds to
CFi_uniq = mod(CFi-1,nd2+1) + 1;
CFj_uniq = mod(CFj-1,nt2) + 1;

% Separate rows corresponding to a constant term (rather than a dvar)
Fi_cnstnt = (CFi_uniq==1);    % rows of F.C corresponding to constant term
CFi1 = CFi(Fi_cnstnt);
CFj1 = CFj(Fi_cnstnt);
CFval1 = CFval(Fi_cnstnt);

% Extract values for the remaining rows
CFi2 = CFi(~Fi_cnstnt) - 1;
CFj2 = CFj(~Fi_cnstnt);
CFval2 = CFval(~Fi_cnstnt);

% Shift the row and column indices to indicate where the coefficients from
% F should be added to the new coefficient matrix
CFi1_shift = floor((CFi1-1)/(nd2+1)) * (ndn+1);
CFi1_new = CFi_uniq(Fi_cnstnt) + CFi1_shift;

CFi2_shift = floor((CFi2-1)/(nd2+1)) * (ndn+1) + nd1;
CFi2_new = CFi_uniq(~Fi_cnstnt) + CFi2_shift;

CFj1_shift = floor((CFj1-1)/(nt2)) * ntn + nt1;
CFj1_new = CFj_uniq(Fi_cnstnt) + CFj1_shift;

CFj2_shift = floor((CFj2-1)/(nt2)) * ntn + nt1;
CFj2_new = CFj_uniq(~Fi_cnstnt) + CFj2_shift;

% Build the new coefficient matrix
Ci = [CEi_new; CFi1_new; CFi2_new];
Cj = [CEj_new; CFj1_new; CFj2_new];
Cval = [CEval; CFval1; CFval2];

C = sparse(Ci,Cj,Cval,(ndn+1)*mdim,ntn*ndim,length(Cval));

% Build the new dpvar, combining non-unique vars, dvars, and monomial terms
G = combine(dpvar(C,degmat,pvarname,dvarname,E.matdim),'extended');

end