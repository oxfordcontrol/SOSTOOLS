function [Enew, Fnew] = common_basis(E,F)
% [Enew,Fnew] = DPcommon_basis(E,F) extends the coeffs of E and F to have
% same left and right bases
% 
% INPUTS:
% E,F : dpvar class objects
% 
% OUTPUTS:
% Enew, Fnew: dpvar objects, that have same dvarnames, varnames, and degmats
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
% Initial coding DJ, MP, SS - 07/07/2021


% Check the density of the coefficient matrices
dnstyE = nnz(E.C)/(size(E.C,1)*size(E.C,2));
dnstyF = nnz(F.C)/(size(F.C,1)*size(F.C,2));

% If the coefficient matrices are sparse, use the sparse implementation
if dnstyE<=0.2 && dnstyF<=0.2
    [Enew, Fnew] = DPcommon_basis_sparse(E,F); % see end of this file
    return
end

% extract standard information
m = E.matdim(1);            n = E.matdim(2);
p = F.matdim(1);            q = F.matdim(2); 
Edegmat = E.degmat;         Fdegmat = F.degmat; % we assume there are no repeats in degmat...
Edvarname = E.dvarname;     Fdvarname = F.dvarname;
Epvarname = E.varname;      Fpvarname = F.varname;
Ecoeff = E.C;               Fcoeff = F.C;
ndE = length(Edvarname);    ndF = length(Fdvarname);
ntE = size(Edegmat,1);      ntF = size(Fdegmat,1);

% Check if any decision var appears as independent var in one of the objects
Ep2d = ismember(Epvarname,Fdvarname);  % pvars in E that appear as dvars in F
Fp2d = ismember(Fpvarname,Edvarname);  % pvars in F that appear as dvars in E

if any(Ep2d) || any(Fp2d)
    error('Some decision variables in one dpvar appear as polynomial variables in the other.')
end

% Combine the independent variable names and expand the degmat columns accordingly
[pvarnames_all,IE_p,~] = union(Epvarname,Fpvarname,'stable');
[~,~,IF_p] = intersect(Fpvarname,pvarnames_all,'stable');
np_all=length(pvarnames_all);

Edegmat_new = spalloc(ntE,np_all,nnz(Edegmat));
Fdegmat_new = spalloc(ntF,np_all,nnz(Fdegmat));
Edegmat_new(:,IE_p) = Edegmat;
Fdegmat_new(:,IF_p) = Fdegmat;

% Combine the degmats (rows), and store mapping info:
% Cdegmat: combined degmat
% IE_m_2: indices such that Ecoeff_newd(IE_d_2,IE_m_2) = Ecoeff
% IF_m_2: indices such that Fcoeff_newd(IF_d_2,IF_m_2) = Fcoeff
[Cdegmat,IE_m,~] = union(Edegmat_new,Fdegmat_new,'rows','stable');
nt_all = size(Cdegmat,1);
[~,IF_m] = ismember(Fdegmat_new,Cdegmat,'rows'); 
IE_m_2 = repmat(IE_m,n,1)+nt_all*kron((0:n-1)',ones(ntE,1));
IF_m_2 = repmat(IF_m,q,1)+nt_all*kron((0:q-1)',ones(ntF,1));

% Combine the dvarnames, and store the mapping info:
% dvarnames_all: combined dvarnames
% IE_d_2: indices such that Ecoeff_newd(IE_d_2,IE_m_2) = Ecoeff
% IF_d_2: indices such that Fcoeff_newd(IF_d_2,IF_m_2) = Fcoeff
[dvarnames_all,IE_d,~] = union(Edvarname,Fdvarname,'stable');
[~,~,IF_d] = intersect(Fdvarname,dvarnames_all,'stable');
nd_all = length(dvarnames_all)+1;
IE_d_2 = repmat([1;IE_d+1],m,1)+nd_all*kron((0:m-1)',ones(ndE+1,1)); 
IF_d_2 = repmat([1;IF_d+1],p,1)+nd_all*kron((0:p-1)',ones(ndF+1,1));

% Build the new coefficient matrices for new dvarnames list and new degmat
% by expanding the row/cols and rearranging
IE_d_3 = repmat(IE_d_2,[length(IE_m_2),1]);
IE_m_3 = kron(IE_m_2,ones(length(IE_d_2),1));
IF_d_3 = repmat(IF_d_2,[length(IF_m_2),1]);
IF_m_3 = kron(IF_m_2,ones(length(IF_d_2),1));
Ecoeff_newd = sparse(IE_d_3,IE_m_3,Ecoeff,(nd_all)*m,nt_all*n);
Fcoeff_newd = sparse(IF_d_3,IF_m_3,Fcoeff,(nd_all)*p,nt_all*q);

% Build the new dpvars
Enew = dpvar(Ecoeff_newd,Cdegmat,pvarnames_all,dvarnames_all,E.matdim);
Fnew = dpvar(Fcoeff_newd,Cdegmat,pvarnames_all,dvarnames_all,F.matdim);

end



function [Enew, Fnew] = DPcommon_basis_sparse(E,F)
% [Enew,Fnew] = DPcommon_basis_sparse(E,F) extends the coeffs of E and F to
% have same left and right bases
% Alternative to the original implementation, found to work better when
% E.C and F.C are (very) sparse
%
% INPUTS:
% E,F : dpvar class objects of any dimensions
% 
% OUTPUTS:
% Enew, Fnew: dpvar objects, that have same dvarnames, varnames, and degmats
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, MP, SS  - 07/07/2021


% % % Extract the necessary parameters % % %

mdimE = E.matdim(1);        ndimE = E.matdim(2);
mdimF = F.matdim(1);        ndimF = F.matdim(2); 
Edegmat = E.degmat;         Fdegmat = F.degmat;
Edvarname = E.dvarname;     Fdvarname = F.dvarname;
Epvarname = E.varname;      Fpvarname = F.varname;

npE = length(Epvarname);    npF = length(Fpvarname);
ndE = length(Edvarname)+1;  ndF = length(Fdvarname)+1;
ntE = size(Edegmat,1);      ntF = size(Fdegmat,1);

Ep2d = ismember(Epvarname,Fdvarname);  % pvars in E that appear as dvars in F
Fp2d = ismember(Fpvarname,Edvarname);  % pvars in F that appear as dvars in E

if any(Ep2d) || any(Fp2d)
    error('Some decision variables in one dpvar appear as polynomial variables in the other.')
end

pvarname_all = [Epvarname;Fpvarname];
dvarname_all = [Edvarname;Fdvarname];
degmat_all = blkdiag(Edegmat,Fdegmat);


% % % Decompose the coefficient matrices % % %

% Extract row and column indices and associated values of nonzero elements of E.C
[CEi,CEj,CEval] = find(E.C);
CEi = vec(CEi);     CEj = vec(CEj);     CEval = vec(CEval);

% Establish the dvar, monomial term, row, and column to which each element
% contributes
Erow = floor((CEi-1)/ndE) + 1;   % row number
Ecol = floor((CEj-1)/ntE) + 1;   % col number
Edvar = CEi - (Erow-1)*ndE;      % dvar number
Emon = CEj - (Ecol-1)*ntE;       % monomial number

% Separate the elements corresponding to a constant term (dvar number = 1)
Ei_cnstnt = (Edvar==1);
CEval1 = CEval(Ei_cnstnt);      CEval2 = CEval(~Ei_cnstnt);
Edvar1 = Edvar(Ei_cnstnt);      Edvar2 = Edvar(~Ei_cnstnt)-1;
Emon1 = Emon(Ei_cnstnt);        Emon2 = Emon(~Ei_cnstnt);
Erow1 = Erow(Ei_cnstnt);        Erow2 = Erow(~Ei_cnstnt);
Ecol1 = Ecol(Ei_cnstnt);        Ecol2 = Ecol(~Ei_cnstnt);

% repeat above steps for F dpvar
% Extract row and column indices and associated values of nonzero elements of F.C
[CFi,CFj,CFval] = find(F.C);
CFi = vec(CFi);     CFj = vec(CFj);     CFval = vec(CFval);

% Establish the dvar, monomial term, row, and column to which each element
% contributes
Frow = floor((CFi-1)/ndF) + 1;
Fcol = floor((CFj-1)/ntF) + 1;
Fdvar = CFi - (Frow-1)*ndF;
Fmon = CFj - (Fcol-1)*ntF;

% Separate the elements corresponding to a constant term (rather than a dvar)
Fi_cnstnt = (Fdvar==1);
CFval1 = CFval(Fi_cnstnt);      CFval2 = CFval(~Fi_cnstnt);
Fdvar1 = Fdvar(Fi_cnstnt);      Fdvar2 = Fdvar(~Fi_cnstnt)-1;
Fmon1 = Fmon(Fi_cnstnt);        Fmon2 = Fmon(~Fi_cnstnt);
Frow1 = Frow(Fi_cnstnt);        Frow2 = Frow(~Fi_cnstnt);
Fcol1 = Fcol(Fi_cnstnt);        Fcol2 = Fcol(~Fi_cnstnt);

% Establish which dvars and monomial terms from the combined object each
% element of F.C contributes to
Fdvar2 = Fdvar2 + ndE - 1;  % dvar from dvarname_all associated with the elements of F.C
Fmon1 = Fmon1 + ntE;        % row from degmat_all associated with the elements of F.C
Fmon2 = Fmon2 + ntE;


% % % Get rid of non-unique pvars % % %

if npE>=1 && npF>=1
    [newdegmat_all,pvarname_unique] = DPVuniquepvar(degmat_all,pvarname_all);
else
    newdegmat_all = degmat_all;
    pvarname_unique = pvarname_all;
end

% % % Get rid of non-unique dvars % % %

if ndE-1>=1 && ndF-1>=1
    % Sort the dvarnames, and establish which are repeated
    [dvarsort,dsortidx] = sortrows_integerTable(char(dvarname_all));    % dvarname_all(dsortidx) = dvarsort
    repeatidx = all(dvarsort(2:end,:) == dvarsort(1:end-1,:),2);
    drvec = cumsum([1; ~repeatidx]);

    % Get rid of repeated dvars
    duniqueidx = dsortidx;
    duniqueidx(repeatidx) = [];
    dvarname_unique = dvarname_all(duniqueidx);

    % Establish mapping from dvarname_all to dvarname_unique
    dsort_array = [dsortidx,drvec];
    new_darray = sortrows_integerTable(dsort_array);
    dindices = new_darray(:,2);         % dvarname_all = dvarname_unique(dindices)

    % Adjust the dvar indices
    Edvar2_new = dindices(Edvar2) + 1;  % dvar from [1;unqiuedvar] associated with the elements of E.C
    Fdvar2_new = dindices(Fdvar2) + 1;  % dvar from [1;unqiuedvar] associated with the elements of F.C
else
    dvarname_unique = dvarname_all;
    Edvar2_new = Edvar2 + 1;
    Fdvar2_new = Fdvar2 + 1;
end
    

% % % Get rid of non-unique monomial terms % % %

if ntE>=1 && ntF>=1
    % Sort the rows of degmat according to the total degree of the monomials
    degmat_temp = [sum(newdegmat_all,2) newdegmat_all];
    [degsort,degsortidx] = sortrows_integerTable(degmat_temp);
    degsort = degsort(:,2:end);

    % Establish which rows are not unique
    repeatidx = ~any(degsort(2:end,:) ~= degsort(1:end-1,:),2);
    degrvec = cumsum([1; ~repeatidx]);

    % Create a unique set of monomials
    uniqueidx = degsortidx;
    uniqueidx(repeatidx) = [];
    degmat_unique = newdegmat_all(uniqueidx,:);

    % Establish mapping from rows in newdegmat_all to degmat_unique
    degsort_array = [degsortidx,degrvec];
    new_degarray = sortrows_integerTable(degsort_array);
    degindices = new_degarray(:,2);   

    % Adjust the monomial indices
    Emon1_new = degindices(Emon1);  % row from degmat_unique associated with the elements of E.C
    Fmon1_new = degindices(Fmon1);  % row from degmat_unique associated with the elements of F.C
    Emon2_new = degindices(Emon2);
    Fmon2_new = degindices(Fmon2);
else
    degmat_unique = newdegmat_all;
    Emon1_new = Emon1;
    Fmon1_new = Fmon1;
    Emon2_new = Emon2;
    Fmon2_new = Fmon2;
end


% % % Build the new coefficient matrices, and the new dpvars % % %

% Re-combine the indices
CEval_new = [CEval1; CEval2];
Edvar_new = [Edvar1; Edvar2_new];
Emon_new = [Emon1_new; Emon2_new];
Erow_new = [Erow1; Erow2];
Ecol_new = [Ecol1; Ecol2];

CFval_new = [CFval1; CFval2];
Fdvar_new = [Fdvar1; Fdvar2_new];
Fmon_new = [Fmon1_new; Fmon2_new];
Frow_new = [Frow1; Frow2];
Fcol_new = [Fcol1; Fcol2];

% Determine the appropriate row and column in the new coefficient matrices
ndu = length(dvarname_unique) + 1;
ntu = size(degmat_unique,1);

CEi_new = (Erow_new-1)*ndu + Edvar_new;
CEj_new = (Ecol_new-1)*ntu + Emon_new;
CFi_new = (Frow_new-1)*ndu + Fdvar_new;
CFj_new = (Fcol_new-1)*ntu + Fmon_new;

% Build the coefficient matrices and dpvars
EC_new = sparse(CEi_new,CEj_new,CEval_new,ndu*mdimE,ntu*ndimE,length(CEi));
FC_new = sparse(CFi_new,CFj_new,CFval_new,ndu*mdimF,ntu*ndimF,length(CFi));

Enew = dpvar(EC_new,degmat_unique,pvarname_unique,dvarname_unique,E.matdim);
Fnew = dpvar(FC_new,degmat_unique,pvarname_unique,dvarname_unique,F.matdim);

end