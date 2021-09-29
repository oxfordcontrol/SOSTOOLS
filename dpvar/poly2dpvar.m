function dpvar_out = poly2dpvar(pvar_in, dvarname)
% This function converts an object of polynomial class to dpvar class
% object
% 
% INPUTS:
% pvar_in: polynomial class object
% dvarnames: cell string of names of variables that need to converted to decision
% variables (this is an optional argument). If dvarnames is not passed then
% the any variable names starting with 'coeff_' is assumed to be a dvaranme
% 
% OUTPUTS:
% dpvar_out: dpvar class equivalent of pvar_in with dvarnames as decision
% variables
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
% Initial coding, DJ, MP, SS - 06/13/2021
% SS - 09/26/2021, added special case of polynomials with one or more zero
% dimensions

if ~isa(pvar_in, 'polynomial')
    error('Input object must be of type polynomial');
end

% remove redundant or unused monomials
pvar_in = combine(pvar_in);

% Extract the full set of variable names
varname_in = pvar_in.varname;

% Extract the set of decision variable names
if nargin==1
    % if no dvarname is specified, take "coeff" vars as dvars
    isdvar = contains(pvar_in.varname,'coeff'); % just output position!
    
    % extract dvarnames
    dvarname = varname_in(isdvar);   
    
    %alternative for dvarname extraction
%    dvarname=pvar_in.varname{dvar_idx}
%    dvarname(cellfun(@isempty,dvarname))=[];    
%    dvarname = cellfun(@(x) x{1}, dvarname,'un',0);

    dindx = 1:length(dvarname);
else
    isdvar = ismember(varname_in,dvarname);         % logical indices of vars that are dvars
    dvarname_old = pvar_in.varname(isdvar);         % list of dvars already present in object
    [~,dindx] = ismember(dvarname_old,dvarname);    % location of already present dvars in dvarname
    % Should we have a unique function to assure no dvar is specified
    % twice? Alternatively, use DPcombine to remove redundant
    % variable names at the end
end

% special case of polynomials with zero dimensions, early exit
if size(pvar_in,1)*size(pvar_in,2)==0
    matdim = pvar_in.matdim;
    pvarname = {};
    dvarname = {};
    degmat = [];
    C = zeros(matdim(1)*(length(dvarname)+1),matdim(2)*(size(degmat,1)));
    dpvar_out = dpvar(C,degmat,pvarname,dvarname,matdim);
    return
end

% Extract the set of polynomial variable names
pvarname=varname_in(~isdvar);

% Separate polynomial and decision variable degmats
degmat_in = pvar_in.degmat;
dvar_degmat = degmat_in(:,isdvar);  
pvar_degmat = degmat_in(:,~isdvar);


% Extract other necessary parameters
[m,n] = size(pvar_in);      
matdim = pvar_in.matdim;    % dimensions of the object do not change
coeff = pvar_in.coeff;      % coefficient array
nd = length(dvarname);      % total number of decision variables
ntt = size(coeff,1);        % number of monomials in the pvar


if ~all((sum(dvar_degmat,2)==0)|(sum(dvar_degmat,2)==1))
    error('Degmat of input polynomial indicates product of decision variables specified in dvarnames. Unable to convert to dpvar class');
end
%toc




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Removed code - alternative implementation to find dpvar degmat
% of the polynomial variables (not decision variables) along with the index
% in which each appears in degmat and varname_in

  %degmat_new=[degmat_in(:,isdvar) degmat_in(:,~isdvar)];
 %degmat_new=[dpvar_degmat pvar_degmat];
% degmat_new=[degmat_in(:,dinds) degmat_in(:,pinds)];
%varname_new=[varname_in(:,dinds) varname_in(:,pinds)]; %not needed - should be equal to
%[dvarname pvarname]
% the pvar structure has not changed, but now the variables are ordered
% first by decision variables, and second by independent variables

%[degmat_new2,term_ind] = sortrows(degmat_new,'descend');
  %[rownum,dvarnum] = find(degmat_new(:,1:nd) ); %vector of decision variable
% coeff_new2=coeff_in(term_ind,:);
% coefficient matrix is now ordered by block according to each decision
% variable. Interesting, but not useful for now

% Create list of index positions in pvar, then create list of corresponding
% index positions in dpvar

  %pvar_degmat=degmat_new(:,(nd+1):(nd+np)); % creates a list of pvar monomials. Will have repeats
% now collapse the dvar part of the new degmat to simply the decision
% variable number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Build a new degmat with contribution of dvars described in first column

[rownum,dvarnum] = find(dvar_degmat);   % dvarnum specifies which degvar occurs in row rownum of degmat (recall that dvars are not multiplied or raised to a power)
dvarnums = spalloc(ntt,1,ntt);
dvarnums(rownum,1) = dindx(dvarnum);    % this vector specifies which dvar (in dvarnames) occurs in each term in degmat_new 
degmat_new2 = [dvarnums pvar_degmat];   % a new degmat with the first column describing which dvar appears in each term


% % Create a list of all possible monomials using the list of all
% possible pvar monomials. Then put in the collapsed pvar format
% [pvar_degmat_unique] = unique(pvar_degmat,'rows'); % set of unique pvar monomials
[pvar_degmat_unique] = unique(pvar_degmat,'rows','stable'); % use stable to retain order
ntu = size(pvar_degmat_unique,1);   % number of unique monomials (terms in the output degmat)


% % Create expanded degmat, of (nd+1) blocks, each consisting on ntu rows and npvar+1 columns
% each row corresponds to one term possibly occuring in the pvar_in object
% the first column describes which dvar the term is multiplied with, using
% dvar0 to account for multiplication with constant factor 1
% the remaining columns describe the power of each polynomial as they occur
% in the considered term
degmat_full_part = [kron((1:nd)',ones(ntu,1)) repmat(pvar_degmat_unique,nd,1)]; % degmat_unique repeated for each possible mutiplcation with a dvar
degmat_full = [zeros(ntu,1) pvar_degmat_unique; degmat_full_part];              % expanded version to include multiplication of degmat_unique with constant factor "1"


% % Now, extract rows in which each of the terms actually occuring in 
% pvar_in appears in our "degmat_full" object (the full list of possible terms)
[~,LOC_DF] = ismember(degmat_new2,degmat_full,'rows');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Removed code
% Now create a larger pvar with lots of zeros
%
% [~,LOC_DF] = ismember(degmat_new,degmat_full,'rows');
%
% This is not much faster than ismember:
% [degmat_new2,IA,IB] = intersect(degmat_full,degmat_new,'rows','stable');
% degmat_new2 = degmat_full(IA,:) % need the positions for each row in
% degmat_new as they appear in degmat full, 
% coeff_new(IA,:)=coeff(IB,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Create a sparse coefficient array "coeff_new" of (nd+1) blocks, each consisting on ntu rows
% each block corresponds to multiplication with one of the (nd+1) dvars: "1",dvarname{1},dvarname{2},...,dvarname{end}
% each row in each block corresponds to multiplication with one of the ntu unique monomial ters
% each column corresponds to one element of the polynomial matrix object
coeff_new = spalloc((nd+1)*ntu,n*m,nnz(coeff));
coeff_new(LOC_DF,:) = coeff;

% % Reshape "coeff_new" to an appropriate coefficient matrix C
% Smaller implementation of the reshaping
% first split into cells describing each column of the polynomial matrix
%B3 = mat2cell(coeff_new,ntu*(nd+1), ones(1,m*n));
B3 = mat2cell(coeff_new,ntu*(nd+1), m*ones(1,n));

% % Now reshape and then transpose each cell (which is a vector)
% [b1]                      [b1']
% [b2]   -> [b1 b2 b3] ->   [b2']
% [b3]                      [b3']
%B3r = cellfun(@(x) reshape(x,[ntu,nd+1]), B3, 'un', 0);
B3r = cellfun(@(x) reshape(x,[ntu,m*(nd+1)]), B3, 'un', 0);
%B3t = cellfun(@transpose, B3r,'un',0);
B3t = cellfun(@transpose, B3r,'un',0);

% % Finally, remerge cells to obtain coeff matrix C
%C = cell2mat(reshape(B3t,m,n));
C = cell2mat(reshape(B3t,1,n));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removed code - alternative implementation to find C
% This part still takes 100 times longer than the rest
% tic
% % the pvar is now of the format coeff_new, degmat_full. Now forming the dpvar should be not too hard. 
% % This is variation of the Bshape command
% 
% % first convert to cell block matrix vector of size ntu x matdim(1)*matdim(2)
% B = mat2cell(coeff_new,(ntu)*ones(1,(nd+1)), ones(1,m*n));
% 
% 
% % c2=reshape(coef_new,ntu*(nd+1)*m*n,1);
% % c3=c2';
% 
% % transpose each element of the block matrix and convert back to a matrix
% B2cell = cellfun(@transpose, B, 'un', 0);
% B2 = cell2mat(B2cell);
% % now reform a new cell block matrix of size 1 x m*n
% B3 = mat2cell(B2,(nd+1), ntu*ones(1,m*n));
% % reshape 1 x mn vector to m x n matrix
% C= cell2mat(reshape(B3,m,n));
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dpvar_out = dpvar(C,pvar_degmat_unique,pvarname,dvarname,matdim);

end
