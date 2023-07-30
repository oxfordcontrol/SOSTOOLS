function [sos,P,Q] = sosquadvar(sos,Z1c,Z2c,mdim_in,ndim_in,option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [sos,P] = sosquadvar(sos,Z1c,Z2c,mdim_in,ndim_in,option)
%
% This program declares a symbolic positive scalar semidefinite matrix P of
% size nxn which is positive semidefinite for all values of the variables.
% The matrix has the form
% P = [Z1{1}    ]^T    [Z2{1}     ]
%     [ Z1{1}   ]    Q [ Z2{1}    ]
%     [  Z1{2}  ]      [  Z2{1}   ]
%     [   Z1{2} ]      [   Z2{2}  ], Q>=0 (if option='pos')
%     [    Z1{2}]      [    Z2{2} ]
%                      [     Z2{2}]
%
% Where Z is a vector of monomials. if option='pos', dims is square and Z1
% and Z2 are of equal length, then Q is a positive semidefinite matrix.
%
% NONCELLULAR INPUT FORMAT
% INPUTS:
% prog - The SOS program to which to attach the variable
% Z1c - column vector of monomials (pvar) to be multiplied on the left
% Z2c - column vector of monomials (pvar) to be multiplied on the right
% mdim_in - a scalar indicating the row dimension of the output.
%           If this entry is left empty, it defaults to 1
% ndim_in - a scalar indicating the column dimension of the output.
%           If this entry is left empty, it defaults to 1
% option - 'pos' if Q should be a positive semidefinite matrix
%        - 'sym' if Q should be a symmetric matrix 
%           (note this may not result in P being symmetric, if Z1c~=Z2c)
% In order to use the positive or symmetric options, the length of the
% monomial vectors Z1 and Z2 should be identical. In addition, we require
% mdim_in=ndim_in.
%
% OUTPUTS:
% prog - updated program structure
% P -   dpvar polynomial structure. This structure is not cellular unless
%       Z1c or Z2c is cellular
% Q -   (optional) cellstring corresponding to the decision variables as
%       in the matrix Q. If Z1c or Z2c is cellular, this object will be
%       a cell of cellstr objects as well, with dimensions corresponding to 
%       the sizes of Z1c and Z2c
%
% ALTERNATIVE CELLULAR INPUT FORMAT
% For a multipartite structure, Z1c and Z2c may be cells of monomial column
% vectors. The number of cells in Z1c and Z2c need not be the same unless
% the 'pos' or 'sym' options are used. Note that mdim_in and ndim_in are
% now expected to be vectors indicating the number of rows and columns to
% be output for each of the cells in Z1c and Z2c. 
%
% CELLULAR INPUTS:
% mdim_in - array of same length as number of cells in Z1c, indicating the row dimensions of the
% output P{i,j} - see output description. If mdim_in is a scalar and Z1c
% has multiple cells, the dimension of each cell will be set to mdim_in. 
% ndim_in - array of same length as number of cells in Z2c, indicating the column dimensions of the
% output P{i,j} - see output description. If ndim_in is a scalar and Z2c
% has multiple cells, the dimension of each cell will be set to ndim_in.
%
% CELLULAR OUTPUTS:
% prog - updated program structure
% P - pxq cellular structure of dpvars. p is the number of cells in Z1c and
% q is the number of cells in Z2c. P{i,j} is a matrix dpvar of dimension
% mdin_in(i) x ndin_in(j) of the form 
%    (I_mdim_in(i) otimes Z1c{i})^T Q{i,j} (I_ndim_in(j) otimes Z2c{j})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 4.00.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018, 2021  
%                                      A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5),
%                                      M. Peet (6), D. Jagt (6)
% (1) Department of Engineering Science, University of Oxford, Oxford, U.K.
% (2) Laboratoire de Signaux et Systmes, CentraleSupelec, Gif sur Yvette,
%     91192, France
% (3) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (4) Aerospace and Engineering Mechanics Department, University of
%     Minnesota, Minneapolis, MN 55455-0153, USA.
% (5) Laboratory for Information and Decision Systems, M.I.T.,
%     Massachusetts, MA 02139-4307
% (6) Cybernetic Systems and Controls Laboratory, Arizona State University,
%     Tempe, AZ 85287-6106, USA.
%
% Send bug reports and feedback to: sostools@cds.caltech.edu
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
% along with this program. If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log and developer notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 07/09/21, DJ: split dims into mdim and ndim, specifying number of blocks
% for each monomial Z1c{i} and Z2c{j} respectively. Also separate diagonal
% and off-diagonal computations under toggle==1
% 07/12/21, DJ: Added code to adjust the order of the rows in the degmats
% of Z1c{i} and Z2c{j} to match the rows in the associated monomials
% 07/14/21, DJ: Adjusted specification of dvars to use fastint2str function
% 04/19/22, DJ: Update to output decision variable matrix (if desired),
%               and allow general polynomial basis functions Z1 and Z2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract the inputs, and check for (some) issues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<=2
    error('quadvar requires at least 3 inputs')
elseif nargin==3
    mdim_in = 1;
    ndim_in = 1;
elseif nargin>=4
    if isempty(mdim_in)
        mdim_in = 1;
    end
    if isempty(ndim_in)
        ndim_in = 1;
    end
end
if nargout>3
    error('The function produces at most 3 outputs')
end

% Output will be a cell iff either of the monomials is specified as cell
celltoggle=1;
if ~iscell(Z1c) && ~iscell(Z2c)
    celltoggle=0;
end
% Internally, we will work with cells either way
if ~iscell(Z1c)
    Z1c={Z1c};
end
if ~iscell(Z2c)
    Z2c={Z2c};
end
if isa(Z1c{1},'sym') || isa(Z2c{1},'sym')
    error('quadvar is incompatible with symbolic variables. Please convert your monomials to polynomial/pvar.')
end

% Determine the dimensions of P{i,j} for each combination of Z1c{i}, Z2c{j}
nZ1c=length(Z1c);                   nZ2c=length(Z2c);
mdim_in = reshape(mdim_in,[],1);    ndim_in = reshape(ndim_in,[],1);
if isa(mdim_in,'cell')
    mdim_in = cell2mat(mdim_in);
end
if length(mdim_in)==1
    if nZ1c>=2
        mdim_in=mdim_in*ones(nZ1c,1);
    end
elseif length(mdim_in)~=nZ1c
    error('length of mdim_in must match number of cells in Z1')
end
if isa(ndim_in,'cell')
    ndim_in = cell2mat(ndim_in);
end
if length(ndim_in)==1
    if nZ2c>=2
        ndim_in=ndim_in*ones(nZ2c,1);
    end
elseif length(ndim_in)~=nZ2c
    error('length of ndim_in must match number of cells in Z2')
end

% Determine the number of basis functions in each element of Z1c, Z2c
nZ1=zeros(nZ1c,1);      nZ2=zeros(nZ2c,1);
for i=1:nZ1c
    if isnumeric(Z1c{i}) && Z1c{i} == 1
        pvar Z1;
        Z1c{i} = Z1*0+1;
    end
    nZ1(i)=length(Z1c{i});
    %    mbdim(i)=mdim*nZ1(i);
end
for i=1:nZ2c
    if isnumeric(Z2c{i}) && Z2c{i} == 1
        pvar Z2;
        Z2c{i} = Z2*0+1;
    end
    nZ2(i)=length(Z2c{i});
    %    nbdim(i)=ndim*nZ2(i);
end

% Do we need Q to be symmetric and positive?
postoggle=0;
symtoggle=0;
if nargin>=6
    if strcmp(option,'pos')
        if ~all(mdim_in==ndim_in)
            error('desired sizes must be square when using the positive option')
        elseif nZ1c~=nZ2c
            error('number of monomial vectors must be the same when using the positive option')
        elseif ~all(nZ1==nZ2)
            error('left and right monomial vectors must be same length when using the positive option')
        else
            postoggle=1;
            symtoggle=1;
        end
    elseif strcmp(option,'sym')
        if ~all(mdim_in==ndim_in)
            error('desired size must be square when using the symmetric option')
        elseif nZ1c~=nZ2c
            error('number of monomial vectors must be the same when using the symmetric option')
        elseif ~all(nZ1==nZ2)
            error('left and right monomial vectors must be same length when using the symmetric option')
        else
            symtoggle=1;
        end
    else
        error('option type not reocgnized')
    end
end

% Determine the dimensions of the components Q{i,j}
mbdim = mdim_in.*nZ1;
nbdim = ndim_in.*nZ2;

% Compute the total number of decision variables (in Q) to add
mbdimT = sum(mbdim);
nbdimT = sum(nbdim);
ndvars_full = mbdimT*nbdimT;
% If the matrix Q is symmetric, the decision variables on the upper
% triangular half must equal those on the lower triangular have. Therefore,
% it suffices to construct only decision variables for one half of the
% matrix:
if ~symtoggle
    ndvars_sym = ndvars_full;   % no symmetry --> need to construct a full set of variables
else
    ndvars_sym = (mbdimT+1)*0.5*nbdimT; % symmetry --> construct only variables on lower (upper) triangular half
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now add the necessary information to the SOS Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add new variable
sos.var.num = sos.var.num+1;
var = sos.var.num;
sos.var.Z{var} = []; %makesparse(Z);
if postoggle
    sos.var.type{var} = 'sos';
    %        [T,ZZ] = getconstraint(Z);
    %        sos.var.ZZ{var} = ZZ;
    %        sos.var.T{var} = T';
    sos.var.ZZ{var} = [];
    sos.var.T{var} = [];
else
    sos.var.type{var} = 'poly';
    sos.var.ZZ{var} = [];%makesparse(Z); %is this correct? -MP
    sos.var.T{var} = [];%speye(size(Z,1));
end
% Specify the size of Q, without exploiting potential symmetry!
sos.var.idx{var+1} = sos.var.idx{var} + ndvars_full;

% Modify existing constraints to make room for new dvars, not adjusting for
% any symmetry!
for i = 1:sos.expr.num
    sos.expr.At{i} = [sos.expr.At{i}; ...
        sparse(ndvars_full,size(sos.expr.At{i},2))];
end

% Modify existing objective to make room for new dvars at end of
% decision variable vector
sos.objective = [sos.objective; sparse(ndvars_full,1)];

% Finally, construct the names of all new decision variables in Q,
% exploiting symmetry
dvars_new = fastint2str(sos.var.idx{var}-sos.var.idx{1}+(1:ndvars_sym)');
% For each element i,j of Q, provide an index lind_dvars(i,j) of the
% corresponding decision variable. 
if ~symtoggle
    % If Q is not symmetric, indices are just 1 through ndvars_full
    lind_dvars = reshape(1:ndvars_full,mbdimT,nbdimT);
    sos.decvartable = [sos.decvartable; dvars_new];
else
    % If Q is symmetric, construct a matrix of which the values on the
    % lower triangular half correspond to indices of decision variables.
    % Then, replace the upper triangular half with a flipped version of
    % the lower triangular half:
    lind_dvars = (1:mbdimT)' + (0:mbdimT:(nbdimT-1)*mbdimT) - cumsum(0:1:nbdimT-1);
    lind_dvars = tril(lind_dvars)+tril(lind_dvars,-1)';
    % Add the unique list of dvars to the decvartable, repeating variables
    % that appear in both the lower and upper triangular half
    sos.decvartable = [sos.decvartable; dvars_new(lind_dvars(:))];
end
% Divide the matrix of indices in accordance with the number of monomials
% Z1c{i} and Z2c{j}
lindT_cell_mat = mat2cell(lind_dvars,mbdim,nbdim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Done with SOSPROGRAM, now to construct the corresponding dpvar P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over the monomial cells Z1c{i}, Z2c{j}, and construct the
% associated dpvar P{i,j}
P = cell(nZ1c,nZ2c);
if nargout>=3
    Q = cell(nZ1c,nZ2c);
end
disp_warning = 1;   % Display warning if the specified basis functions are inappropriate
for cell_idx=1:(nZ1c*nZ2c)
    [i,j] = ind2sub([nZ1c,nZ2c],cell_idx);
    if nnz(Z1c{i}.coeff)==nZ1(i) && nnz(Z2c{j}.coeff)==nZ2(j) && ...
            all(any(Z1c{i}.coeff,1)) &&  all(any(Z2c{j}.coeff,1))
        % % Basis functions are (scaled) monomials
        
        % First, rearrange the degmat of Z1{i} so that row p in the
        % degmat corresponds to the pth monomial:
        % Z1{i}(p) ~ Z1c{i}.coeff(p,:)*Z1c{i}.degmat(:,:)
        %          = Z1i_coeffs(p)*Z1i_degmat(p,:)
        Z1i_varname = Z1c{i}.varname;
        Z1i_degmat = Z1c{i}.degmat;
        [tmp_i,tmp_j,Z1i_coeffs] = find(Z1c{i}.coeff);
        Z1i_degmat(tmp_j,:) = Z1i_degmat(tmp_i,:);  % Sort the degmat in same order as the basis functions
        Z1i_coeffs(tmp_j) = Z1i_coeffs;             % Sort coefficients accordingly
        
        % Similarly, rearrange the degmat of Z2{j} so that row q in the
        % degmat corresponds to the qth monomial:
        Z2j_degmat = Z2c{j}.degmat;
        Z2j_varname = Z2c{j}.varname;
        [tmp_i,tmp_j,Z2j_coeffs] = find(Z2c{j}.coeff);
        Z2j_degmat(tmp_j,:) = Z2j_degmat(tmp_i,:);  % Sort the degmat in accordance with the monomials
        Z2j_coeffs(tmp_j) = Z2j_coeffs;             % Sort coefficients accordingly
        
        % Next, extract the decision variables in Q{i,j}
        % - lindT_vec(k) gives the index of decision variable Q{i,j}(k),
        %       for linear position k\in\{1,...,mbdim(i)*nbdim(j)}
        % - lindL_vec is the list of unique decision variable indices
        % - lindL_I is the index vector such that
        %       lindL_vec(lindL_I) = lindT_vec
        % - dvarnamesL is a list of unique decision variables in Q{i,j}(k)
        lindT_mat = lindT_cell_mat{i,j};
        lindT_vec = lindT_mat(:);
        [lindL_vec, ~,lindL_I] = unique(lindT_vec);
        ndvarsL = length(lindL_vec);
        dvarnamesL = dvars_new(lindL_vec);
        
        
        % Now, for each of the mbdim(i)*nbdim(j) elements of Q{i,j},
        % establish associated row and column positions in Q{i,j}
        % Q{i,j}(k) = Q{i,j}(IposQij(k),JposQij(k))
        IposQij = repmat((1:mbdim(i))',nbdim(j),1);
        JposQij = kron((1:nbdim(j))',ones(mbdim(i),1));
        
        % For each element k of Q{i,j}, establish also which position
        % P{i,j}(Ipos(k),Jpos(k)) it map to...
        Ipos = ceil(IposQij./nZ1(i));
        Jpos = ceil(JposQij./nZ2(j));
        % ... as well as which monomials Z1c{i}(Z1ind(k)) and
        % Z2c{j}(Z2ind(k)) it maps:
        Z1ind = mod(IposQij-1,nZ1(i))+1;
        Z2ind = mod(JposQij-1,nZ2(j))+1;
        
        % Establish also what factor the monomials are multiplied with:
        valCij = repmat(Z1i_coeffs,[mdim_in(i)*nbdim(j),1]).*repmat(kron(Z2j_coeffs,ones(mbdim(i),1)),[ndim_in(j),1]);
        
        % % We now have
        % % P{i,j}(Ipos(k),Jpos(k)) = valCij(k) * Z1c{i}(Z1ind(k)) * Q{i,j}(k) * Z2c{j}(Z2ind(k))
        
        % We still need to sync independent variables in Z1{i} and Z2{j}
        [BvarnameL] = union(Z1i_varname,Z2j_varname);
        nvarsL = length(BvarnameL);
        % Reorganize columns of Z1i_degmat to match new varnames
        [~,idx1a,idx2a] = intersect(Z1i_varname,BvarnameL);
        Z1T = spalloc(nZ1(i),nvarsL,nnz(Z1i_degmat));
        Z1T(:,idx2a) = Z1i_degmat(:,idx1a);
        % Reorganize columns of Z2i_degmat to match new varnames
        [~,idx1b,idx2b] = intersect(Z2j_varname,BvarnameL);
        Z2T = spalloc(nZ2(j),nvarsL,nnz(Z2j_degmat));
        Z2T(:,idx2b) = Z2j_degmat(:,idx1b);
        
        % Next, create an amalgamated list of monomials Z12T = Z2 otimes Z1
        Z12T = kron(Z2T,ones(nZ1(i),1))+sprepmat(Z1T,nZ2(j),1);
        % in this unreduced list of nZ1*nZ2 terms, we have
        % Z1c{i}(p) * Z2c{j}(q) = Z12T(p+(q-1)*nZ1). Defining
        Z12Tind = Z1ind + (Z2ind-1)*nZ1(i);
        % we then have
        % Z1c{i}(Z1ind(k)) * Z2c{j}(Z2ind(k)) = Z12T(Z12Tind(k))
        % Get rid of redundant monomials...
        [ZT,~,IC] = unique(Z12T,'rows');    % ZT(IC,:) = Z12T
        nZT = size(ZT,1);
        % ... and update the indices accordingly
        ZTind = IC(Z12Tind);
        
        % Now, Z12T(k) = ZT(IC(k)), so we have
        % P{i,j}(Ipos(k),Jpos(k)) = valCij(k) * Z1c{i}(Z1ind(k)) * Q{i,j}(k) * Z2c{j}(Z2ind(k))
        %                         = valCij(k) * Q{i,j}(k) * Z12T(Z12Tind(k))
        %                         = valCij(k) * Q{i,j}(k) * ZT(ZTind)
        
        % Now we construct the position in the C{i,j} matrix for each
        % decision variable k
        IindCij = lindL_I + 1 + (Ipos-1)*(ndvarsL+1);
        JindCij = ZTind + (Jpos-1)*nZT;
        
        % With that, we can build the coefficient matrix, and the
        % actual dpvar object
        C = sparse(IindCij,JindCij,valCij,mdim_in(i)*(ndvarsL+1),ndim_in(j)*nZT);
        dvarname = dvarnamesL;
        varname = BvarnameL;
        degmat = ZT;
        P{i,j} = dpvar(C,degmat,varname,dvarname,[mdim_in(i) ndim_in(j)]);
        
        % If desired, also output the names of the decision variables
        % in Q{i,j} parameterizing P{i,j}
        if nargout>=3
            Q{i,j} = dvars_new(lindT_mat);
        end
        
    else
        % % Z1c and/or Z2c are arbitrary polynomial basis functions
        
        % Extract the variable names and degree matrices
        Z1i_varname = Z1c{i}.varname;
        Z1i_degmat = Z1c{i}.degmat;
        ndZ1 = size(Z1i_degmat,1);
        Z2j_degmat = Z2c{j}.degmat;
        Z2j_varname = Z2c{j}.varname;
        ndZ2 = size(Z2j_degmat,1);
        ndZT = ndZ1*ndZ2;
        
        % If the basis functions are composed of multiple monomials, we may
        % end up with more monomials than decision variables, potentially
        % resulting in an ill-posed problem...
        if disp_warning && postoggle && (ndZ1 > nZ1(i) || ndZ2 > nZ2(j))
            warning(['The number of monomials defining your positive (SOS) operator appears to exceed ',...
                     'the number of decision variables parameterizing it; your problem is likely ill-posed. ',...
                     'Consider just using (scaled) monomials as your SOS basis functions.'])
            disp_warning = 0;
        end
        
        % Next, extract the coefficients, and consider all possible
        % combinations of monomials describing Z1c{i} and Z2c{j}
        Z1i_coeffs = Z1c{i}.coeff;
        Z1i_coeffs = sprepmat(Z1i_coeffs,[ndZ2,mdim_in(i)*nbdim(j)]);
        Z1i_coeffs = reshape(Z1i_coeffs,[],1);
        
        Z2j_coeffs = Z2c{j}.coeff;
        Z2j_coeffs = sprepmat(kron(Z2j_coeffs,ones(ndZ1,mbdim(i))),[1,ndim_in(j)]);
        Z2j_coeffs = reshape(Z2j_coeffs,[],1);
        
        % Compute the product of the coefficents for each combination
        % of monomials.
        % Element (p1 + p2*nZ1 + k1*nZ1*nZ2 + k2*nZ1*nZ2*nZ1 + l2*nZ1*nZ2*nZ1*mdim_in(1) + l2**nZ1*nZ2*nZ1*mdim_in*nZ2)
        % corresponds to contribution of monomial p1 of function Z1(k1)
        % and monomial p2 of Z2(k2) to element (l1,l2) of the
        % matrix-valued P{i,j}
        valCij = Z1i_coeffs.*Z2j_coeffs;
        
        % Next, extract the decision variables in Q{i,j}
        % - lindT_vec(k) gives the index of decision variable Q{i,j}(k),
        %       for linear position k\in\{1,...,mbdim(i)*nbdim(j)}
        % - lindL_vec is the list of unique decision variable indices
        % - lindL_I is the index vector such that
        %       lindL_vec(lindL_I) = lindT_vec
        % - dvarnamesL is a list of unique decision variables in Q{i,j}(k)
        lindT_mat = lindT_cell_mat{i,j};
        lindT_vec = lindT_mat(:);
        [lindL_vec, ~,lindL_I] = unique(lindT_vec);
        ndvarsL = length(lindL_vec);
        dvarnamesL = dvars_new(lindL_vec);
        % NOTE: Although each index in lindL_I is linked to a
        % particular combination of basis functions
        % Z1c{i}(p)*Z2c{j}(q), each basis function may be comprised of
        % nZ1 and respectively nZ2 monomials. To account for this, we
        % will later repeat each index nZ1*nZ1 times to account for
        % every possible combination of monomials for each decision
        % variable.
        
        % Now, for each of the mbdim(i)*nbdim(j) elements of Q{i,j},
        % establish associated row and column positions in Q{i,j}
        % Q{i,j}(k) = Q{i,j}(IposQij(k),JposQij(k))
        IposQij = repmat((1:mbdim(i))',nbdim(j),1);
        JposQij = kron((1:nbdim(j))',ones(mbdim(i),1));
        
        % For each element k of Q{i,j}, establish also which position
        % P{i,j}(Ipos(k),Jpos(k)) it map to.
        Ipos = ceil(IposQij./nZ1(i));
        Jpos = ceil(JposQij./nZ2(j));
        % NOTE: Each basis function is still composed of multiple
        % monomials! To account for this, we will later repeat each of
        % the matrix positions nZ1*nZ2 times to account for each
        % possible combination of monomials contributing to each
        % position of the matrix P{i,j}
        
        % Now, we focus on the monomials Z1m, Z2m defining the basis
        % functions Z1{i} and Z2{j}.
        % First, we sync independent variables in Z1m and Z2m
        [BvarnameL] = union(Z1i_varname,Z2j_varname);
        nvarsL = length(BvarnameL);
        % Reorganize columns of Z1i_degmat to match new varnames
        [~,idx1a,idx2a] = intersect(Z1i_varname,BvarnameL);
        Z1T = spalloc(ndZ1,nvarsL,nnz(Z1i_degmat));
        Z1T(:,idx2a) = Z1i_degmat(:,idx1a);
        % Reorganize columns of Z2i_degmat to match new varnames
        [~,idx1b,idx2b] = intersect(Z2j_varname,BvarnameL);
        Z2T = spalloc(ndZ2,nvarsL,nnz(Z2j_degmat));
        Z2T(:,idx2b) = Z2j_degmat(:,idx1b);
        
        % Next, create an amalgamated list of monomials
        % Z12T = Z2m otimes Z1m
        Z12T = kron(Z2T,ones(ndZ1,1))+sprepmat(Z1T,ndZ2,1);
        % in this unreduced list of ndZ1*ndZ2 terms, we have
        % Z1m(p) * Z2m(q) = Z12T(p+(q-1)*ndZ1). Defining
        Z1ind = repmat((1:ndZ1)',[ndZ2,1]);
        Z2ind = kron((1:ndZ2)',ones(ndZ1,1));
        Z12Tind = Z1ind + (Z2ind-1)*nZ1(i);
        % we then have
        % Z1m(Z1ind(k)) * Z2m(Z2ind(k)) = Z12T(Z12Tind(k))
        % Get rid of redundant monomials...
        [ZT,~,IC] = unique(Z12T,'rows');    % ZT(IC,:) = Z12T
        nZT = size(ZT,1);
        % ... and update the indices accordingly
        ZTind = IC(Z12Tind);
        % Then Z1m(Z1ind(k)) * Z2m(Z2ind(k)) = ZT(ZTind(k))
        
        % Now we determine which row in the coefficient matrix C is
        % associated to each particular decision variable and each
        % particular row position in the matrix P{i,j}. Note that each
        % index is repeated ndZT=nZ1*nZ2 times, to account for the fact
        % that every possible combination of monomials in Z1m and Z2m might
        % appear in the product Z1c{i}(p)*Z2c{i}(q) associated to this
        % particular decision variable and position in P{i,j}.
        IindCij = kron(lindL_I + 1 + (Ipos-1)*(ndvarsL+1),ones(ndZT,1));
        % Similarly, we determine which column in the coefficient
        % matrix C is associated to each particular monomial in ZT and
        % each particular column position in the matrix P{i,j}. Note
        % that the monomial indices are repeated mbdim*nbdim times,
        % accounting for each possible combination of basis function
        % Z1c{i}(p) and Z2c{i}(q) that we map, and the column positions
        % are repeated ndZT=nZ1*nZ2 times, accounting for every possible
        % combination of the monomials appearing in these basis functions.
        JindCij = repmat(ZTind,[mbdim(i)*nbdim(j),1]) + kron((Jpos-1)*nZT,ones(ndZT,1));
        
        % With that, we can build the coefficient matrix, and the
        % actual dpvar object
        C = sparse(IindCij,JindCij,valCij,mdim_in(i)*(ndvarsL+1),ndim_in(j)*nZT);
        dvarname = dvarnamesL;
        varname = BvarnameL;
        degmat = ZT;
        P{i,j} = dpvar(C,degmat,varname,dvarname,[mdim_in(i) ndim_in(j)]);
        
        % If desired, also output the names of the decision variables
        % in Q{i,j} parameterizing P{i,j}
        if nargout>=3
            Q{i,j} = dvars_new(lindT_mat);
        end
        
    end
end
if celltoggle==0
    P = P{1,1};
    if nargout>=3
        Q = Q{1,1};
    end
end
    
    
end