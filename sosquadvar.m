function [sos,P] = sosquadvar(sos,Z1c,Z2c,mdim_in,ndim_in,option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [sos,P] = sosquadvar(sos,Z1c,Z2c,mdim_in,ndim_in,option)
%
% This program declares a symbolic positive scalar semidefinite matrix P of
% size nxn which is positive semidefinite for all values of the variables.
% The matrix has the form
% P=[Z1{1}    ]^T    [Z2{1}     ]
%   [ Z1{1}   ]    Q [ Z2{1}    ]
%   [  Z1{2}  ]      [  Z2{1}   ]
%   [   Z1{2} ]      [   Z2{2}  ], Q>0
%   [    Z1{2}]      [    Z2{2} ]
%                    [     Z2{2}]
%
%Where $Z$ is a vector of monomials. if option='pos', dims is square and Z1
% and Z2 are of equal length, then Q is a positive semidefinite matrix.
%
% NONCELLULAR INPUT FORMAT
% INPUTS:
% prog - The SOS program to which to attach the variable
% Z1c - column vector of monomials (pvar) to be multiplied on the left
% Z2c - column vector of monomials (pvar) to be multiplied on the right
% mdim_in - a scalar indicating the row dimension of the output.
% If this entry is left empty, it defaults to 1
% ndim_in - a scalar indicating the column dimension of the output.
% If this entry is left empty, it defaults to 1
% option - 'pos' if Q should be a positive semidefinite matrix
%        - 'sym' if Q should be a symmetric matrix (note this may not result in P being symmetric)
% In order to use the positive or symmtric options, the length of the
% monomial vectors Z1 and Z2 should be identical. In addition, we require
% mdim_in=ndim_in.
%
% OUTPUTS:
% prog - updated program structure
% P - dpvar polynomial structure. This structure is not cellular unless Z1c
% or Z2c were cellular
%
% ALTERNATIVE CELLULAR INPUT FORMAT
% For a multipartite structure, Z1c and Z2c may be cells of monomial column
% vectors. The number of cells in Z1c and Z2c need not be the same unless
% the 'pos' or 'sym' options are used. Note that mdim_in and ndim_in are
% now expected to be vectors indicating the number of rows and columns to
% be output for each of the cells in Z1c and Z2c. 
%
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
% SAFETY FEATURE. If you encounter an error or unexpected result, set
% toggle=1 in order to call a simpler and presumably more reliable, but slower construction of the
% variable.
%
%
% NOTE: REVERSIONARY

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log and developer notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 07/09/21, DJ: split dims into mdim and ndim, specifying number of blocks
% for each monomial Z1c{i} and Z2c{j} respectively. Also separate diagonal
% and off-diagonal computations under toggle==1
% 07/12/21, DJ: Added code to adjust the order of the rows in the degmats
% of Z1c{i} and Z2c{j} to match the rows in the associated monomials
% 07/13/21, DJ: Adjusted specification of dvars to use int2coeff function
% 07/14/21, DJ: Adjusted specification of dvars to use fastint2str function

toggle=2;

if nargin<=2
    error('quadvar requires at least 3 inputs')
elseif nargin==3
    mdim_in = 1;
    ndim_in = 1;
end

if nargin>=4
    if isempty(mdim_in)
        mdim_in = 1;
    end
    if isempty(ndim_in)
        ndim_in = 1;
    end
end


celltoggle=1;
if ~iscell(Z1c) && ~iscell(Z1c)
    celltoggle=0;
end

if ~iscell(Z1c)
    Z1c={Z1c};
end
if ~iscell(Z2c)
    Z2c={Z2c};
end

if isa(Z1c{1},'sym') || isa(Z2c{1},'sym')
    error('quadvar is incompatible with symbolic variables. Please convert your monomials to polynomial/pvar.')
end

nZ1c=length(Z1c);
nZ2c=length(Z2c);
nZ1=zeros(nZ1c,1);
nZ2=zeros(nZ2c,1);
%mbdim=zeros(1,nZ1c);
%nbdim=zeros(1,nZ2c);

%mdim_in = reshape(mdim_in,1,[]);
%ndim_in = reshape(ndim_in,1,[]);

mdim_in = reshape(mdim_in,[],1);
ndim_in = reshape(ndim_in,[],1);

if length(mdim_in)==1
    if nZ1c>=2
        mdim_in=mdim_in*ones(nZ1c,1);
    end
elseif length(mdim_in)~=nZ1c
    error('length of mdim_in must match number of cells in Z1')
end

if length(ndim_in)==1
    if nZ2c>=2
        ndim_in=ndim_in*ones(nZ2c,1);
    end
elseif length(ndim_in)~=nZ2c
    error('length of ndim_in must match number of cells in Z2')
end


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


postoggle=0;
symtoggle=0;
if nargin>3
    if length(mdim_in)==1
        mdim_in = mdim_in*ones(nZ1c,1);
    elseif length(mdim_in)~=nZ1c
        error('mdim declaration is invalid')
    end
    if length(ndim_in)==1
        ndim_in = ndim_in*ones(nZ2c,1);
    elseif length(ndim_in)~=nZ2c
        error('ndim declaration is invalid')
    end
    if nargin==6
        if option == 'pos'
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
        elseif option == 'sym'
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
else
    mdim_in = ones(nZ1c,1);
    ndim_in = ones(nZ2c,1);
end
%mdim = sum(mdim_in);
%ndim = sum(ndim_in);


mbdim=mdim_in.*nZ1;
nbdim=ndim_in.*nZ2;


mbdimT=sum(mbdim);
nbdimT=sum(nbdim);
ndvars=mbdimT*nbdimT;

if toggle==1
    BZ1=cell(1,nZ1c);
    BZ2=cell(1,nZ2c);
    P=cell(nZ1c,nZ2c);
    for i=1:nZ1c
        BZ1{i}=kron(eye(mdim_in(i)),Z1c{i});
    end
    for i=1:nZ2c
        BZ2{i}=kron(eye(ndim_in(i)),Z2c{i});
    end
    if postoggle
        [sos, Q]=DPsosposmatr(sos,mbdimT);
        for i=1:nZ1c
            Qtemp=Q(1+sum(mbdim(1:i-1)):sum(mbdim(1:i)),1+sum(nbdim(1:i-1)):sum(nbdim(1:i)));
            P{i,i}=BZ1{i}.'*Qtemp*BZ2{i};
            for j=i+1:nZ2c
                Qtemp=Q(1+sum(mbdim(1:i-1)):sum(mbdim(1:i)),1+sum(nbdim(1:j-1)):sum(nbdim(1:j)));
                P{i,j}=BZ1{i}.'*Qtemp*BZ2{j};
                P{j,i}=BZ1{j}.'*Qtemp.'*BZ2{i};
                %                P{j,i}=P{i,j}.'; % need a varswap here? at least for poslpivar, yes
            end
        end
    else
        if symtoggle
            [sos, Q]=DPsospolyvar_mat(sos,1,[mbdimT nbdimT],'symmetric');
            for i=1:nZ1c
                for j=1:nZ2c
                    P{i,j}=BZ1{i}.'*Q(1+sum(mbdim(1:i-1)):sum(mbdim(1:i)),1+sum(nbdim(1:j-1)):sum(nbdim(1:j)))*BZ2{j};
                end
            end
        else
            [sos, Q]=DPsospolyvar_mat(sos,1,[mbdimT nbdimT]);
            for i=1:nZ1c
                for j=1:nZ2c
                    P{i,j}=BZ1{i}.'*Q(1+sum(mbdim(1:i-1)):sum(mbdim(1:i)),1+sum(nbdim(1:j-1)):sum(nbdim(1:j)))*BZ2{j};
                end
            end
            
        end
    end
    if celltoggle==0
        P=P{1,1};
    end
end

if toggle==2
    
    % Modify SOSPROGRAM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% now add the necessary information to the SOS Program
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Disabling all the extra variable information until we decide how this should be handled.
    %    [~,idx1,idx2] = intersect(Bvarname,sos.vartable);
    %    Z = sparse(size(Z12T,1),length(sos.vartable));
    %    Z(:,idx2) = sparse(Z12(:,idx1));
    %    Z = sparse(nZT,length(sos.vartable));
    %    Z(:,idx2) = sparse(ZT(:,idx1));
    %    lenZ = size(Z,1);
    
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
    sos.var.idx{var+1} = sos.var.idx{var}+ndvars;
    
    % Construct the names of all new decision variables in Q
    dvars_new = fastint2str(sos.var.idx{var}-sos.var.idx{1}+(1:ndvars)');
    sos.decvartable = [sos.decvartable; dvars_new];
    
    % Modify existing constraints to make room for new dvars at end of
    % decision variable vector
    for i = 1:sos.expr.num
        sos.expr.At{i} = [sos.expr.At{i}; ...
            sparse(ndvars,size(sos.expr.At{i},2))];
    end
    
    % Modify existing objective to make room for new dvars at end of
    % decision variable vector
    sos.objective = [sos.objective; sparse(ndvars,1)];
    
    %%%%%%%%%%%%%%%%%%%
    % done with SOSPROGRAM, now to construct the corresponding dpvar
    %%%%%%%%%%%%%%%%%%%
    % list of decision variable numbers in the SOSPROGRAM
    lindSOS = 1:ndvars;
    
    % time to map decision variable number to monomial and i,j coordinate.
    % These will then be use to map to positions in the C matrix
    % This is slightly complicated by the fact that if we are using a
    % symmetric matrix of decision variables and only want to retain the
    % upper triangular part.
    if symtoggle
        % this is the matrix of decision variable numbers in Q
        temp=reshape(lindSOS,mbdimT,nbdimT);
        % this is the
        %partition into submatrices of SOS decision variable numbers in the shape of Q{i,j}
        % since this is a symmetric matrix, we reassign the lower triangular variables to the upper-triangular ones.
        lindT_cell_mat=mat2cell(triu(temp)+triu(temp,1)',mbdim,nbdim);
    else
        %partition into submatrices of SOS decision variable numbers in the shape of Q{i,j}
        lindT_cell_mat=mat2cell(reshape(lindSOS,mbdimT,nbdimT),mbdim,nbdim);

    end
    
    
    
    
    
    % apparently foor loops are faster than cellfuns... lame
    P=cell(nZ1c,nZ2c);
    for i=1:nZ1c
        % make sure the rows in degmat match the rows in the monomial Z1{i}
        deg1_new = Z1c{i}.coef\Z1c{i}.degmat;
        Z1i = polynomial(speye(nZ1(i)),deg1_new,Z1c{i}.varname,[nZ1(i),1]);
        for j=1:nZ2c
            % make sure the rows in degmat match the rows in the monomial Z2{j}
            deg2_new = Z2c{j}.coef\Z2c{j}.degmat;
            Z2j = polynomial(speye(nZ2(j)),deg2_new,Z2c{j}.varname,[nZ2(j),1]);
            
            % Alternative:
            %Z2j = Z2c{j};
            %[row2,col2] = find(Z2j.coef);
            %indcs2 = [row2,col2];
            %indcs2 = sortrows_integerTable(indcs2);
            %col2 = indcs2(:,2);
            %[~,new_order] = sortrows(col2);
            %deg2_new = Z2j.degmat(new_order,:);
            %Z2j = polynomial(speye(nZ2(j)),deg2_new,Z2c{j}.varname,[nZ2(j),1]);
            
            % preprocessing the monomials to be used in each cell of P{i,j}
            lindT_mat=lindT_cell_mat{i,j};
            lindT_vec= lindT_mat(:);
            [lindL_vec, ~,lindL_I]=unique(lindT_vec);
            % lindL_vec is the list of decision variable numbers
            % lindL_I is the vector of local variable index numbers in
            % LindL_vec corresponding to lindT_vec
            ndvarsL=length(lindL_vec);
            % list of local decision variable positions
            %            temp=(1:ndvarsL);
            % ordered list of decision variable names
            dvarnamesL=dvars_new(lindL_vec);
            
            
            
            % List of (I,J) position numbers corresponding to IC
            IposQij = repmat((1:mbdim(i))',nbdim(j),1);
            JposQij = kron((1:nbdim(j))',ones(mbdim(i),1));
            %            JposQij = kron(ones(mbdim(i),1),(1:nbdim(j))');
            
            % Go modular arithmatic!
            % now convert Iind and Jind to i,j position in P{i,j} matrix
            Ipos=ceil(IposQij./nZ1(i));
            Jpos=ceil(JposQij./nZ2(j));
            % now convert Iind and Jind to monomial numbers
            Z1ind=mod(IposQij-1,nZ1(i))+1;
            Z2ind=mod(JposQij-1,nZ2(j))+1;
            
            
            % need to sync independent variables in Z1{i} and Z2{j}
            [BvarnameL] = union(Z1i.varname,Z2j.varname);
            nvarsL=length(BvarnameL);
            [~,idx1a,idx2a] = intersect(Z1i.varname,BvarnameL);
            Z1T = sparse(nZ1(i),nvarsL);
            Z1T(:,idx2a) = Z1i.degmat(:,idx1a);
            %                Z1T(:,idx2a) = sparse(Z1{i}.degmat(:,idx1a));
            %    Z1T = sparse(Z1.degmat(:,idx1));
            
            [~,idx1b,idx2b] = intersect(Z2j.varname,BvarnameL);
            Z2T = sparse(nZ2(j),nvarsL);
            Z2T(:,idx2b) = Z2j.degmat(:,idx1b);
            %   Z2T = sparse(Z2.degmat(:,idx1b));
            
            
            %    creating amalgamated list of monomials using Z2 otimes Z1
            Z12T = kron(Z2T,ones(nZ1(i),1))+sprepmat(Z1T,nZ2(j),1);
            % in this unreduced list of nZ1*nZ2 terms, the position of monomial
            % Z12T(i+(j-1)*nZ1)=Z1(i)Z2(j) - the other direction of conversion will
            % not be necessary
            % create list of unique monomials
            [ZT,~,IC] = unique(Z12T,'rows');    % ZT(IC,:) = Z12T
            %    [ZT,IA,IC] = unique(Z12T,'rows');
            nZT=size(ZT,1);
            % Now, Z12T(i)=ZT(IC(i)) and ZT(i)=Z12T(IA(i)) so we have
            % ZT(IC(i+(j-1)*nZ1))= Z12T(i+(j-1)*nZ1)=Z1(i)Z2(j)
            
            % now convert Z2 and Z1 index to Ztot index
            ZTind=IC(Z1ind+(Z2ind-1)*nZ1(i)); % this is the monomial in ZT associated with
            % Now we construct the position in the C{i,j} matrix for each decision
            % variable
            IindCij=lindL_I+1+(Ipos-1)*(ndvarsL+1);
            JindCij=ZTind+(Jpos-1)*nZT;
            C=sparse(IindCij,JindCij,1,mdim_in(i)*(ndvarsL+1),ndim_in(j)*nZT,mbdim(i)*nbdim(j));
            
            % Ordered names of decision vars which actually will appear in the expression:
            dvarname=dvarnamesL;%dvars_new(dvarlist{i,j});
            varname=BvarnameL;
            degmat=ZT;
            P{i,j} = dpvar(C,degmat,varname,dvarname,[mdim_in(i) ndim_in(j)]);
            
        end
    end
    if celltoggle==0
        P=P{1,1};
    end
    % Error handling needs to be added here, e.g. if Z is not a
    % valid monomial vector.
    
    
end
