function [sos,P] = sospolyvar_mat(sos,ZSym,dims,matrixstr,wscoeff)
% SOSPOLYVAR_MAT --- Declare a new polynomial matrix variable in
%       an SOS program
%
% OUTPUTS:
% SOS: the modified sum of squares program structure.
% P: the new polynomial matrix variable. The output type is dpvar.
%
% INPUTS
% SOS: The sum of squares program structure to be modified.
% ZSym: The vector of monomials to be contained in P. Decision
% variables corresponding to those monomials will be assigned
% automatically by SOSPOLYVAR_mat.
% dims: The desired dimension of the output matrix P: n(1) x n(2)
% matrixstr (optional): a char string with the option 'symmetric' when required
% wscoeff (optional): If wscoeff='symmetric', sospolyvar will create
% the decision variables corresponding to VAR (i.e., coeff_xxx)
% also in MATLAB workspace.
%
% NOTE: this is an internal file which is only called for symmetric
% polynomial matrices when the output type is dpvar. This is to efficiently
% handle the case which sosquadvar misses.

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
% -MP 6/27/2021: updated to dpvar output only. Disabled use of symbolic variables.
% Also, removed for-loop in both coefficient name declaration and degmat matrix declaration
% DJ, 07/25/21: Small adjustment to fix issue with ndim=1 case

if nargin<3
    mdim=1;ndim=1;
else
    mdim=dims(1); ndim=dims(2);
end
symtoggle=0;
if nargin==4
    if strcmp(matrixstr,'symmetric')
        if ndim~=mdim
            error(['''symmetric''' ' option used, matrix must be square.'])
        else
            symtoggle=1;
        end
    end
    
end

if isnumeric(ZSym) & ZSym==1
    pvar ZSym;
    ZSym = 0*ZSym+1;
end
% locating independent varibles in SOS program structure
[~,idx1,idx2] = intersect(ZSym.varname,sos.vartable);
Z = sparse(size(ZSym.degmat,1),length(sos.vartable));
Z(:,idx2) = sparse(ZSym.degmat(:,idx1));
lenZ = size(Z,1);
if symtoggle==1
    if mdim~=ndim
        error('matrix dimension is not square. Symmetric matrices must be square')
    end
    %        ndvars=lenZ*(ndim^2-ndim);
    [lind_p] = find(kron(triu(ones(ndim)),ones(1,lenZ))); % positional index
    ndvars=length(lind_p);
    [Iind1,Jind1] = ind2sub([ndim,ndim],lind_p);
    Iind1 = reshape(Iind1,[],1);    Jind1 = reshape(Jind1,[],1); % for when ndim==1
    lind1=(1:ndvars)'; % these are the positions in dvarname vector corresponding to each linear index position
    % in the symmetric case, we have
    lind2=[lind1;lind1];
    
    Jind_t=ceil(Jind1./lenZ); % j position of lind
    monind=mod(Jind1-1,lenZ)+1; % My better attempt at modular arithmatic - monomial position of lind
    Iind2=[Iind1; Jind_t]; % I and j positions are swapper for second set of assignments
    Jind2=[Jind1; (Iind1-1)*lenZ+monind];
    % now calculate the corresponsing position in the C matrix
    IindC=lind2+1+(Iind2-1)*(ndvars+1);
    JindC=Jind2;
    %                [Iind1,Jind1] = find(kron(triu(ones(ndim)),ones(1,lenZ)));
    %                ndvars=length(lind1);
else
    ndvars=lenZ*mdim*ndim;
    % this is the position in the dvartable
    [lind1] = find(ones(mdim,ndim*lenZ));
    % this is the position
    [Iind1,Jind1] = ind2sub([mdim,ndim*lenZ],lind1);
    
    IindC=lind1+1+(Iind1-1)*(ndvars+1);
    JindC=Jind1;
end

% Error handling needs to be added here, e.g. if Z is not a
% valid monomial vector.

% Add new variable
sos.var.num = sos.var.num+1;
var = sos.var.num;
sos.var.type{var} = 'poly';
sos.var.Z{var} = makesparse(Z);
sos.var.ZZ{var} = makesparse(Z); %is this correct? -MP
sos.var.T{var} = speye(size(Z,1));
sos.var.idx{var+1} = sos.var.idx{var}+ndvars;

% Modify existing constraints to make room for new dvars at end of
% decision variable vector
for i = 1:sos.expr.num
    sos.expr.At{i} = [sos.expr.At{i}; ...
        sparse(ndvars,size(sos.expr.At{i},2))];
    %            sparse(size(sos.var.T{var},1),size(sos.expr.At{i},2))];
end

% Modify existing objective to make room for new dvars at end of
% decision variable vector
sos.objective = [sos.objective; sparse(ndvars,1)];        % 01/07/02
%    sos.objective = [sos.objective; sparse(sos.var.idx{var+1}-sos.var.idx{var},1)];        % 01/07/02


% eliminating for-loop for coefficient declaration
dvars_new=cellstr([repmat('coeff_',ndvars,1), strjust(int2str(sos.var.idx{var}-sos.var.idx{1}+(1:ndvars)'),'left')]); %SS
%   dvars_out=dvars_new(lind1);
sos.decvartable = [sos.decvartable; dvars_new];

% Option 1: seems to be fastest
%    C=sparse(2:1+lenZ,1:lenZ,ones(lenZ,1),lenZ+1,lenZ);

%Iinds=kron((0:mdim-1)*(1+lenZ*ndim),ones(1,lenZ*ndim))+repmat(2:1+lenZ*ndim,1,mdim);
%Jinds=repmat(1:lenZ*ndim,1,mdim);
if symtoggle==1
    C=spones(sparse(IindC,JindC,1,(ndvars+1)*mdim,lenZ*ndim,lenZ*ndim*mdim+lenZ*ndim));
else
    C=sparse(IindC,JindC,1,(ndvars+1)*mdim,lenZ*ndim,lenZ*ndim*mdim+lenZ*ndim);
end
% option 2:
%   C=sparse(1+lenZ,lenZ);
%   C(2:end,:)=speye(lenZ);
% option 3:
%   C=[zeros(1,lenZ); speye(lenZ)];
%    varname=sos.vartable;
varname=ZSym.varname;
degmat=ZSym.degmat;
dvarname=dvars_new;
P = dpvar(C,degmat,varname,dvarname,dims);

if nargin > 4 & strcmp(wscoeff,'wscoeff')
    var = sos.var.num;
    for i = sos.var.idx{var}:sos.var.idx{var+1}-1
        pvar(['coeff_',int2str(i-sos.var.idx{1}+1)]);
        assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],eval(['coeff_',int2str(i-sos.var.idx{1}+1)]));
    end
end


end
