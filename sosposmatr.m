function [sos,P] = sosposmatr(sos,n,sp_pat,wscoeff,PVoption)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [sos,P] = sosposmatr(sos,n,sp_pat,wscoeff,PVoption)
% This program declares a positive scalar (no independent variables)
% semidefinite matrix P of size nxn of the dpvar class format
%
% INPUTS:
% sos - The SOSprogram to which to append the variable
% n - size of symmetric matrix variable
% sp_pat - a desired sparsity pattern in the output matrix. sp_pat must be
% of size nxn with ones in places a variable is needed and 0's everywhere
% else. sp_pat may be a sparse matrix.
% wscoeff - if wscoeff='wscoeff', elements of the matrix will be declared
% as pvars in the workspace (not currently enabled)
% PVoption (reversionary) - if PVoption='pvar', the output is a pvar
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
% NOTES:
% MMP 6/2/2013 - restructured to reduce overhead. Eliminating the need for
% Z by direct assignment of variables without call to sossosvar.
% DJ, 07/19/21: Adjusted to use fastint2str for coefficient declaration

%%%%%%%%%%%%%%%%%% REVERSIONARY CODE %%%%%%%%%%%%%%%%%
if nargin==5 && strcmp(PVoption,'pvar')
    if ~isempty(sp_pat)
        sp_toggle=1;
        if ~all(size(sp_pat)==n)
            error('sparsity pattern must be nxn')
        end
        if ~all(all(triu(sp_pat)==tril(sp_pat)'))
            error('sparsity pattern must be symmetric')
        end
    else
        sp_toggle=0;
    end
    % Add new variable
    sos.var.num = sos.var.num+1;
    var = sos.var.num;
    sos.var.type{var} = 'sos';
    sos.var.Z{var} = [];%makesparse(Z);
    sos.var.ZZ{var} = [];
    sos.var.T{var} = [];
    sos.var.idx{var+1} = sos.var.idx{var}+n^2;
    
    
    % Augment existing equality constraints by adding new variables
    for i = 1:sos.expr.num
        sos.expr.At{i} = [sos.expr.At{i}; ...
            sparse(n^2,size(sos.expr.At{i},2))];
    end
    
    sos.objective = [sos.objective; sparse(sos.var.idx{var+1}-sos.var.idx{var},1)];        % 01/07/02
    
    if sp_toggle==0
        [lind1] = find(triu(ones(n)));
    else
        [lind1] = find(triu(sp_pat));
    end
    
    temp=cellstr(char(ones(n^2,1)));
    temp1=cellstr([repmat('coeff_',n^2,1), strjust(int2str(sos.var.idx{var}-sos.var.idx{1}+(1:n^2)'),'left')]); %SS
    temp=temp1(lind1);
    sos.decvartable = [sos.decvartable; temp1]; %  modified to add ALL decision variables to decvartable
    
    if nargin > 3 && strcmp(wscoeff,'wscoeff')
        var = sos.var.num;
        for i = sos.var.idx{var}:sos.var.idx{var+1}-1
            pvar(['coeff_',int2str(i-sos.var.idx{1}+1)]);
            assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],eval(['coeff_',int2str(i-sos.var.idx{1}+1)]));
        end
    end
    
    nterms=length(lind1);
    
%    begin_n=sos.var.idx{length(sos.var.idx)-1};
    nTT=n;
    
    matdim=[nTT,nTT];
    
%    degmat=speye((nTT^2/2+nTT/2));
 degmat=speye(nterms);
 
    varname = temp;
    [Ilist, Jlist] = ind2sub([n n], lind1);
    lind2=sub2ind([n n], Jlist, Ilist);

    coeff = spones(sparse(repmat([1:nterms]',2,1),[lind1;lind2],1,nterms,(n)^2));
    P=polynomial(coeff,degmat,varname,matdim);
else
    %%%%%%%%%%%%%%%%%% DPVAR CODE %%%%%%%%%%%%%%%%%
    if exist('sp_pat','var')
        sp_toggle=1;
        if ~all(size(sp_pat)==n)
            error('sparsity pattern must be nxn')
        end
        if ~all(all(triu(sp_pat)==tril(sp_pat)'))
            error('sparsity pattern must be symmetric')
        end
    else
        sp_toggle=0;
    end
    % Add new variable
    sos.var.num = sos.var.num+1;
    var = sos.var.num;
    sos.var.type{var} = 'sos';
    sos.var.Z{var} = [];%makesparse(Z);
    %    [T,ZZ] = getconstraint(Z);
    sos.var.ZZ{var} = [];
    sos.var.T{var} = [];
    sos.var.idx{var+1} = sos.var.idx{var}+n^2; % This declares all n^2 variables to exist, but only the non-zeros are sent to the workspace
    % P= [x1 0]   [y1 y2]
    %    [0 x2] = [y3 y4]>0
    % so that clearly x1=y1 and x2=y4. hence y1 and y4 are sent to the
    % workspace
    
    
    % Modify existing constraints to make room for new dvars at end of
    % decision variable vector
    for i = 1:sos.expr.num
        sos.expr.At{i} = [sos.expr.At{i}; ...
            sparse(n^2,size(sos.expr.At{i},2))];
    end
    
    % Modify existing objective to make room for new dvars at end of
    % decision variable vector
    sos.objective = [sos.objective; sparse(sos.var.idx{var+1}-sos.var.idx{var},1)];        % 01/07/02
    
    % Modify decision variable table
    %    oldlen = length(sos.decvartable);
    
    %
    % Now that the positive matrix has been added, we would like to only
    % declare the relevant decision variables. The problem is that the size of
    % the decision variable table should probably match the length of At(,2).
    % Now sossolve doesn't use the size of the decvartable. Neither does
    % sosgetsol.
    
    if sp_toggle==0
        [lind1] = find(triu(ones(n)));
    else
        [lind1] = find(triu(sp_pat));
    end
    
    %temp1=cellstr([repmat('coeff_',n^2,1), strjust(int2str(sos.var.idx{var}-sos.var.idx{1}+(1:n^2)'),'left')]); %SS
    %temp1 = int2coeff_sequence(sos.var.idx{var}-sos.var.idx{1}+1, sos.var.idx{var}-sos.var.idx{1}+n^2);
    temp1 = fastint2str(sos.var.idx{var}-sos.var.idx{1}+(1:n^2)');
    sos.decvartable = [sos.decvartable; temp1]; %  modified to add ALL decision variables to decvartable
    
    % decision vars which actually will appear in the expression:
    temp=temp1(lind1);
    
    % I should probably optimize this too. But really, who uses this option?
    if nargin > 3 & wscoeff == 'wscoeff'
        var = sos.var.num;
        for i = sos.var.idx{var}:sos.var.idx{var+1}-1
            pvar(['coeff_',int2str(i-sos.var.idx{1}+1)]);
            assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],eval(['coeff_',int2str(i-sos.var.idx{1}+1)]));
        end;
    end;
    
    % Lets try and create this directly
    % Lets eliminate uneccesary variables by enforcing the symmetry constraint.
    %
    
    [Ilist, Jlist] = ind2sub([n n], lind1);
    nterms=length(lind1);
    Ilist_temp=(Ilist-1)*(nterms+1)+[2:nterms+1]';
    Jlist_temp=(Jlist-1)*(nterms+1)+[2:nterms+1]';
    
    % This implementation is slower. Probably because of assignments
    % %lind4=sub2ind([n*nterms n], Jlist.*[1:nterms]', Ilist);
    % lind3=sub2ind([n*(nterms+1) n], (Ilist-1)*(nterms+1)+[2:nterms+1]', Jlist);
    % lind4=sub2ind([n*(nterms+1) n], (Jlist-1)*(nterms+1)+[2:nterms+1]', Ilist);
    % C=spalloc(n*(nterms+1),n,2*nterms-n);
    % C(lind4)=ones(nterms,1);
    % C(lind3)=ones(nterms,1);
    
    C= spones(sparse([Ilist_temp;Jlist_temp],[Jlist;Ilist],1,n*(nterms+1),n,n^2+n));
    
    tt=polynomial(1);
    varname=tt.varname;
    degmat=tt.degmat;
    dvarname = temp;
    matdim=[n,n];
    P=dpvar(C,degmat,varname,dvarname,matdim);
    
end