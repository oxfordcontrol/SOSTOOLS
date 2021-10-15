function [A,b,K,z,dv2x,Nfv,feas,zrem, removed_rows] = sospsimplify(A,b,K,z,dv2x,Nsosvarc, tol)
% function [A,b,K,z,dv2x,Nfv,feas,zrem, removed_rows] = sospsimplify(A,b,K,z,dv2x,Nsosvarc)
% 
% INPUT
% A - a constraint matrix that needs to be simplified (A*x = b)
% b - a constraint vector that needs to be simplified (A*x = b)
% K - describe the PSD cone in the Sedumi format
% z - vector of monomials
% dv2x - map between Sedumi variable and A matrix
% Nsosvarc - Number of PD constraints
% tol - desired tolerance for the simplification procedure
%
% OUTPUT
% A - simplified matrix A
% b - simplified vector b
% K - simplofied K
% dv2x - map reduced between environment variables and monomials
% Nfv  - the number of free variables
% feas - 0 if the problem is clearly infeasible, 1 otherwise
% zrem - removed monomials
% removed_rows - removed rows of matrix A
%
% DESCRIPTION 
%   This function performs a simplification procedure on the SOS problem.
%   First, it tries to detect the sign of optimization variables based on
%   simple constraints. Second, it searches for monomials that can be 
%   removed from each SOS constraint. This search is based on diagonal 
%   entries of the Gram matrix that are forced to be zero and it is 
%   equivalent to the Newton Polytope method. These two steps are repeated
%   until no new sign information can be detected.  The removed monomials
%   are stored in zrem.
%
%   In the code, the information about the sign of the optimization 
%   variables is stored in xsign where:
%      xsign(i)=NaN    if x(i) has unknown sign
%      xsign(i)=+1     if x(i)>=0
%      xsign(i)=-1     if x(i)<=0
%      xsign(i)=0      if x==0
%

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 4.00.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018, 2021  
%                                      A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5),
%                                      M. Peet (6), D. Jagt (6), A. Talitckii (6),
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
% 9/21/08  PJS Initial Coding
% 12/14/10 PJS Bug related to re-indexing of dv2x when removing free vars  
% 9/11/21 AT - add new output 'removed rows' to aid in reconstruction in SOSTOOLS 
%              removed equality constraints for remaining dec vars known to be zero
% 09/11/21 AT Bug related to A( : , abs(xsign)< tol ) = 0;
%             Change all element-wise operations to 'spfun' for sparse matrix


%--------------------------------------------------------------------
% Grab problem dimensions
%--------------------------------------------------------------------
Nfv = K.f;
Nlp = K.l; 
Nx = size(A,2);


% XXX -- We need a tolerance here.  How should we choose tol?
% Should It be the same tolerance as sedumi solver?. 
% tol = 1e-12;
%--------------------------------------------------------------------
% Find non-negative optimization variables
%--------------------------------------------------------------------
xsign = NaN([Nx,1]);

% Process LP constraints: A*d+y=b, y>=0
if Nlp>0
    % LP slack vars are >= 0
    xsign(Nfv+1:Nfv+Nlp) = +1;
end

% Process SOS Inequality Constraints: Ad*d + Aq*Q(:) = b and Q>=0
ptr = Nfv+Nlp;
for i1=1:length(K.s)
    % Diag entries of Q are >= 0    
    lz = K.s(i1);
    diagidx = (0:lz-1)*lz+(1:lz);
    xsign(ptr+diagidx) = +1;
    ptr = ptr+lz^2;
end

%--------------------------------------------------------------------
% Find dec vars and monomials that can be removed
%--------------------------------------------------------------------

go = 1;
zrem = cell(length(K.s),1);
xsignprev = xsign;
while (go == 1)
    % Use simple constraints to determine sign of optimization vars. 
    xsign = LOCALxsignupdate(xsign,A,b, tol);
    
    % AT: This line can make unfeasible problem to feasible, look at
    % test_dpvar_SOS_nonlinear_stability (n, deg, rng) = (3, 3, 70)
    %     A( : , abs(xsign)< tol ) = 0;
    
    % Monomial reduction procedure 
    % (This is equivalent to the Newton Polytope method)
    ptr = Nfv+Nlp;
    for i1=1:length(K.s)
        % Find diag entries of Q that are forced to be zero
        lz = K.s(i1);
%         if lz <15
%             continue
%         end
        blkidx = ptr+(1:lz^2);
        Qsign = diag(mat(xsign(blkidx)));
        loc = find(abs(Qsign)< tol);        
        if ~isempty(loc)
            % Corresponding rows/cols of zero diag entries are also zero
            tmp = sparse(lz,lz);
            tmp(loc,:) = 1;
            tmp(:,loc) = 1;
            rmidx = find(tmp)+ptr;
            xsign(rmidx) = 0;    
            
            % Remove vars/monoms associated with zero Gram matrix entries
            A(:,rmidx) = [];
%             AA(:,rmidx) = [];
            xsign(rmidx) = [];
            zrem{i1} = [zrem{i1}; z{i1}(loc)];
            z{i1}(loc) = [];
            %K.s(i1) = size( z{i1}, 1); %length( z{i1} );
            if isempty(z{i1})
                K.s(i1) = 0;
            else
                K.s(i1) = length( z{i1} );
            end
            
%             Update the mapping of dec vars into the optimization vars
            if i1<=Nsosvarc
                % Number of optim vars to remove
                Nremove = length(rmidx);
                
                % Optim vars currently in this block (before removal)
                blkidx = ptr+(1:lz^2);

                % Mark removed dec vars
                dv2x( ismember(dv2x,rmidx) ) = 0;
                          
                % Relabel the remaining dec vars in this block
                % AT 09/10/2021 sostools uses all variable and not just a
                % low triangle part like sosp
                
                idx = find( ones(K.s(i1)) );                            
                dv2x( ismember(dv2x,blkidx) ) = ptr+idx;
                                
                % Relabel remaining dec vars in subsequent blocks
                idx = find( dv2x>blkidx(end) );
                dv2x(idx) = dv2x(idx) - Nremove;                                               
            end        
        end
        
        % Update pointer
        ptr = ptr+K.s(i1)^2;
    end

    % Continue if xsign has been updated
    go = ~isequalwithequalnans(xsign,xsignprev);
    xsignprev = xsign;    
end

%--------------------------------------------------------------------
% Clean up 
%--------------------------------------------------------------------

% Mark removed free decision vars
rmidx = find( abs(xsign(1:K.f))< tol );
dv2x( ismember(dv2x,rmidx) ) = 0;
idx = find(dv2x<=K.f & dv2x>0);
A(:,rmidx) = [];
xsign(rmidx) = [];
Nremf = length(rmidx);
K.f = K.f - Nremf;  
Nfv = K.f;
dv2x(idx) = 1:Nfv;

idx2 = find(dv2x>K.f & dv2x>0);
dv2x(idx2) = dv2x(idx2)-Nremf;

% Remove any constraints of the form 0=0 
ridx = find( sum(abs(A)>tol,2)==0 & abs(b)<max(tol,tol*max(abs(b))) );
A(ridx,:) = [];

% just a new output
removed_rows = ridx;
b(ridx) = [];

% Check for infeasible problems of the form 0 = bi where bi is not equal 
% to zero (Our simplify code should flag infeasible problems because
% Sedumi can error out on problems that are trivially infeasible)
if isempty(A)
    feas = 1;
    return    
else
    ridx = find( sum(abs(A)>tol,2)==0 & abs(b)>tol*max(abs(b)) );
end
feas = 1;
if ~isempty(ridx)
    feas = 0;
end

% Add equality constraints for remaining dec vars known to be zero
% Now sossolve set the values equal 0.
% idx = find( xsign==0 );
% lidx = length(idx);
% A(end+1:end+lidx,idx) = speye(lidx);
% b = [b; sparse(lidx,1)];



%--------------------------------------------------------------------
% Local function to update sign of optimization var 
%--------------------------------------------------------------------
function xsign = LOCALxsignupdate(xsignOld,A,b, tol)
% AT: change a lot of elementwise operations to 
%     spfun(@function, sparse matrix)
% The tolerance is used for sign function 
% Now the sign function is used without tolerance
% Initialize output
xsign = xsignOld;

% Process constraints of the form:  aij*xj = bi
ridx = find(  sum(spones(A), 2)==1  );
if ~isempty(ridx)
    [cidx,tmp]=find( A(ridx,:)' );
    idx = sub2ind(size(A),ridx,cidx);
    % change element wise operation for sparse matrix.
    signA = spfun(@sign, A(idx) );
    signb = spfun(@sign, b(ridx) );    
    xsignUpdate = signA.*signb;
    xsignUpdate =  signb;
 
%     % XXX PJS 12/07/09: If cidx = [2;2] and xsignUpdate = [1; NaN] then 
%     % the next line will replace xsign(2) with NaN because the last index 
%     % in a subsasgn wins. This caused problems on a GSOSOPT problem.
%     %
%     % xsign(cidx) = LOCALupdate(xsign(cidx),xsignUpdate);
%     
%     % The correct code (also below) is below.  I'll try to vectorize
%     % if speed becomes an issue.
% %     for i1 =1:length(cidx)
% %         xsign(cidx(i1)) = LOCALupdate(xsign(cidx(i1)),xsignUpdate(i1));
% %     end

    xsign(cidx) = LOCALupdate(xsign(cidx),xsignUpdate);
end
% 
% % Process constraints of the form:  aij*xj + aik*xk = bi
ridx = find(  sum(spones(A), 2) == 2  );
if ~isempty(ridx)
    [cidx, tmp] = find( A(ridx,:)'  ); %remove tolerance, AT 10-15-21
    cidx  = reshape(cidx,[2 length(ridx)])';
    cidx1 = cidx(:,1);
    idx1  = sub2ind(size(A),ridx,cidx1);
    cidx2 = cidx(:,2);
    idx2  = sub2ind(size(A),ridx,cidx2);

    c1 = sign(b(ridx)./A(idx1));
    c2 = sign(-A(idx2)./A(idx1));
    xsignUpdate = NaN([length(ridx) 1]);
    xsignUpdate( c1<=0 & (c2.*xsign(cidx2)<=0) ) = -1;
    xsignUpdate( c1>=0 & (c2.*xsign(cidx2)>=0) ) = +1;
    
    xsign(cidx1) = LOCALupdate(xsign(cidx1),xsignUpdate);
    idx_update = ((xsignUpdate == 1) | (xsignUpdate == -1));
    xsign(cidx1(idx_update)) = LOCALupdate(xsign(cidx1(idx_update)),xsignUpdate(idx_update));
%     for i1 =1:length(cidx1)
%         if ((xsignUpdate == 1) | (xsignUpdate == -1))
%             xsign(cidx1(i1)) = LOCALupdate(xsign(cidx1(i1)),xsignUpdate(i1));
%         end
%     end
    
    c1 = sign(b(ridx)./A(idx2));
    c2 = sign(-A(idx1)./A(idx2));
    xsignUpdate = NaN([length(ridx) 1]);
    xsignUpdate( c1<=0 & (c2.*xsign(cidx1)<=0) ) = -1;
    xsignUpdate( c1>=0 & (c2.*xsign(cidx1)>=0) ) = +1;
    xsign(cidx2) = LOCALupdate(xsign(cidx2),xsignUpdate);
    idx_update = ((xsignUpdate == 1) | (xsignUpdate == -1));
    xsign(cidx1(idx_update)) = LOCALupdate(xsign(cidx1(idx_update)),xsignUpdate(idx_update));
%     for i1 =1:length(cidx2)
%         if ((xsignUpdate == 1) | (xsignUpdate == -1))
%             xsign(cidx1(i1)) = LOCALupdate(xsign(cidx1(i1)),xsignUpdate(i1));
%         end
%     end        
end

% % Process constraints of the form:  aij*xj + aik*xk + ail*xl= 0
% % where aij*xj, aik*xk, ail*xl all have the same sign.  
% % This implies that each of the three vars  = 0
ridx = find(  sum(spones(A),2)==3 & spones(b) < tol );
if ~isempty(ridx)
    [cidx,tmp]=find( A(ridx,:)' );
    cidx = reshape(cidx,[3 length(ridx)])';
    cidx1 = cidx(:,1);
    idx1 = sub2ind(size(A),ridx,cidx1);
    cidx2 = cidx(:,2);
    idx2 = sub2ind(size(A),ridx,cidx2);
    cidx3 = cidx(:,3);
    idx3 = sub2ind(size(A),ridx,cidx3);

    % All terms are non-neg
    rsign = (sign(A(idx1)).*xsign(cidx1)>=tol) & (sign(A(idx2)).*xsign(cidx2)>=tol) ...
                & (sign(A(idx3)).*xsign(cidx3)>=tol);
    idx = find(rsign==1);    
    for i1=idx
        xsign(cidx1(i1)) = 0;
        xsign(cidx2(i1)) = 0;
        xsign(cidx3(i1)) = 0;        
    end
    
    % All terms are non-pos
    rsign = (sign(A(idx1)).*xsign(cidx1)<=-tol) & (sign(A(idx2)).*xsign(cidx2)<=-tol) ...
                & (sign(A(idx3)).*xsign(cidx3)<=-tol);
    idx = find(rsign==1);    
    for i1=idx
        xsign(cidx1(i1)) = 0;
        xsign(cidx2(i1)) = 0;
        xsign(cidx3(i1)) = 0;        
    end
end

%--------------------------------------------------------------------
% Local function to update sign of optimization var
%--------------------------------------------------------------------
function xsignNew = LOCALupdate(xsign,xsignUpdate) 
% Find constraints that force ai=0, ai<0 and ai>0
zidx = find( (xsignUpdate==0) | (xsign==-1 & xsignUpdate>0.5) |  (xsign==+1 & xsignUpdate<-0.5));
nidx = find( isnan(xsign) & xsignUpdate<-0.5);
pidx = find( isnan(xsign) & xsignUpdate>0.5 );

% Update xsign
xsignNew = xsign;
xsignNew( zidx ) = 0;
xsignNew( nidx ) = -1;
xsignNew( pidx ) = +1;

 
