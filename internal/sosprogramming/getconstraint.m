function [A,ZZ] = getconstraint(Z)
% GETCONSTRAINT --- Find constraint for sum of squares decomposition.
%
% [A,ZZ] = getconstraint(Z)
%
% Z is a monomial vector description.
% This function computes the constraint matrix A and the polynomial
% vector ZZ, such that if q satisfies
%
%    A*q = F'
%
% Then
%
%    Z'*Q*Z = F*ZZ = q'*A'*ZZ      (where q=Q(:) is the vector form of Q)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log and developer notes
% 12/01/01 - SP
% Update: 12/12/01 -- Use sparse form. SP.
% Update: 24/12/01 -- Extra ZZ output. SP.
% 12/17/22 - DJ: Completely new version, use "uniquerows_integerTable".

% Zero/One monomial case (won't occur a lot...)
nZ = size(Z,1);
if nZ<=1
    ZZ = 2*Z;
    A = speye(nZ);
    return
end

toggle = 1;
if toggle==1
    % % Letting q = Q(:) and Z = [Z1; ...; Zn], we can represent:
    % %
    % % Z'*Q*Z = [Z1*Z; Z2*Z; ...; Zn*Z]' * Q(:)
    % %        = ZZ_full' * q 
    % %        = q' * ZZ_full
    % % 
    % % However, Zi*Zj = Zj*Zi. Therefore, we build a reduced set of
    % % monomials, ZZ_red, including only degrees Zi*Zj for i>=j:
    % % ZZ_red = [ [Z1*Z1; Z1*Z2; ...; Z1*Zn] ]
    % %          [ [       Z2*Z2; ...; Z2*Zn] ]
    % %          [                      :     ]
    % %          [                    [Zn*Zn] ]

    % Establish logical indices associated to i>=j.
    full_idcs = (1:nZ^2)';
    log_retain = tril(true(nZ,nZ));
    log_vec = log_retain(:);

    % Construct ZL := [ [Z1; Z1; ...; Z1]; [Z2; ...; Z2]; ...; Zn] ];
    ZL = repelem(Z,(nZ:-1:1)',1);

    % Construct ZR := [ [Z1; Z2; ...; Zn]; [Z2; ...; Zn]; ...; Zn] ];
    ZR_idcs = tril(toeplitz((1:nZ)) + (0:nZ-1));
    ZR = Z(ZR_idcs(log_vec),:);

    % Construct ZZ_red = ZL.*ZR;
    ZZ_red = ZL + ZR;
    nZ_red = size(ZZ_red,1);

    % Determine which row in ZZ_full corresponds to which row in ZZ_red;
    red_idcs = (1:nZ_red)';
    new_idcs = zeros(nZ,nZ);
    new_idcs(log_vec) = red_idcs;
    new_idcs = new_idcs + triu(new_idcs',1);
    new_idcs = new_idcs(:);     % ZZ_full(k,:) = ZZ_red(new_idcs(k),:);

    % Build permutation matrix R s.t.   ZZ_full = R' * ZZ_red;
    R = sparse(new_idcs,full_idcs,1,nZ_red,nZ^2,nZ^2);


    % % Finally, discard non-unique monomials from ZZ_red, finding ZZ and A
    % % such that:
    % % Z'*Q*Z = [Z*Z1; ...; Z*Zn]' * Q(:)
    % %        = ZZ_full' * q 
    % %        = q' * ZZ_full
    % %        = q' * R' * ZZ_red = q' * A' * ZZ
    [A,ZZ] = uniquerows_integerTable(R,ZZ_red,'transpose'); % A'*ZZ = R'*ZZ_red

else
% Original implementation (much slower...)
    %
    % We write Z'*Q*Z as
    %
    % Z'*Q*Z = [Z1; ...; Zn]' * Q * Z
    %        = (Z1*e1)'*Q*Z  + ... + (Zn*en)'*Q*Z
    %        = e1'*Q*(Z*Z1)  + ... + en'*Q*(Z*Zn)
    %        = e1'*Q*(R1*ZZ) + ... + en'*Q*(Rn*ZZ)
    %        = (e1'*Q*R1     + ... + en'*Q*Rn) * ZZ
    %
    % where ej is jth standard basis vector.
    % We compute a common basis ZZ, corresponding to unique monomials in
    % [Z1*Z; ...; Zn*Z], and determine permutation matrices Rj s.t. 
    % Z*Zj = Rj*ZZ.
    R = cell(nZ,1);  % {R1; ...; Rn}
    
    % Initial guess ZZ = Z1*Z;
    ZZ = Z + sprepmat(Z(1,:),nZ,1); 
    R{1} = speye(nZ);    % Z1*Z = I*ZZ;
    
    % Construct R2 through Rn
    for j = 2:nZ
        % Determine Z_tmp = Zj*Z;
        Z_tmp = Z + sprepmat(Z(j,:),nZ,1);
        % Merge ZZ with Z_tmp, returning R1, R2 s.t. 
        %   ZZ_old = R_old*ZZ_new,     Ztemp = R_tmp*ZZ_new;
        [R_old,R_tmp,ZZ] = findcommonZ(ZZ,Z_tmp);     
        % Update Rk for k=1,...j-1:   Zk*Z = Rk*ZZ_old = (Rk*R_old)*ZZ_new;
        R(1:j-1) = cellfun(@(Rk) Rk*R_old, R(1:j-1),'UniformOutput',false);
        % Store new Rj s.t. Zj*Z = Rj*ZZ_new
        R{j} = R_tmp;
    end
    
    % Construct the constraint equations:
    % Z'*Q*Z = q'*A'*ZZ
    Q = sparse([],[],[],nZ,nZ,1);
    A = sparse(size(ZZ,1),nZ^2);         % 12/12/01
    % Loop over all elements of vector q = Q(:)
    for i = 1:(nZ)^2
        Q(i) = 1; 
        % Each row in Q corresponds to a particular permutation matrix Rj.
        % (ei*A')*ZZ = ... + (ej*Q*Rj)' * ZZ 
        [j,~] = find(Q);                    % 12/12/01 
        A(:,i) = A(:,i) + R{j}' * Q(j,:)';
        Q(i) = 0;
    end

end