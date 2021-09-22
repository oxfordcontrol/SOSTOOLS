function p = sosgetsol(sos,V,digit)
% 
% SOSGETSOL --- Get the solution from a solved SOS program 
%
% SOL = sosgetsol(SOSP,VAR,DIGIT) 
%
% SOL is the solution from a solved sum of squares program SOSP,
% obtained through substituting all the decision variables
% in VAR by the numerical values which are the solutions to
% the corresponding semidefinite program. 
%
% The third argument DIGIT (optional) will determine the 
% accuracy of SOL in terms of the number of digits. Default 
% value is 5.
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
% 12/27/01 - SP
% 02/21/02 - SP -- Symbolic polynomial
% 03/01/02 - SP -- New syntax
% 03/15/02 - SP -- Fast
% 19/06/14 - GV, extends getsol for matrix variable under pvar
% 06/09/13 - MP -- faster implementation for matrix-valued polynomials
% 06/25/20 - Sachin -- fixed bug related constant polynomials

if nargin == 2
    digit = 5;   % Default
end

if isfield(sos,'symvartable')
    
    p = mysymsubs(V,sos.symdecvartable,sos.solinfo.RRx(1:length(sos.symdecvartable)),digit);

else
    if isa(V,'dpvar')
        V=dpvar2poly(V);
    end
    
    [~,idxdecvar1,idxdecvar2] = intersect(V.varname,sos.decvartable);
    idxvar = setdiff(1:length(V.varname),idxdecvar1);
    coeffs = [V.degmat(:,idxdecvar1) 1-sum(V.degmat(:,idxdecvar1),2)]*[sos.solinfo.RRx(idxdecvar2);1];
                                                %MMP 6/9/2013 updated to
                                                %allow for matrix-valued
                                                %polynomials and to fix
                                                %problem with terms with no
                                                %decision variables
    coeffs = V.coefficient.*repmat(coeffs,1,size(V.coefficient,2));         % 01/31/02
    varname = V.varname(idxvar);
    if isempty(idxvar)
        degmat = [];
    else
        degmat = V.degmat(:,idxvar);
    end
 %   p = set(V,'varname',varname,'degmat',degmat,'coefficient',coeffs);
%     degmat 
%     coeffs
%     size(V)
    if isempty(degmat)
      coeffs=sum(coeffs,1); % modified by sachin - 6/25/2020 original version sum(coeffs)
      p=polynomial(coeffs);
      p=reshape(p,size(V));
  else
    p=polynomial(coeffs,degmat,varname,size(V));
   end
    p=combine(p);
end
