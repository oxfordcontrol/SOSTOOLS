function [sos,P] = sossosmatrixvar(sos,ZSym,n,matrixstr,PVoption)
% SOSMATRIXVAR --- Declare a matrix variable P in the sos program, all of
% whose entries are SOS variables
%
% [SOSP,P] = sossosmatrixvar(SOSP,ZSym,n,matrixstr,PVoption)
%
% OUTPUTS:
% SOS: the modified sum of squares program structure.
% P: the new polynomial matrix variable with SOS entries. If SOS contains 
% symbolic variables, the output type of P will be sym. Otherwise, the 
% output type defaults to
% dpvar unless the caller specifies PVoption='pvar' as a fifth input.
%
% INPUTS
% SOS: The sum of squares program structure to be modified.
% ZSym: The vector of monomials to be contained in P. Decision
% variables corresponding to those monomials will be assigned
% automatically by SOSSOSMATRIXVAR.
% matrixstr: a char string with the option 'symmetric' when required
% PVoption (optional): a char string with the PVoption 'pvar' when pvar 
% output type is desired (reversionary option)


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

% JA&GV - 10/01/2013
% MP - 8/6/2021 - added dpvar option. This implementation is inefficient,
% but is little used. Contact developers if a more efficient implementation
% is desired.


% Original Code
if nargin > 3 && ~isempty(matrixstr)
    if strcmp(matrixstr,'symmetric')
        if n(1)==n(2)
            if isfield(sos,'symvartable')
                P = sym(zeros(n(1),n(2)));
                for i = 1:n(1)
                    for j = i:n(1)
                        [sos,var] = sossosvar(sos,ZSym);
                        P(i,j) = var;
                        P(j,i) = var;
                        clear var;
                    end
                end
                return
            elseif nargin==5 && strcmp(PVoption,'pvar')
                % Code for multipoly: PJS 9/9/2013
                P = polynomial(zeros(n(1),n(2)));
                
                for i = 1:n(1)
                    for j = i:n(1)
                        [sos,var] = sossosvar(sos,ZSym,'pvar');
                        P(i,j) = var;
                        P(j,i) = var;
                        clear var;
                    end
                end
                return
            else
                P = dpvar(zeros(n(1),n(2)));
                
                for i = 1:n(1)
                    for j = i:n(1)
                        [sos,var] = sossosvar(sos,ZSym);
                        P(i,j) = var;
                        P(j,i) = var;
                        clear var;
                    end
                end
                return
            end
        else
            disp(['''symmetric''' ' option used, matrix must be square.']);
            P = [];
            return
        end
    else
        disp(['Matrix structure ' matrixstr ' is not defined.' ]);
        P = [];
        return
        
    end
else
    if isfield(sos,'symvartable')
        P = sym(zeros(n(1),n(2)));
        for i = 1:n(1)
            for j = 1:n(2)
                [sos,var] = sossosvar(sos,ZSym);
                P(i,j) = var;
                clear var;
            end
        end
        return
    elseif nargin==5 && strcmp(PVoption,'pvar')
        % Code for multipoly: PJS 9/9/2013
        P = polynomial(zeros(n(1),n(2)));
        for i = 1:n(1)
            for j = 1:n(2)
                [sos,var] = sossosvar(sos,ZSym,'pvar');
                P(i,j) = var;
                clear var;
            end
        end
        return
    else
        P = dpvar(zeros(n(1),n(2)));
        for i = 1:n(1)
            for j = 1:n(2)
                [sos,var] = sossosvar(sos,ZSym);
                P(i,j) = var;
                clear var;
            end
        end
        return
    end
end

