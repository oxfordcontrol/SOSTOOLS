function [sos,P] = sospolymatrixvar(sos,ZSym,n,matrixstr,PVoption)
% SOSPOLYMATRIXVAR --- Declare a polynomial matrix variable P in the sos program
% SOS of
%
% [SOSP,P] = sospolymatrixvar(SOSP,ZSym,n,matrixstr,output)
%
% OUTPUTS:
% SOS: the modified sum of squares program structure.
% P: the new polynomial matrix variable. If SOS contains symbolic variables, the
% output type of P will be sym. Otherwise, the output type defaults to
% dpvar unless the caller specifies output='pvar' as a fifth input.
%
% INPUTS
% SOS: The sum of squares program structure to be modified.
% ZSym: The vector of monomials to be contained in P. Decision
% variables corresponding to those monomials will be assigned
% automatically by SOSPOLYVAR.
% n: The desired dimension of the output matrix P: n(1) x n(2)
% matrixstr (optional): a char string with the option 'symmetric' when required
% PVoption (optional): a char string with the option 'pvar' when pvar output type is desired (reversionary option)

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
%
% JA&GV - 6/13/2013
%
% -MP 6/27/2021: updated to dpvar output. Disabled use of symbolic
% variables. 
% -DJ, 12/10/2021: Adjusted call to sospolyvar when using 'pvar' option


% Original Code
if nargin > 3 && ~isempty(matrixstr)
    if strcmp(matrixstr,'symmetric')
        if n(1)==n(2)
            if isfield(sos,'symvartable')
                P = sym(zeros(n(1),n(2)));
                for i = 1:n(1)
                    for j = i:n(1)
                        [sos,var] = sospolyvar(sos,ZSym);
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
                        [sos,var] = sospolyvar(sos,ZSym,[],'pvar'); % DJ, 12/10/1997
                        P(i,j) = var;
                        P(j,i) = var;
                        clear var;
                    end
                end
                return
            else
                % quadvar has no obvious way to construct the symmetric matrix we
                % want here, so using the older version which does not call
                % quadvar
                [sos,P]=sospolyvar_mat(sos,ZSym,n,'symmetric');
                % Code for multipoly: PJS 9/9/2013
                %                P = dpvar(zeros(n(1),n(2)));
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
                [sos,var] = sospolyvar(sos,ZSym);
                P(i,j) = var;
                clear var;
            end
        end
        
    elseif nargin==5 && strcmp(output,'pvar')
        P = polynomial(zeros(n(1),n(2)));
        for i = 1:n(1)
            for j = 1:n(2)
                [sos,var] = sospolyvar(sos,ZSym,[],'pvar'); % DJ, 12/10/1997
                P(i,j) = var;
                clear var;
            end
        end
        
    else
        
        [sos,P]=sosquadvar(sos,1,ZSym,n(1),n(2));
        
        %         if nargin > 4 & wscoeff == 'wscoeff'
        %             var = sos.var.num;
        %             for i = sos.var.idx{var}:sos.var.idx{var+1}-1
        %                 pvar(['coeff_',int2str(i-sos.var.idx{1}+1)]);
        %                 assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],eval(['coeff_',int2str(i-sos.var.idx{1}+1)]));
        %             end
        %         end
    end
    
end
