function sos = sosconstr(sos,Type,symexpr)
% SOSCONSTR --- Add a new constraint (equality/inequality)
%    to an SOS program
%
% SOSP = sosconstr(SOSP,TYPE,EXPR)
%
% SOSP is the sum of squares program.
% The new constraint is described by TYPE as follows:
%   TYPE = 'eq'   : constraint is an equality, viz., f(x) = 0
%   TYPE = 'ineq' : constraint is an inequality, viz., f(x) >= 0 (an SOS)
% SYMEXPR is the expression in the left hand side of the constraint, i.e., f(x).
%   SYMEXPR may be of dpvar, pvar, or symbolic format and may be either
%   matrix, vector, or scalar valued. However, matrix and vector-valued
%   inputs are only compatible with 'eq' type constraints.
% 

% NOTE: REVERSIONARY

% NOTE: There is no existing SOSTOOLS function which calls DPsosconstr with
% matrix or vector-valued symexpr when the type is pvar or sym.

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
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log and developer notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 12/24/01 - SP
% 03/01/02 - SP -- New syntax
% 04/06/03 - PJS -- Handle poly objects w/out for-loops
% 7/29/21  - MP -- Added support for dpvar expression types.
% 02/14/22 - DJ -- Add "Type" to getequation input for dpvar case.
% 02/21/22 - DJ -- Add check to convert "double" to appropraite data class.

sos.expr.num = sos.expr.num+1;
expr = sos.expr.num;
sos.expr.type{expr} = Type;


if isfield(sos,'symvartable')
    
    if isdouble(symexpr)    % DJ, 02/21/22
        symexpr = sym(symexpr);
    end
    
    charvartable = converttochar([sos.vartable]);
    isvectorvar = size(symexpr,2)==1;
    for i = 1:size(symexpr,1)
        if isvectorvar
            degcheck = evalin(symengine,char(symexpr(i)));
        else
            degcheck = evalin(symengine,char(symexpr(i,i)));
        end
        degcheck = feval(symengine,'expand',degcheck);
        degcheck = feval(symengine,'collect',degcheck,charvartable);
        degcheckmon = feval(symengine,'poly2list',degcheck,charvartable);
        mondeg = zeros(length(degcheckmon),1);
        for j = 1:length(degcheckmon)
            test = degcheckmon(j);
            mondeg(j) = sum(double(test(2)));
        end
        mondegmax = max(mondeg);
        mondegmin = min(mondeg);
%        if mod(mondegmax,2)~=0||mod(mondegmin,2)~=0
%            error('SOSTOOLS: Leading and trailing degree of SOS expression can not be odd. I suggest you check your polynomial expression.');
%        end
    end
    
    
    [sos.expr.At{expr},sos.expr.b{expr},sos.expr.Z{expr}] = ...
        getequation(char(symexpr),sos.vartable,sos.decvartable,sos.varmat.vartable);
else
%     

if isdouble(symexpr)    % DJ, 02/21/2022
    symexpr = dpvar(symexpr);
elseif isa(symexpr,'polynomial')
    % Convert polynomial to dpvar to avoid issue in sossolve
    symexpr = poly2dpvar(symexpr,sos.decvartable);
    % Not efficient, and not desirable if people want reversionary
    % implementation...
end
if isa(symexpr,'dpvar')    
    symexpr = compress(combine(symexpr,'extended'));    % Get a minimal representation
        [sos.expr.At{expr},sos.expr.b{expr},sos.expr.Z{expr}] = ...
            getequation(symexpr,sos.vartable,sos.decvartable,sos.varmat.vartable,Type);
        
% elseif isa(symexpr,'polynomial')
%         
%         [sos.expr.At{expr},sos.expr.b{expr},sos.expr.Z{expr}] = ...
%             getequation(symexpr,sos.vartable,sos.decvartable,sos.varmat.vartable);
else
    error('symexpr type not recognized')
end

end
