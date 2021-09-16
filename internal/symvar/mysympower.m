function ZTemp = mysympower(vartable,Z)
% MYSYMPOWER --- Fast symbolic power.
%
% ZSYM = mysympower(VARTABLE,ZNUM)
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
% 03/01/02 - SP

ZTemp = '[';
for j = 1:size(Z,1)
    idx = find(Z(j,:));
    charexpr = [];
    for i = idx
        if Z(j,i)>1
            charexpr = [charexpr,char(vartable(i)),'^',int2str(full(Z(j,i))),'*'];
        else
            charexpr = [charexpr,char(vartable(i)),'*'];    
        end;
    end;
    if isempty(idx)
        charexpr = '1*';
    end;
    charexpr = [charexpr(1:end-1),';'];
    ZTemp = [ZTemp,charexpr];
end;    
ZTemp = [ZTemp(1:end-1),']'];
ZTemp = str2sym(ZTemp);