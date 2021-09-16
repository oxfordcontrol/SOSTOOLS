function [X,Y] = MosekSol2SedumiSol(K,res)
% MosekSol2SedumiSol --- Convert Mosek output to Sedumi output format.
%
% [X,Y] = MosekSol2SedumiSol(K,res)
%
% K: is the cone structure from the Sedumi input format.
% res: is the output structure from Mosek.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
% 
% 04/24/14 - MP - initial coding

x = res.sol.itr.xx; % the LP solution
y = res.sol.itr.y;  % the dual LP variables

xx=[];
yy=[];
idx=0;
for i=1:length(K.s)
    I=find(tril(ones(K.s(i))));
    nX1=zeros(K.s(i));
    nX1(I)=res.sol.itr.barx((idx+1):(idx+length(I)));
    nS1=zeros(K.s(i));
    nS1(I)=res.sol.itr.bars((idx+1):(idx+length(I)));
    nX1=nX1+nX1'-diag(diag(nX1));
    nS1=nS1+nS1'-diag(diag(nS1));
    idx=idx+length(I);
    xx=[xx;vec(nX1)];
    yy=[yy;vec(nS1)];
end
X=[x;xx];
Y=[y;yy];
%INFO=res.info;

