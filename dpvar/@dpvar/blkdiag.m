function Fdiag = blkdiag(varargin)
% F = blkdiag(varargs) generates a block diagonal dpvar object using
% varargs as the diag blocks
% 
% INPUTS:
%   varargs : dpvar or poly object. Atleast one object should be dpvar

% 
% OUTPUTS:
% Fdiag: a dpvar object, Fdiag = [varargs(1)                        ]
%                                              varargs(2)                                                 
%                                                         varargs(n)]
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, MP, SS - 09/27/2021

if nargin==1
    Fdiag=varargin{1};
else % sequentially creating a block diagonal matrix
    A = varargin{1}; B = varargin{2};
    
    % ensure classes are same type
    if ~isa(A,'dpvar')
        A = dpvar(A);
    end
    if ~isa(B,'dpvar')
        B = dpvar(B);
    end
    
    mdim = A.matdim+B.matdim;
    
    % First, ensure the left and right multiplication monomials are the
    % same
    [Anew,Bnew] = common_basis(A,B);
    
    Cnew = blkdiag(Anew.C,Bnew.C);
    Fdiag = dpvar(Cnew, Anew.degmat, Anew.varname, Anew.dvarname, mdim); 
    if nargin>2 % repeat when there are more than two block elements
        Fdiag = blkdiag(Fdiag,varargin{3:end});
    end
end
end