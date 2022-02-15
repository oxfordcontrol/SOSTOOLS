classdef (InferiorClasses={?polynomial})dpvar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This defines the class dpvar, which is a variant of the multipoly class
    % where decision variables (for which dependence is always affine) are
    % isolated from true independent variables.
    %
    % Polynomials of the dpvar class are convertible to polynomials of
    % the multipoly class. Alternatively, polynomials of the multipoly class
    % can be converted to dpvars by specifying the names of the decision
    % variables in the polynomial.
    %
    % obj = dpvar(C, degmat, varname, dvarname, matdim)
    %
    % A polynomial of the dpvar class is formed as
    % D(d;p) = (Z_1(d)^T \otimes I_m)C (Z_2(p)\otimes I_n)
    % where
    %
    % D.matdim: an 1x2 array with entries [m,n] is the dimension of the
    % matrix. If not specified, dpvar will default to a scalar - [1 1]
    %
    % D.dvarname: A vector of names of the nd decision variables (d). Z_1(d) is then
    % simply the vector of all decision variables with 1 in the first position
    %
    % D.varname: This format is the same as in multipoly and consists of the
    % length np vector of all independent variables (p) in the polynomial and the ordering of
    % these variables names is synchronized with degment
    %
    % D.degmat: This is likewise the same as in multipoly and allows us to
    % formulate the vector of monomials Z_2(p) as
    %
    % Z(p) = varname{1}.^degmat(:,1).*...*varname{p}.^degmat(:,p)
    %
    % D.C: This is the matrix of elements used to define the
    % polynomial D(d;p) as above. So that each element adds a term to
    %
    % D(d;p)_{i,j} = C_{k+i*m,l+j*n} * d_{k}*varname{1}.^degmat(l,1).*...*varname{p}.^degmat(l,np)
    %
    % NOTES:
    % For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
    % S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS 2021b - dpvar
    %
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % If you modify this code, document all changes carefully and include date
    % authorship, and a brief description of modifications
    %
    % Initial coding DJ, MP, SS - 07/12/2021
    % 12/15/21 - DJ -- Add single input single output case (y=dpvar('x'))
    % 02/14/22 - DJ -- Add option to input cell of strings (y=dpvar({'x1','x2'}))
    
    
    
    properties (Access=public)
        dvarname = {}; % this is cell array with decision varnames as strings
        varname = {}; % this is cell array with polynomial varnames as strings
        C = sparse([]); % this is the coeff that multiplies monomials,
        degmat = sparse([]); % this specifies the monomial set
        matdim = [0 0];
    end
    properties (SetAccess=private)
        chkval;
       % private properties can only be set in class methods and class submembers
    end
    methods
        function obj = dpvar(varargin)
            if nargin == 0
                % Return dpvar object obj(p;d) with default values
                return
            elseif iscellstr(varargin) % if input is a cell array of strings create dpvars with those names
                if nargout==0
                    for i=1:nargin
                        Cf = sparse([0;1]);
                        dmat = zeros(1,0);
                        vname = {};
                        dvname = {varargin{i}};
                        matdim = [1,1];
                        assignin('caller', varargin{i}, dpvar(Cf,dmat,vname,dvname,matdim));
                        clear obj;
                    end
                elseif nargin==1 && nargout==1                   
                    Cf = sparse([0;1]);
                    dmat = zeros(1,0);
                    vname = {};
                    dvname = {varargin{1}};
                    matdim = [1,1];
                    obj = dpvar(Cf,dmat,vname,dvname,matdim);
                else
                    error('For cell char inputs, number of outputs should be at most 1');
                end
            elseif nargin==1 % single input is either dpvar, double or polynomial
                if isa(varargin{1},'dpvar')
                    % If input is already a dpvar
                    obj = varargin{1};
                elseif isa(varargin{1},'polynomial')
                    % If input is a poly, convert to dpvar
                    obj = poly2dpvar(varargin{1});
                elseif isa(varargin{1},'double') && isreal(varargin{1}) && ismatrix(varargin{1})
                    % If input is a real double, then convert to a poly
                    szc = size(varargin{1});
                    obj.C = varargin{1}; 
                    obj.degmat = sparse(1,0);
                    obj.varname  = {};
                    obj.dvarname = {};
                    obj.matdim = szc;
                    obj.chkval = 1; % this skips error check
                elseif iscellstr(varargin{1})
                    [nd1,nd2] = size(varargin{1});
                    nd = nd1*nd2;
                    if nargout==1   
                        Cf = sparse(reshape([repmat([0;1;zeros(nd,1)],nd-1,1);0;1],nd1*(nd+1),nd2));
                        dmat = zeros(1,0);
                        vname = {};
                        dvname = varargin{1}(:);
                        matdim = [nd1,nd2];
                        obj = dpvar(Cf,dmat,vname,dvname,matdim);
                    else
                        error('For cellstr input, number of outputs should be 1')
                    end
                else
                    error(['For single input, argument must be a dpvar, '...
                        'real double, or a cell of strings']);
                end
            elseif nargin<4
                errstr1 = 'Invalid number of inputs for the "dpvar" command.';
                error([errstr1]);
            elseif nargin<=6 % more than one input means fields of the dpvar object
                if nargin==4   % adding an option to skip the matdim for scalars
                    varargin{5} = [1 1];
                    varargin{6} = 0;
                elseif nargin==5
                    varargin{6} = 0;
                end
                % convert row cells to columns cells and ensure inputs are
                % row/column vectors
                if size(varargin{3},1)==1
                    varargin{3} = varargin{3}';
                end
                if size(varargin{4},1)==1 % convert row cells to columns cells
                    varargin{4} = varargin{4}';
                end

                obj.varname = varargin{3};
                obj.dvarname = varargin{4};
                obj.C = varargin{1};
                obj.degmat = varargin{2};
                obj.matdim = varargin{5}; % why (:)' -DJ????
                obj.chkval = varargin{6}; % this is optional which specifies if error check needs to be skipped
                
                % check if inputs have correct format and dimensions
                if ~obj.chkval
                    [logval, errmsg] = dpvarconstructmsg(obj);
                    if logval~=0
                        error(errmsg);
                    else
                        obj.chkval=1;
                    end
                end
            else 
                error('Too many inputs for "dpvar" class object');
            end
        end
    end
    methods (Access=private)
        function [logval, errmsg] = dpvarconstructmsg(obj)
            sizeC = size(obj.C);
            % error checking,
            % logval 0, no problem
            % logval 1, incorrect data type
            % logval 2, incorrect dimensions
            
            % check datatype
            logval=0;
            errmsg = '';
            if ~isa(obj.C,'double')
                logval=1;
                errmsg = "Coefficient matrix should be matrix/sparse matrix double";
                return
            end
            if ~isa(obj.degmat, 'double')
                logval=1;
                errmsg = "Coefficient matrix should be matrix/sparse matrix double";
                return
            end
            if any(mod(obj.matdim,1)~=0)
                logval=1;
                errmsg = "matdim should be an integer array of length 2";
                return
            end
            if ~iscellstr(obj.varname)
                logval=1;
                errmsg = "varnames should be cell array of chars";
                return
            end
            if ~iscellstr(obj.dvarname)
                logval=1;
                errmsg = "dvarnames should be cell array of chars";
                return
            end
            
            
            % check dimensions
            if length(obj.varname)~=numel(obj.varname)
               logval = 2;
               errmsg = "varnames should be a 'row/col' cell array";
               return
            end
            if length(obj.dvarname)~=numel(obj.dvarname)
               logval = 2;
               errmsg = "dvarnames should be a 'row/col' cell array";
               return
            end
            if numel(obj.matdim)~=2
                logval=2;
                errmsg = "matdim should be an integer array of length 2";
                return
            end
            if (sizeC(1)~=obj.matdim(1)*(length(obj.dvarname)+1)) ||...
                    (sizeC(2)~=obj.matdim(2)*(size(obj.degmat,1)))
                logval = 2;
                errmsg = "Incorrect C dimension. dpvar coefficient matrix should have dimension"+...
                    " rows*(n_dvars+1) x cols*(rowsize_degmat). Don't forget that the first row of C has no associated decision variable.";
                return
            end
            if size(obj.degmat,2)~=length(obj.varname)
                logval=2;
                errmsg = "degmat column dimension should match number of variable names";
                return
            end
        end
    end
end
