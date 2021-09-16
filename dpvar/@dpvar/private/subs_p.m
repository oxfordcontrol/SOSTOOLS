function b = subs_p(a,old,new)
% b = subs_p(a,old,new) substitutes the variables "old" in polynomial "a"
% by values or alternative variables "new" and returns the resulting
% polynomial b
%
% DESCRIPTION 
%  Symbolic Substitution.
%  This function differs from the command "subs" in the symbolic toolbox in
%  three ways.  First, the output of this function is always a polynomial.
%  Use the command "double" to depromote a constant polynomial to a
%  double. Second, the input "Old" cannot be an arbitrary polynomial
%  expression.  Third, the order of "Old" and "New" cannot be switched.
% 
% INPUTS:
% a: polynomial
% old: polynomial variable or a cell array of polynomial variables
% new: polynomial variable or a cell array of polynomial variables
% 
% OUTPUTS:
% b: polynomial object with variables old assigned a value new
%
% SYNTAX 
%   B = subs(A);
%       -Replaces all variables in A with values in the workspace.
%   B = subs(A,Old,New);
%       -Replaces Old with New in the symbolic expression A.
%       -If Old and New are cell arrays of the same size, each element 
%       of Old is replaced by the corresponding element of New. 
%       -If A and Old are scalars and New is an array or cell array, 
%       the scalars are expanded to produce an array result.
%   B = subs(A,New);
%       -Equivalent to B=subs(A,A.varname{1},New);
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
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% 3/25/2003: PJS  Initial Coding  
% 6/20/2010: MMP  Matrix Reshaping Correction
% 5/31/2013: MMP complete rewrite of polynomial substitution to improve
% scalability. partial rewrite of numerical substitution to improve
% efficiency, eliminate internal errors. replaced strmatch with strcmp.
% Also added comments.

% Promote a to polynomial 
  a = polynomial(a);

% Convert old to a cell array of strings and new to a cell 
% array of doubles and/or polynomials
if nargin == 1
  
  new = cell(0);
  nva = length(a.varname);
  cnt = 1;
  for i1 = 1:nva
    
    if evalin('caller', ['isa(' a.varname{i1} ',''double'')'] );
      old{cnt} = a.varname{i1};
      new{cnt} = evalin('caller',a.varname{i1});
      cnt = cnt +1;   % counts number of 'variables' in a of type 'double' or 'polynomial'
    end
    
    if evalin('caller', ['isa(' a.varname{i1} ',''polynomial'')'] );
      temp = evalin('caller',a.varname{i1});    % retrieve polynomial representation of variable
      if temp~=polynomial(1,1,a.varname(i1),[1 1])
 	old{cnt} = a.varname{i1};    % Converting old to string format
 	new{cnt} = temp;      % Converting new to polynomial format
 	cnt = cnt+1;      % counts number of 'variables' in a of type 'double' or 'polynoial'
      end
    end
    
  end

  if length(new)==0     % nothing to substitute, so polynomial is unchanged
    b = a;
    return;
  end
  
elseif nargin == 2    % nothing to change to, so polynomial is unchanged 
  
  b=subs(a,a.varname{1},old);
  return;
            
elseif nargin == 3

  if ~iscell(old)
    old = {old};% cellular version of old
  end  
  szo = size(old); % number of old variables
  for i1 = 1:szo(1)
    for i2 = 1:szo(2)
      if isa(old{i1,i2},'polynomial');
	old{i1,i2} = old{i1,i2}.varname{1};  % converting old to string format
      elseif ~ischar(old{i1,i2})
	error(['Old should be a polynomial variable or a '  ...
	       'cell array of polynomial variables']);
      end
    end
  end    

  if ~iscell(new)
    new = {new};
  end
  szn = size(new);
  for i1 = 1:szn(1)
    for i2 = 1:szn(2)
      if ~isa(new{i1,i2},'polynomial') & ...
	    ~isa(new{i1,i2},'double')  % ensuring that new is a set of doubles or polynomials
	error(['New should be a polynomial or double or a ' ... 
	      'cell array of polynomials and/or doubles.']);
      end
    end
  end    
  
end

% Some quick error checking
sza = size(a);
szn = size(new{1});
if any(cellfun('size',new,1)~=szn(1)) | ...
      any(cellfun('size',new,2)~=szn(2)) 
  error('All entries of new must have the same dimensions');
end
if any(szn~=1) & any(sza~=1)
  error('Either A or New must be 1x1');
end

% Make New a single array if A and Old are 1x1
%b=a;% start with b=a  % Lets not start with b=a, we'll do this later
szb=size(a);
szo = size(old);
if all(szo==[1 1]) & all(szb==[1 1]) & ...
      all(cellfun('isclass',new,'double')) % detect case of numeric substitution

  new = {cell2mat(new)};   % convert to double
  szn = size(new{1});
  
elseif all(szo==[1 1]) & all(szb==[1 1])                % detect case with single substitution and polynomial is
                % scalar

  tempnew = polynomial;
  for i1 = 1:size(new,1);
    temppoly = polynomial;
    for i2 = 1:size(new,2);
      temppoly = [temppoly new{i1,i2}];   % creates a polynomial matrix version of new 
    end
    tempnew = [tempnew; temppoly];
  end  
  new = {tempnew};
  szn = size(new{1});
  
end

  
old = old(:);
new = new(:);
[nrb,ncb]=size(a);
if all(cellfun('isclass',new,'double')) % Numeric Substitution
  ntb = size(a.degmat,1);
  if all(szb==[1 1])
    bmatdim = szn;
    tempsz = szn(1)*szn(2);
    bcoef = sprepmat(a.coefficient,1,tempsz);    
  else
    tempsz = szb(1)*szb(2);   % number of coefficients for each monomial
    bcoef = reshape(a.coefficient,ntb,tempsz);
%    [nrb,ncb]=size(a);
    bmatdim = [nrb*ncb 1];   
 %    bmatdim=a.matdim;
  end
%  b=a;
  bdegmat=a.degmat;bvarname=a.varname;
  for i1 = 1:length(old)  
%    idx = strmatch(old{i1},bvarname); 
    idx = strcmp(old{i1},bvarname); 
    if idx==0
      error(['Undefined variable ' old{i1}]);
    end
    tempdeg = bdegmat(:,idx);

    bdegmat(:,idx) = [];
    bvarname(idx)=[];  
    
    tempdeg = sprepmat(tempdeg,1,tempsz);       % repeats tempdeg once for each coefficient (columns)
    % for each coefficient in coeff, lists the associated degree for that
    % monomial.
    if all(szb==[1 1])
      tempnew = sprepmat(new{i1}(:)',ntb,1);
    else      
      tempnew = sprepmat(new{i1}(:)',ntb,tempsz);
                % repeats new value for each term (cols) and then for each coefficeint (rows)
                % for each term in coeff, repeats the value of new.
    end
    bcoef = bcoef.*(tempnew.^tempdeg); 
  end

  if isempty(bdegmat)
      bcoef=sum(bcoef);
      b=polynomial(bcoef);
  else
    b=polynomial(bcoef,bdegmat,bvarname,bmatdim);
  end

  
  
else % Symbolic Substitution - MMP 5/31/2013
% The default for pvar.coefficent is that the ith row number corresponds to the
% monomial indicated by the ith row of degmat. For a given row i, the jth column of degmat refers
% to the power of the jth variable in the ith monomial. The ith row of
% pvar.coefficient is the vector form of the matrix coefficient for the
% ith monomial.


% To substitute a single variable xi in to g(x), for each power of xi, we will
% construct the polynomial which multiplies that power g(x) = a + b xi+cxi^2+.... Then we will
% substitute in for each h_k=(xi->new)^k and construct g'(x) = a + b*h_1 + c*h_2 + ....

%  nvb = size(b.degmat,2);   %unused
b=a;
%[nrb,ncb]=size(b)
  for i1 = 1:length(old);
    bmatdim = [nrb*ncb 1];   
    ntb = size(b.degmat,1);
    bcoefficient = reshape(b.coefficient,ntb,nrb*ncb); % This should do nothing unless something is wrong....
    bdegmat = b.degmat;     % At present, there is no apparant ordering to the monomials in degmat
    bvarname = b.varname;
    temp = strcmp(old{i1},b.varname);
 %   temp = strmatch(old{i1},b.varname); 
    if temp==0
%       error(['Undefined variable ' old{i1}]);
      disp(['Undefined variable' old{i1}]);
      return %%Check this modification
    else
      idx = temp;        % gives the column of degmat associated with variable old(i1)
    end
    
  
  % Now we create a polynomial without the variable old(i1)
  bdegmat(:,idx) = [];
  bvarname(idx)=[];  
    % Need to handle the case where the coefficent is 0 or a constant (e.g. no other variables)
    
 % now for each power, we create a new polynomial which multiplies that power
    temp_degs = b.degmat(:,idx);
    p=polynomial(0); 
    for id=1:(max(temp_degs)+1)
          temp_idx = find(temp_degs==(id-1)); % these are the terms in degmat with a power of id-1 in variable old(i1)
            % Need to handle the cases where the coefficent is 0 or a constant (e.g. no other variables)
          if ~isempty(temp_idx)   % What happens if there are no monomials with power id-1? Then skip this term
              bdegmat2 = bdegmat(temp_idx,:);
          bcoefficient2 = bcoefficient(temp_idx,:);
              if isempty(bvarname) % if there are no other variables, then the coefficient is scalar
                  btemp2 = polynomial(bcoefficient2);
              else
                  btemp2 = polynomial(bcoefficient2,bdegmat2,bvarname,bmatdim);
              end
            % Now create new term of new^(id-1)
          pow = power(new{i1},id-1);
          p=p+pow*btemp2; % multiply the power by the associated term
          end
    end
%    b=p;
    b = combine(p);
    end
end


  

  
% Combine down
b = combine(b);
if isempty(a.degmat)
    b=a;
else
    b=reshape(b,nrb,ncb);     % 6/20/2010: MMP  Matrix Reshaping Correction
end












