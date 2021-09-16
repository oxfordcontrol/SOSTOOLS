function [y,ndx] = sortrows_integerTable(x)

%sortrows_integerTable: A specialized version of sortrows for 
%sparse, integer-valued, nonnegative data inputs.
%
% Calling [Y,I] = sortrows_integerTable(X) produces the same
% output as [Y,I] = sortrows(X), provided that X is a sparse
% matrix containing non-negative integer values only, but is
% faster.
%
% Author : P. Goulart, 01/07/2013

if(size(x,2) <= 1 || size(x,1) <= 1)
    [y,ndx] = sortrows(x);
    return;
end

%find one plus the maximum entry the whole matrix
d = max(x(:)) + 1;

%protect against empty dataset
if(d == 1)
    y = x;
    ndx = (1:size(x,1))';
    return;
end

%Will assemble a hash-like key using base-d values,
%but can only represent a limited number of columns at ]
%a time this way
K = 1/(2*eps); 
r = floor(log(K)/log(d)) - 1;  %maximum number of base-d digits

%Will need to separate x into 'chunks' of at most 
%r columns.  How many chunks?
hashcols = ceil(size(x,2)/r);

%create an matrix of multipliers for the data,
%such that x*m will produce a (possibly
%multicolumn) hashcode
m = d.^(r:-1:1)';            %multipliers assuming 1 column hashcode
m = kron(speye(hashcols),m); %block diagonal repitition of m 
m = m(1:size(x,2),:);        %drops unneeded values for the last chunk

%sort chunk-wise right to left, forming
%only one column at a time.  It would seem
%faster to do a sortrows on x*m, but this
%is too memory intensive and turns out to
%be slower anyway

ndx = (1:size(x,1))';
for k = size(m,2):-1:1
    %the next hashed block from the right
    v = x*m(:,k);
    [~,ind] = sort(v(ndx,1));
    ndx = ndx(ind,1);
end

y = x(ndx,:);


