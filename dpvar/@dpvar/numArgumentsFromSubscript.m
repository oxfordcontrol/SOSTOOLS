function n = numArgumentsFromSubscript(obj,s,indexingContext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This method is automatically called when one of the following statements
% are executed:
% a) obj.s — Number of elements referenced in a statement (indexingContext = Statement) 
% b) (obj.s) — Number of elements returned in an expression (indexingContext = Expression) 
% c) [obj.s] = rhs — Number of values assigned with a comma-separated list (indexingContext = Assignment)
% The function retunrs the expected number of outputs n for each of these
% three types of calls.
% Since we overloaded the function 'numel', we need this function to avoid
% breaking 'subsref' and 'subsasgn'.
% For now, this function just returns 1, expecting 1 output in all cases,
% though something more sophisticated could be implemented as well.

   n=1;
end