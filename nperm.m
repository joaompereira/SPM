function [n] = nperm(v)
%NMULTPERM Number of permutations of vector with repeated elements
% Given a vector v, with possibly repeated elements, gives the number of
% different permutations of v. The number is given by the multinomial of n
% over the amount of times each unique element of v appears in v.

persistent factorial_store

v = sort(v(:));

m = length(v);

% pre-calculating factorial seemed to improve performance
if isempty(factorial_store) || length(factorial_store)<m
    factorial_store = factorial(1:m);
end

ind = find(v(1:end-1)<v(2:end));

rep = [ind;m]-[0;ind];

if length(rep) == m
    n = factorial_store(m);
else
    n = factorial_store(m)/prod(factorial_store(rep));
end

end
    