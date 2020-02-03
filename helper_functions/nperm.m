function [n] = nperm(v)
%NMULTPERM Number of permutations of vector with repeated elements
% Given a vector v, with possibly repeated elements, gives the number of
% different permutations of v. The number is given by the multinomial of n
% over the amount of times each unique element of v appears in v.

persistent mm

if isempty(mm)
    mm = memoize(@multinomial);
end

v = sort(v(:));

m = length(v);

ind = [0;find(v(1:end-1)<v(2:end));m];

rep = sort(diff(ind));
n = mm(rep);

end
   
function n = multinomial(v)

    csv = cumsum(v);
    n = prod(arrayfun(@nchoosek,csv(2:end),v(2:end)));

end
