function [M] = generate_lowrank_tensor_general(A, W, n)

if nargin<3 || isempty(n)
    n = length(size(W{1}));
    if n==2 && any(size(W{1})==1)
        n = 1;
    end 
end

M = 0;

for i=1:length(A)
    Mi = tucker_product(A{i},W{i},n);
    M = M + Mi;
end