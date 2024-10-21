function [T] = generate_lowrank_tensor(lambda, varargin)

factors = varargin;

if ~isempty(lambda) || size(lambda,1) > 1
    factors = [{lambda},factors];
    lambda = [];
end

order = length(factors);
rank = size(factors{1},2);

if isempty(lambda)
    lambda = ones(1,rank);
end

T = reshape(factors{1}.*lambda, [], 1, rank);
dims = size(factors{1},1);

for i = 2:order-1
    dims = [dims, size(factors{i},1)];
    T = reshape(T .* reshape(factors{i}, 1, [], rank), [], 1, rank);
end

dims = [dims, size(factors{order},1)];
T = reshape(T, [], rank) * factors{order}';
T = reshape(T, dims);

end
    
