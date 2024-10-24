function [T, symvec] = generate_lowrank_tensor(lambda, varargin)

factors = varargin;

if ~isempty(lambda) && size(lambda,1) > 1
    factors = [{lambda},factors];
    lambda = [];
end

nfacts = length(factors);

if size(factors{nfacts}, 1) == 1
    nsyms = factors{nfacts};
    nfacts = nfacts - 1;
    factors = factors(1:nfacts);
    assert(size(nsyms, 2) == nfacts)
else
    nsyms = ones(1, nfacts);
end   

rank = size(factors{1},2);

if isempty(lambda)
    lambda = ones(1,rank);
end

udims = cellfun(@(M) size(M,1), factors);

symvec = repelem(1:nfacts, nsyms);
order = length(symvec);

o2 = floor(order/2);
T = khatri_rao_product(lambda, factors{symvec(1:o2)}) * ...
        khatri_rao_product(factors{symvec(o2+1:end)})';
T = reshape(T, udims(symvec));

end

function krp = khatri_rao_product(varargin)
    if nargin == 1
        krp = varargin{1};
    else
        m2 = ceil(length(varargin) / 2);
        A = khatri_rao_product(varargin{1:m2});
        B = khatri_rao_product(varargin{m2+1:end});
        r = size(A, 2);
        krp = reshape(A, [], 1, r) .* reshape(B, 1, [], r);
        krp = reshape(krp, [], r);
    end
end
    
