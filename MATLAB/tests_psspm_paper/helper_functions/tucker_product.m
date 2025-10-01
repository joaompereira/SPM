function T = tucker_product(M, T, n)
% Returns symmetric multilinear multiplication of T by M

if nargin<3 || isempty(n)
    n = length(size(T));
    if n==2 && any(size(T)==1)
        n = 1;
    end 
end

if n==1
    % Special case n=1
    T = M*T;
else
    [L, K] = size(M);
    
    % If T is a collection of tensors, this algorithm calculates the 
    % symmetric tucker product of all tensors by M, just need that the
    % the last dimension is the number of tensors
    L2 = numel(T)/K^n;
    
    % Replace M by transpose
    M = M.';

    for k = 1:n
        % In order to multiply by each dimension we reshape T and multiply
        % by tranpose of M
        T = reshape(T,K,[]);
        T = T.'*M;        
    end
    
    if L2>1
        T = reshape(T,L2,[]);
        T = T.';
    end
    
    T = reshape(T,[L*ones(1,n), L2]);
end

end

