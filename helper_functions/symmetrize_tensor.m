function [SymT] = symmetrize_tensor(M, L, n)
% Implementation of the symmetrizing operator

if ~exist('L','var') || isempty(L) || ~exist('n','var') || isempty(n)
    n = length(size(M));
    L = size(M,1);
    assert(all(L==size(M)),'M not a symmetric tensor');
end
    
    if n>1
        sizeM = size(M);
        % This implementation lets symmetrizing K tensors at the same time
        K = numel(M)/L^n;
        M = reshape(M,[L*ones(1,n),K]);
        for k=2:n
            % Recursively symmetrize only k-th first dimensions
            % This takes k operations, for a total operation time of
            % n(n+1)/2
            iterM=M;
            for j=1:k-1
                M = M + permute(iterM,[j+1:k 1:j k+1:n+1]);
            end
        end
    end
    
    SymT=reshape(M,sizeM)./factorial(n);

end
