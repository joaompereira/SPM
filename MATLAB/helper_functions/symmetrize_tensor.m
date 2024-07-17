function [M] = symmetrize_tensor(M, L, n)
% Implementation of the symmetrizing operator

if nargin<3 || isempty(L) || isempty(n)
    n = length(size(M));
    L = size(M,1);
    assert(all(L==size(M)),'M not a symmetric tensor');
end
    
    if n==1
        return
    end
        
    sizeM = size(M);
    % This implementation lets symmetrizing K tensors at the same time
    M = reshape(M,L^n,[]);
    K = size(M, 2);

    if issparse(M)        
        [i, j, v] = find(M);
        
        [cellind{1:n}] = ind2sub(L*ones(1,n),i);
        n_ind = horzcat(cellind{:});   

        sort_i = 1 + sort(n_ind-1, 2, "descend") * L.^(0:n-1)';

        M = sparse(sort_i,j, v, L^n,K);
        
        for k=2:n
            [i, j, v] = find(M);

            [cellind{1:n}] = ind2sub(L*ones(1,n),i);
            n_ind = horzcat(cellind{:});

            indk = repmat(i, 1, k);
            
            for p=2:k
                perm = [p:k 1:p-1 k+1:n]'-1;
                indk(:, p) = 1+ (n_ind-1) * L.^perm;
            end

            indk = reshape(indk, [], 1);
            M = sparse(indk, repmat(j, k, 1), repmat(v, k, 1) / k);

        end        

    else
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
        M = reshape(M,sizeM)./factorial(n);
    end
end

