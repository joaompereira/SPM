function [x_aligned, perm, O, mincost] = align_to_reference_multiPCA(x, xref)
% Input: x and xref both of size NxK
% Output: x_aligned of size NxK: contains the same data as x, with columns
% permuted and circularly shifted (individually) to match xref as closely
% as possible (in some sense defined by the code.)

    assert(all(size(x) == size(xref)), 'x and xref must have identical size');
    
    [L, K, R] = size(x);
    
    E = zeros(R);
    
    for i = 1:R
        for j = 1:R
            [~, S, ~] = svd(x(:,:,i)'*xref(:,:,j));
            E(j,i) = norm(diag(S)-1,1);
        end
    end
    [perm, mincost] = munkres(E);
    
    O = zeros(K,K,R);
    
    for k=1:R
        [x_aligned(:,:,k), O(:,:,k)] = ...
                align_subspaces(x(:,:,perm(k)), xref(:,:,k));
    end

end
