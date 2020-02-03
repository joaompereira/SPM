function [S_aligned, O, cost] = align_subspaces(S, Sref)
% Input: x and xref both of size NxK
% Output: x_aligned of size NxK: contains the same data as x, with columns
% permuted and circularly shifted (individually) to match xref as closely
% as possible (in some sense defined by the code.)

    assert(all(size(S) == size(Sref)), 'S and Sref must have identical size');
    
    [U, D, V] = svd(Sref'*S);
    
    O = V * U';
    S_aligned = S * O;
    if nargout > 2
        cost = norm(diag(D)-1);
    end

end
