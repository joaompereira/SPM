function [err, A_align, perm, E] = rank_decomposition_error(A_est, A_true)
% Returns the error of the X_est in terms of X_true
% Code inspired by alignment algorithms in
%       https://github.com/NicolasBoumal/HeterogeneousMRA

    assert(all(size(A_est) == size(A_true)),...
                    'A_est and A_true must have identical size');
    
    R = size(A_est, 2);
    
    % Calculate distance between all pairs of vectors.
    % Use abs because rank decompositions are unique up to sign flip
    E = vecnorm(A_est).^2+vecnorm(A_true)'.^2-2*abs(A_true'*A_est);
    
    % Use an out-sourced implementation of Hungarian algorithm to find
    % right permutation of indices that minimize the distances above
    perm = munkres(E);
    
    % A_align is a permutation of the columns of A_est
    A_align = A_est(:, perm);
    
    for k = 1:R
        % Since the rank decomposition are unique up to sign flip, make
        % sure A_align and A_true are not the symmetric of each other
        if A_align(:, k)'*A_true(:, k)<0
           A_align(:, k) = -A_align(:, k);
        end
    end
    
    % Frobenius norm of the entries
    err = norm(A_align(:)-A_true(:));
    

end