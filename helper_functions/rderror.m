function [err, perm, A_align, lambda_align, E] = rderror(varargin)
% Returns the error of the X_est in terms of X_true
% Code inspired by alignment algorithms in
%       https://github.com/NicolasBoumal/HeterogeneousMRA
    
    A_est = varargin{1};
    A_true = varargin{2};
    
    [L, R] = size(A_est);
    
    assert(all([L,R] == size(A_true)),...
                    'A_est and A_true must have identical size');
                
    if nargin>=5
        lambda_est = varargin{3};
        lambda_true = varargin{4};
        n = varargin{5};
    else
        lambda_est = ones(1,R);
        lambda_true = ones(1,R);
        n = varargin{3};
    end
    
    lambda_est = lambda_est.*vecnorm(A_est).^n;
    normA_true = vecnorm(A_true);
    lambda_true = lambda_true.*normA_true.^n;
    A_est = A_est./vecnorm(A_est);
    A_true = A_true./normA_true;
    
    % Calculate distance between all pairs of vectors.
    % Use abs because rank decompositions are unique up to sign flip
    E = lambda_est.^2 + lambda_true'.^2 ...
          -2*lambda_true'.*lambda_est.*(A_true'*A_est).^n;
    
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
   
    % lambda_align is a permutation of lambda_est
    lambda_align = lambda_est(perm);
    normA_true = normA_true(perm);
    
    % Stable error calculation to avoid negative errors due to numerical
    % error
    dot = sum(A_true.*A_align);
    A_ = A_true.*dot-A_align;
    err = sqrt(norm(lambda_true-dot.*lambda_align)^2 + lambda_align *...
            (lambda_true.*sum(A_.*A_true).^n - ...
            lambda_align.*sum(A_.*A_align).^n)');
        
    lambda_align = lambda_align./normA_true.^n;
    A_align = A_align .* normA_true;
    
end