
function [A,B,filler] = ffdiag_nonortho(T, eps, max_iter)
    % Blind source separation using FFDIAG. This version does not require that
    % the estimated mixing matrix be orthogonal.
    %
    % INPUT:
    % ------
    % R_tau : cell array
    %     Cell array of covariance matrices.
    %
    % eps : double, optional
    %     Convergence criterion for matrix updates (default: 1.0e-10).
    %
    % max_iter : int, optional
    %     Maximum number of iterations/updates (default: 100).
    %
    % OUTPUT:
    % -------
    % invV : matrix
    %     Inverse of the estimated diagonalizing matrix.
    %
    % V : matrix
    %     Estimated diagonalizing matrix.

    if nargin < 2
        eps = 1.0e-10;
    end
    if nargin < 3
        max_iter = 100;
    end

    k = size(T,3);
    R_tau = cell(1, k);
    % Populate the cell array with slices of the 3D array
    for i = 1:k
        R_tau{i} = T(:, :, i);
    end
    R_tau_old = R_tau;

    dim = size(R_tau{1}, 1); % Dimension of the matrices
    n_lags = length(R_tau); % Number of lags
    W = zeros(dim, dim);
    V = eye(dim);
    C = R_tau;
    niter = 0;
    theta = 0.9;
    iter_eps = 1.0;

    while iter_eps > eps && niter < max_iter
        niter = niter + 1;
        Vn1 = V;

        % Update C
        for tau = 1:n_lags
            C{tau} = (eye(dim) + W) * C{tau} * (eye(dim) + W)';
        end

        % Update term
        W = ffdiag_update(C, false);
        if norm(W, 'fro') > theta
            W = (W * theta) / norm(W, 'fro');
        end

        % Update V
        V = (eye(dim) + W) * V;

        % Compute delta
        delta = 0;
        for i = 1:dim
            for j = 1:dim
                if i ~= j
                    delta = delta + (V(i, j) - Vn1(i, j))^2;
                end
            end
        end
        iter_eps = delta / (dim * (dim - 1));
    end

    invV = inv(V);
    A = invV./vecnorm(invV);
    B = zeros(k,size(T,1));
    for i = 1:k
        B(i,:)= diag(pinv(A)*R_tau_old{i}*pinv(A)');
    end
    filler = [niter,iter_eps];
end