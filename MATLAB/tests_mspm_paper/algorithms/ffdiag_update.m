function W = ffdiag_update(R_tau, ortho)
    % R_tau cell of symmetric matrices
    % Single update for the non-orthogonal FFDIAG algorithm.
    % Set ortho = true to do the orthogonal version.

    % Initialize variables
    n_lags = length(R_tau); % Number of lags
    dim = size(R_tau{1}, 1); % Dimension of the matrices
    Dk = cell(1, n_lags);
    Ek = cell(1, n_lags);

    % Compute Dk and Ek
    for tau = 1:n_lags
        Dk{tau} = diag(diag(R_tau{tau}));
        Ek{tau} = R_tau{tau} - Dk{tau};
    end

    W = zeros(dim, dim);

    if ~ortho
        % Non-orthogonal case
        z = zeros(dim, dim);
        y = zeros(dim, dim);

        for i = 1:dim
            for j = 1:dim
                for tau = 1:n_lags
                    z(i, j) = z(i, j) + Dk{tau}(i, i) * Dk{tau}(j, j);
                    y(i, j) = y(i, j) + Dk{tau}(j, j) * Ek{tau}(i, j);
                end
            end
        end

        % Compute W
        for i = 1:dim
            for j = i+1:dim
                denom = z(j, j) * z(i, i) - z(i, j) * z(i, j);
                W(i, j) = (z(i, j) * y(j, i) - z(i, i) * y(i, j)) / denom;
                W(j, i) = (z(i, j) * y(i, j) - z(j, j) * y(j, i)) / denom;
            end
        end
    else
        % Orthogonal case
        num = zeros(dim, dim);
        den = zeros(dim, dim);

        for i = 1:dim
            for j = i+1:dim
                for tau = 1:n_lags
                    num(i, j) = num(i, j) + Ek{tau}(i, j) * (Dk{tau}(i, i) - Dk{tau}(j, j));
                    den(i, j) = den(i, j) + (Dk{tau}(i, i) - Dk{tau}(j, j))^2;
                end
                if i ~= j
                    W(i, j) = num(i, j) / den(i, j);
                    % W must be skew-symmetric (W = -W^T)
                    W(j, i) = -W(i, j);
                end
            end
        end
    end
end
