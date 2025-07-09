function [A, B, stat] = spm_21sym(T, R, varargin)
% Decompose symmetric even order tensor using subspace power method
%   ** Usage **
%       X = subspace_power_method(T, L, n, R, opts)
%       [X, lambda] = subspace_power_method(T, L, n, R, opts)
%       [X, lambda, stat] = subspace_power_method(T, L, n, R, opts)
%
%   ** INPUT **
%        T: Tensor of dimension L^n
%        L: Length of tensor (optional, it can be obtained from T if it has
%           the correct shape)
%        n: Tensor order (optional, it can be obtained from T if it has the
%           correct shape)
%        R: Tensor rank (optional, it can be estimated using the 
%           eigenvalues of mat(T))
%     opts: Various SPM options, given as a struct or Parameter/Value pairs
%           Options include
%            ntries: Maximum number of iterations of power method
%           gradtol: Gradient tolerance (the power method finishes if the
%                    norm of the gradient is smaller than this value)
%            eigtol: Tolerance for selecting the rank of T using the
%                    eigenvalues (when R is not provided)
%              ftol: If the function value (in the power method) is less
%                    than this value, then restart x. This is useful when
%                    rank(T)<L and the first guess for x is almost
%                    orthogonal to the span of the a_i. Without this check
%                    the convergence when this happened would be very slow
%          adaptive: Flag indicating if using adaptive shifts (depending on
%                    current function value) or fixed. Defaults to true.
%                    
%
%   ** OUTPUT **
%        X: L x R matrix where the columns are the rank decomposition of T
%           If lambda is also returned the columns of X have norm 1.
%   lambda: Scaling factors (optional). If not returned, X is scaled
%           appropriately.
%     stat: Various statistics of SPM.
%
%   NOTE : Make sure 'helper_functions/' are added to path

% Reference:
% J. Kileel, J. M. Pereira, Subspace power method for symmetric tensor
%                           decomposition
% https://github.com/joaompereira/SPM
% 
% version 1.2 (07/17/2024) - MIT License
    
    timer = tic;
    
    if nargin<2; R = []; end
    
    [m_, m, n] = size(T);

    assert(m_ == m);

    %% Set options here
    try
        opts = option_parser(varargin, {'maxiter', 5000, @(x) x>0},...
                                       { 'ntries', 5, @(x) x>0},...
                                       {'gradtol', 1e-15, @(x) x>0},...
                                       {'ranksel', 1e-4},...
                                       {   'ftol', 1e-2, @(x) x>0});        
    
    catch ME
        if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
            setup
            error("Folders added to path. Please re-run the code.")
        end
        rethrow(ME)
    end

    % Flatten T
    T = reshape(T, m,[]);

    [U, D, V] = svd(T, 'econ');
    D = diag(D);
    
    % Determine tensor rank by the eigenvalues of mat(T)
    if isempty(R)
        typical = mean(abs(D).^2) / mean(abs(D));
        R = sum(abs(D) > opts.ranksel * typical);
    end
    
    D1 = diag(1./D(1:R));

    V = V(:,1:R);
    
    U = U(:,1:R);
    
    A = zeros(m, R);
    B = zeros(n, R);

    lap = toc(timer);
    stat.extracttime = lap;
    stat.powertime = 0;
    stat.deflatetime = 0;
    stat.avgiter = 0;
    stat.nrr = 0;

    for k = R:-1:1
        
        for tries = 1:opts.ntries
          
          % Initialize Xk
          Ak = randn(m,1);
          Ak = Ak/norm(Ak);
          Bk = randn(n,1);
          Bk = Bk/norm(Bk);

          for iter = 1:opts.maxiter

            % Calculate contraction of V with x^(n2-1)
            VAk = reshape(Ak' * reshape(V, m, n * k), n, k);
            Ck = (Bk'*VAk)';

            Bk = VAk*Ck;
            Bk = Bk / norm(Bk);

            Ak_new = reshape(pagemtimes(reshape(V, m, n, k), Bk), m, k) * Ck;

            f = Ak' * Ak_new;

            Ak_new = Ak_new / norm(Ak_new);
            err = norm(Ak - Ak_new);
            Ak = Ak_new;

            if err < opts.gradtol
                % Algorithm converged
                break
            end
            
          end 
          
          stat.avgiter = stat.avgiter + iter;
          
          if 1-f<opts.ftol
             break
          elseif tries==1 || f>f_
              stat.nrr = stat.nrr + 1;
              f_ = f;
              Ak_ = Ak;
              Bk_ = Bk;
          else
              stat.nrr = stat.nrr + 1;
              Ak = Ak_;
              Bk = Bk_;
          end
          
        end
        
        timenow = toc(timer);
        stat.powertime = stat.powertime + timenow - lap;
        lap = timenow;
        
        alphaU = (Ak'*U)';
        alphaV = (reshape(Ak*Bk', 1, []) * V);

        % Solve for lambda
        D1alphaU = D1*alphaU;
        D1alphaV = (alphaV * D1)';
        lambdak = norm(alphaU)*norm(alphaV)/(alphaV*D1alphaU);

        if k > 1
            % Calculate the new matrix D and the new subspace

            % Use Householder reflection to update V and D
            y = (sign(D1alphaU(k))/norm(D1alphaU))*D1alphaU;
            xk = sqrt(1+y(k));
            x = [y(1:k-1)/xk;xk];

            D1 = LHR(D1,x);
            V = RHR(V,x);
            
            y = (sign(D1alphaV(k))/norm(D1alphaV))*D1alphaV;
            xk = sqrt(1+y(k));
            x = [y(1:k-1)/xk;xk];
            
            D1 = RHR(D1,x);
            U = RHR(U,x);

        end
        
        A(:, k) = Ak;
        B(:, k) = lambdak * Bk;
        
        timenow = toc(timer);
        stat.deflatetime = stat.deflatetime + timenow - lap;
        lap = timenow;

    end

    stat.avgiter = stat.avgiter/R;
    stat.totaltime = toc(timer);    
        
end

function A = LHR(A,x)

A = A(1:end-1, :) + x(1:end-1) * (-x'*A);

end

function A =  RHR(A,x)

A = A - (A*x)*x';
A = A(:, 1:end-1);

end