function [A, B, stat] = MSPM_21sym(T, R, varargin)
% MSPM_21sym - Decompose a 3D tensor (m x m x n) symmetric in the first two modes
%              using the Multi Subspace Power Method (MultiSPM).
%
%   This function calculates the CP decomposition of a 3D tensor that is symmetric
%   in the first two modes, using the Multi Subspace Power Method.
%
%   ** Usage **
%       [A, B, stat] = MSPM_21sym(T, R, varargin)
%
%   ** INPUT **
%        T: Tensor of size m x m x n, symmetric in the first two modes.
%        R: Tensor rank (optional, can be estimated if not provided).
%     varargin: Optional parameters provided as name-value pairs:
%         - 'maxiter': Maximum number of iterations for the power method
%                      (default: 5000).
%         - 'ntries': Maximum number of attempts for the power method
%                    (default: 5).
%         - 'gradtol': Gradient tolerance for convergence (default: 1e-15).
%         - 'ranksel': Tolerance for rank selection using singular values.
%                      You may provide a function that takes the singular values
%                      as input and decides the rank, or set it to 'return_sv'.
%                      (default: 1e-4).
%         - 'ftol': Function value tolerance for restarting the power method
%                  (default: 1e-2).
%
%   ** OUTPUT **
%        A: Factor matrix corresponding to the first two modes (size m x R).
%        B: Factor matrix corresponding to the third mode (size n x R).
%     stat: Struct containing various statistics of the decomposition process:
%           - extracttime: Time spent on extracting tensor properties.
%           - powertime: Time spent on the power method.
%           - deflatetime: Time spent on deflation.
%           - avgiter: Average number of iterations per rank.
%           - nrr: Number of restarts during the power method.
%           - totaltime: Total runtime of the decomposition.
%
%   ** NOTE **
%   - This function is specialized for tensors of size m x m x n that are symmetric
%     in the first two modes. For general tensors, use the standard MultiSPM function.
%   - Ensure that the required helper functions are added to the MATLAB path.
%
%   ** Reference **
%   K. Wang, J. M. Pereira, J. Kileel, A. Seigal, "Multi-subspace power method
%   for decomposing all tensors",

%   https://github.com/joaompereira/SPM
%
%   ** Version **
%   - Version 1.0 (10/15/2025) - MIT License
    
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
    
    % Determine tensor rank by the singular values of mat(T)
    if isempty(R)
        R = rank_selector(D, opts.ranksel);
        if isempty(R)
            A = D;
            return
        end
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
                        
            Bk = VAk*(Bk'*VAk)';
            Bk = Bk / norm(Bk);

            VBk = reshape(pagemtimes(reshape(V, m, n, k), Bk), m, k);

            Ak_new = VBk*(Ak'*VBk)';

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

function r = rank_selector(S, rank_sel)
% Select rank by looking at singular values
    if isa(rank_sel,'function_handle')
        r = rank_sel(S);
    elseif isscalar(rank_sel)
        typical = sum(S.^2) / sum(abs(S));
        r = sum(abs(S) > rank_sel * typical);
    elseif rank_sel == "return_sv"
        r = [];
    else
        error('Rank selector option not implemented yet')
    end

end

function A = LHR(A,x)

A = A(1:end-1, :) + x(1:end-1) * (-x'*A);

end

function A =  RHR(A,x)

A = A - (A*x)*x';
A = A(:, 1:end-1);

end