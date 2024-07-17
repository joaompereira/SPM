function [A, varargout] = subspace_power_method(T, L, n, R, varargin)
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
    
    if nargin<2 || isempty(L); L = size(T,1); end
    if nargin<3 || isempty(n); n = round(log(numel(T))/log(L)); end
    if nargin<4; R = []; end
    
    n2 = ceil(n/2);
    d = L^n2;
    
    assert(n > 2);

    %% Set options here
    try
        opts = option_parser(varargin, {'maxiter', 5000, @(x) x>0},...
                                       { 'ntries', 3, @(x) x>0},...
                                       {'gradtol', 1e-14, @(x) x>0},...
                                       {'ranksel', 1e-4},...
                                       {   'ftol', 1e-2, @(x) x>0},...
                                       {'frestol', 1e-2/sqrt(L), @(x) x>0},...
                                      {'adaptive', true, @(x) x>0});        
    
    catch ME
        if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
            setup
            error("Folders added to path. Please re-run the code.")
        end
        rethrow(ME)
    end


    timer = tic;

    % Flatten T
    T = reshape(T,d,[]);

    %%
    % In this block of code we exploit the fact that each column and row of   
    % mat(T) is a symmetric tensor of order d and therefore has a lot of 
    % repeated entries. Instead of calculating the eigen decomposition of
    % mat(T), which has L^(2d) entries, we calculate an equivalent eigen
    % decomposition of a matrix with roughly (L^d/d!)^2 entries, this way 
    % getting a speed up of approximately (d!)^3.

    if mod(n,2)
        [symind_l, findsym_l, symindscale_l] = symmetric_indices(L, n2);
        [symind_r, findsym_r, symindscale_r] = symmetric_indices(L, n2-1);
        
        symindscale_l = sqrt(symindscale_l);
        symindscale_r = sqrt(symindscale_r);
        findsym_l = reshape(findsym_l,[],1);
        findsym_r = reshape(findsym_r,[],1);
        symind_l = (symind_l-1)*(L.^(0:n2-1)')+1;
        symind_r = (symind_r-1)*(L.^(0:n2-2)')+1;
        
        [symV, D, symU] = svd(symindscale_l.*T(symind_l,symind_r).*...
                              symindscale_r', 'econ');
        D = diag(D);
        
        % Determine tensor rank by the eigenvalues of mat(T)
        if isempty(R)
            typical = mean(abs(D).^2) / mean(abs(D));
            R = sum(abs(D) > opts.ranksel * typical);
        end
        
        D1 = diag(1./D(1:R));

        V = symV(:,1:R)./symindscale_l;
        V = V(findsym_l, :);
        
        U = symU(:,1:R)./symindscale_r;
        U = U(findsym_r, :);
        
        
    else
        [symind, findsym, symindscale] = symmetric_indices(L, n2);

        symindscale = sqrt(symindscale);
        findsym = reshape(findsym,[],1);
        symind = (symind-1)*(L.^(0:n2-1)')+1;
        %%
        % Eigen decomposition
        [symV, D] = eig2(symindscale.*T(symind,symind).*symindscale');

        % Determine tensor rank by the eigenvalues of mat(T)
        if isempty(R)
            typical = mean(abs(D).^2) / mean(abs(D));
            R = sum(abs(D) > opts.ranksel * typical);
        end

        D1 = diag(1./D(1:R));

        V = symV(:,1:R)./symindscale;
        V = V(findsym, :);
    end

    % Pre-allocation of X and lambda
    A = zeros(L,R);
    if nargout>1
        lambda = zeros(1,R);
    end

    % C_n from Lemma 4.7
    cn = sqrt((n2-1)/n2);
    
    lap = toc(timer);
    stat.extracttime = lap;
    stat.powertime = 0;
    stat.deflatetime = 0;
    stat.avgiter = 0;
    stat.nrr = 0;

    for k = R:-1:1
            
        V_ = reshape(V,[],L*k);

        for tries = 1:opts.ntries
          
          % Initialize Xk
          Ak = randn(L,1);
          Ak = Ak/norm(Ak);
        
          for iter = 1:opts.maxiter

            % Calculate power of Xk
            Apow = Ak;
            for i=2:n2-1
                Apow = reshape(Apow*Ak',[],1);
            end

            % Calculate contraction of V with x^(n2-1)
            VAk = reshape(Apow'*V_,L,k);
                        
            Ak_new = VAk*(Ak'*VAk)';

            f = Ak_new'*Ak;

            % Determine optimal shift
            % Sometimes due to numerical error f can be greater than 1
            if ~opts.adaptive
                clambda = 1;
            elseif f > 2/3
                clambda = sqrt(2*f*max(1-f, 0));
            else
                clambda = 1 - f/2;
            end
            shift = clambda * cn;

            if f < opts.frestol
                % Xk was not a good initialization
                % Initialize it again at random
                Ak = randn(L,1);
                Ak = Ak/norm(Ak);
            else
                % Shifted power method
                Ak_new = Ak_new + shift*Ak;
                Ak_new = Ak_new/norm(Ak_new);

                if norm(Ak - Ak_new) < opts.gradtol
                    % Algorithm converged
                    Ak = Ak_new;
                    break
                else
                    Ak = Ak_new;
                end
            end
          end 
          
          stat.avgiter = stat.avgiter + iter;
          
          if 1-f<opts.ftol
             break
          elseif tries==1 || f>f_
              stat.nrr = stat.nrr + 1;
              f_ = f;
              Ak_ = Ak;
          else
              stat.nrr = stat.nrr + 1;
              Ak = Ak_;
          end
          
        end
        
        timenow = toc(timer);
        stat.powertime = stat.powertime + timenow - lap;
        lap = timenow;
        
        if mod(n,2)

            % Calculate power of Xk
            Apow = Ak;
            for i=2:n2-1
                Apow = reshape(Apow*Ak',[],1);
            end

            % Calculate projection of Xpow in subspace
            alphaU = (Apow'*U)';
            
            Apow = reshape(Apow*Ak',[],1);
            
            alphaV = (Apow'*V);

            % Solve for lambda
            D1alphaU = D1*alphaU;
            D1alphaV = (alphaV * D1)';
            lambdak = 1/(alphaV*D1alphaU);

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
            
        else
            
            % Calculate power of Xk
            Apow = Ak;
            for i=2:n2
                Apow = reshape(Apow*Ak',[],1);
            end

            % Calculate projection of Xpow in subspace
            alpha = (Apow'*V)';

            % Solve for lambda
            D1alpha = D1*alpha;
            lambdak = 1/(alpha'*D1alpha);

            if k > 1
                % Calculate the new matrix D and the new subspace

                % Use Householder reflection to update V and D
                y = (sign(D1alpha(k))/norm(D1alpha))*D1alpha;
                xk = sqrt(1+y(k));
                x = [y(1:k-1)/xk;xk];

                D1 = RHR(LHR(D1,x),x);

                V = RHR(V,x);

            end
            
            
        end
        
        if nargout <= 1
            if mod(n, 2) || lambdak >= 0
                A(:,k) = Ak*nthroot(lambdak, n);
            end
            
        else
            A(:,k) = Ak;
            lambda(k) = lambdak;
        end
        
        timenow = toc(timer);
        stat.deflatetime = stat.deflatetime + timenow - lap;
        lap = timenow;

    end
    
    stat.avgiter = stat.avgiter/R;
    stat.totaltime = toc(timer);
    
    if nargout==3; varargout{2} = stat; end
    if nargout>=2; varargout{1} = lambda; end
        
        
end

function A = LHR(A,x)

A = A(1:end-1, :) + x(1:end-1) * (-x'*A);

end

function A =  RHR(A,x)

A = A - (A*x)*x';
A = A(:, 1:end-1);

end