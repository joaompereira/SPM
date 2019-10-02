function [X, lambda] = subspace_power_method(T, L, n, R, options)
% Decompose symmetric even order tensor using subspace power method
%   ** Usage **
%       X = even_tensor_decomposition(T, L, n, R)
%       [X, lambda] = even_tensor_decomposition(T, L, n, R)
%
%   ** INPUT **
%        T: Tensor of dimension L^n
%        L: Length of tensor (optional, it can be obtained from T if it has
%           the correct shape)
%        n: Tensor order (optional, it can be obtained from T if it has the
%           correct shape)
%        R: Tensor rank (optional, it can be estimated using the 
%           eigenvalues of mat(T))
%     opts: Struct with various SPM options, such as
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
%                    
%
%   ** OUTPUT **
%        X: L x R matrix where the columns are the rank decomposition of T
%           If lambda is also returned the columns of X have norm 1.
%   lambda: Scaling factors (optional) 
%

    %% Set options here
    if ~exist('opts', 'var') || isempty(opts)
        options = struct();
    end
    
    maxtries = setdefault(options, 'maxtries', 1000);
     gradtol = setdefault(options,  'gradtol', 1e-14);
      eigtol = setdefault(options,  ' eigtol', 1e-3);
        ftol = setdefault(options,     'ftol', 1e-2/sqrt(L));
    
    if ~exist('L','var') || isempty(L)
        L = size(T,1);
    end

    if ~exist('n','var') || isempty(n)
        n = round(log(numel(T))/log(L));
    end

    assert(mod(n,2)==0 && n > 0,'n is not even an even positive integer');

    n2 = n/2;
    d = L^n2;

    % Flatten T
    T = reshape(T,d,d);

    %%
    % In this block of code we exploit the fact that each column and row of   
    % mat(T) is a symmetric tensor of order d and therefore has a lot of 
    % repeated entries. Instead of calculating the eigen decomposition of
    % mat(T), which has L^(2d) entries, we calculate an equivalent eigen
    % decomposition of a matrix with roughly (L^d/d!)^2 entries, this way 
    % getting a speed up of approximately (d!)^3.

    % Number of rows and columns of new matrix
    dsym = nchoosek(L+n2-1,n2);

    % Different indices without permutations
    symind = nchoosek(1:L+n2-1,n2)-(0:n2-1);
    symind = symind(:,end:-1:1);
    symindscale = zeros(dsym,1);

    % Map from all indices to indices without permutations
    findsym = zeros(ones(1,n2)*L);
    findsymscale = zeros(ones(1,n2)*L);

    for i=1:dsym
        % Scaling to keep have the same dot product in the new space
        scale = sqrt(nperm(symind(i,:)));
        symindscale(i) = scale;
        scale = 1/scale;
        for j = perms(1:n2)'
            commalist = num2cell(symind(i,j));
            findsym(commalist{:}) = i;
            findsymscale(commalist{:}) = scale;
        end
    end

    findsym = reshape(findsym,[],1);
    findsymscale = reshape(findsymscale,[],1);
    symind = (symind-1)*(L.^(0:n2-1)')+1;
    %%
    % Eigen decomposition
    [symV, D] = eig2(symindscale.*T(symind,symind).*symindscale');

    % Determine tensor rank by the eigenvalues of mat(T)
    if ~exist('R','var') || isempty(R)
        R = sum(abs(D) > eigtol);
    end

    D = D(1:R);

    V = findsymscale.*symV(findsym,1:R);

    % Pre-allocation of X and lambda
    X = zeros(L,R);
    if nargout>1
        lambda = zeros(R,1);
    end

    for k = R:-1:1

    % Initialize Xk
    Xk = randn(L,1);
    Xk = Xk/norm(Xk);

    for tries = 1:maxtries

        % Calculate power of Xk
        Xpow = Xk;
        for i=2:n2
            Xpow = kron(Xpow, Xk);
        end

        % Calculate projection of Xk in subspace
        alpha = (Xpow'*V)';
        PXk = reshape(V*alpha,L*ones(1,n2));

        %Apply the power method on projection PXk
        Xk_new = PXk;
        for i=2:n2
            Xk_new = shiftdim(sum(Xk_new.*Xk));
        end

        f = Xk_new'*Xk;

        % Determine optimal shift
        % Sometimes due to numerical error f can be greater than 1
        lambdak = max(min(f,1),.5);
        shift = (n2-1)*sqrt(lambdak*(1-lambdak));

        if f < ftol
            % Xk was not a good initialization
            % Initialize it again at random
            Xk = randn(L,1);
            Xk = Xk/norm(Xk);
        else
            % Shifted power method
            Xk_new = Xk_new + shift*Xk;
            Xk_new = Xk_new/norm(Xk_new);

            if norm(Xk - Xk_new) < gradtol
                % Algorithm converged
                Xk = Xk_new;
                break
            else
                Xk = Xk_new;
            end
        end
    end

    % Calculate power of Xk
    Xpow = Xk;
    for i=2:n2
        Xpow = kron(Xpow, Xk);
    end

    % Calculate projection of Xpow in subspace
    alpha = (Xpow'*V)';

    % Solve for lambda
    lambdak = (sum(alpha.^2.*D.^-1))^-1;

    if nargout == 1
        X(:,k) = Xk*lambdak^(1/n);
    else
        X(:,k) = Xk;
        lambda(k) = lambdak;
    end

    if k > 1
        % Calculate the new matrix D and the new subspace
        [V2,D] = eig2(diag(D) - alpha*alpha'*lambdak);

        % Update eigenvalues and orthonormal basis for the subspace
        D = D(1:k-1);
        V = V * V2(:,1:k-1);    
    end

    end

end


  

