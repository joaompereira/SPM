function [A, lambda, stat] = subspace_power_method(T, L, n, R, options)
% Decompose symmetric even order tensor using subspace power method
%   ** Usage **
%       X = subspace_power_method(T, L, n, R)
%       [X, lambda] = subspace_power_method(T, L, n, R)
%       [X, lambda, stat] = subspace_power_method(T, L, n, R)
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
%   lambda: Scaling factors (optional). If not returned, X is scaled
%           appropriately.
%     stat: Various statistics of SPM.
%
%   NOTE : Make sure 'helper_functions/' are added to path

    if ~exist('L','var') || isempty(L)
        L = size(T,1);
    end
    
    if exist('K','var') && ~isempty(K)
        assert(mod(R,K)==0);
    end

    if ~exist('n','var') || isempty(n)
        n = round(log(numel(T))/log(L));
    end
    
    n2 = n/2;
    d = L^n2;
    
    assert(mod(n,2)==0 && n > 0,'n is not even an even positive integer');

    %% Set options here
    if ~exist('options', 'var') || isempty(options)
        options = struct();
    end
    
    maxiter = setordefault(options,'maxiter', 2000);
     ntries = setordefault(options, 'ntries', 3);
    gradtol = setordefault(options,'gradtol', 1e-14);
     eigtol = setordefault(options, 'eigtol', 1e-6);
       ftol = setordefault(options,   'ftol', 1e-2/sqrt(L));
       
    timer = tic;

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
    
    D1 = diag(1./D(1:R));

    V = findsymscale.*symV(findsym,1:R);
    
    % Pre-allocation of X and lambda
    A = zeros(L,R);
    if nargout>1
        lambda = zeros(1,R);
    end
    
    % C_n from Lemma 4.7
    if n2<=4
      cn = sqrt(2*(n2-1)/n2);
    else
      cn = (2-sqrt(2))*sqrt(n2);
    end
    
    lap = toc(timer);
    stat.extracttime = lap;
    stat.powertime = 0;
    stat.deflatetime = 0;
    stat.avgiter = 0;

    for k = R:-1:1
            
        V_ = reshape(V,[],L*k);

        for tries = 1:ntries
          
          % Initialize Xk
          Ak = randn(L,1);
          Ak = Ak/norm(Ak);
        
          for iter = 1:maxiter

            % Calculate power of Xk
            Apow = Ak;
            for i=2:n2-1
                Apow = reshape(Apow.*Ak',[],1);
            end

            % Calculate contraction of V with x^(n2-1)
            VAk = reshape(Apow'*V_,L,k);
                        
            Ak_new = VAk*(Ak'*VAk)';

            f = Ak_new'*Ak;

            % Determine optimal shift
            % Sometimes due to numerical error f can be greater than 1
            f_ = max(min(f,1),.5);
            clambda = sqrt(f_*(1-f_));
            shift = cn*clambda;

            if f < ftol
                % Xk was not a good initialization
                % Initialize it again at random
                Ak = randn(L,1);
                Ak = Ak/norm(Ak);
            else
                % Shifted power method
                Ak_new = Ak_new + shift*Ak;
                Ak_new = Ak_new/norm(Ak_new);

                if norm(Ak - Ak_new) < gradtol
                    % Algorithm converged
                    Ak = Ak_new;
                    break
                else
                    Ak = Ak_new;
                end
            end
          end 
          
          stat.avgiter = stat.avgiter + iter;
          
          if 1-f<eigtol
             break
          elseif tries==1 || f>f_
              f_ = f;
              Ak_ = Ak;
          else
              Ak = Ak_;
          end
          
        end
        
        timenow = toc(timer);
        stat.powertime = stat.powertime + timenow - lap;
        lap = timenow;

        % Calculate power of Xk
        Apow = Ak;
        for i=2:n2
            Apow = reshape(Apow.*Ak',[],1);
        end

        % Calculate projection of Xpow in subspace
        alpha = (Apow'*V)';

        % Solve for lambda
        D1alpha = D1*alpha;
        lambdak = 1/(alpha'*D1alpha);

        if nargout == 1
            A(:,k) = Ak*lambdak^(1/n);
        else
            A(:,k) = Ak;
            lambda(k) = lambdak;
        end
        
        

        if k > 1
            % Calculate the new matrix D and the new subspace

            
            % Use Householder reflection to update V and D
            y = (sign(D1alpha(1))/norm(D1alpha))*D1alpha;
            xk = sqrt(1+y(k));
            x = [y(1:k-1)/xk;xk];
            
            D1 = RHR(LHR(D1,x),x);
            
            V = RHR(V,x);
            %V(:,end) = V*x;
            %Q = sparse(1:k-1, 1:k-1, 1,k, k-1,2*(k-1));
            %Q(k,:) = -x(1:end-1);
            %V = V*Q;
        end
        
        timenow = toc(timer);
        stat.deflatetime = stat.deflatetime + timenow - lap;
        lap = timenow;

    end
    
    stat.avgiter = stat.avgiter/R;
    stat.totaltime = toc(timer);

end

function A_ =  LHR(A,x)

A_ = A(1:end-1,:) + x(1:end-1) * (-x'*A);

end

function A_ =  RHR(A,x)

A_ = A(:,1:end-1) + (A*(-x))*x(1:end-1)';

end




