function [A, W, stat] = subspace_power_method_general(T, L, n, R, K, varargin)
% Decompose symmetric even order tensor using generalized subspace power method
%   ** Usage **
%       [A, W, stat] = subspace_power_method_general(T, L, n, R, K, options)
%
%   ** INPUT **
%        T: Tensor of dimension L^n
%        L: Length of tensor (optional, it can be obtained from T if it has
%           the correct shape)
%        n: Tensor order (optional, it can be obtained from T if it has the
%           correct shape)
%        R: Tensor rank (optional, it can be estimated using the 
%           eigenvalues of mat(T))
%        K: Subspace dimension (optional, it can be calculated using the
%           singular values of the jacobian)
%  options: Struct with various SPM options, such as
%            ntries: Maximum number of iterations of power method
%           gradtol: Gradient tolerance (the power method finishes if the
%                    norm of the gradient is smaller than this value)
%            eigtol: Tolerance for selecting the rank of T using the
%                    eigenvalues (when R is not provided)
%           eigtol2: Tolerance for selecting the size of each subspace
%                    using the eigenvalues (when K is not provided)
%              ftol: If the function value (in the power method) is less
%                    than this value, then restart x. This is useful when
%                    rank(T)<L and the first guess for x is almost
%                    orthogonal to the span of the a_i. When this happens 
%                    if x is not restarted the power method will take
%                    longer to convergence
%
%   ** OUTPUT **
%        A: A cell containing L x K_i matrices with the generalized rank
%           decomposition subspaces.
%   Lambda: A cell containing the corresponding scaling factors 
%     stat: Various statistics of SPM. (optional) 
%
%   NOTE : Make sure 'helper_functions/' are added to path

% Reference:
% J. Kileel, J. M. Pereira, Subspace power method for symmetric tensor
%                           decomposition and generalized PCA
% https://github.com/joaompereira/SPM
% 
% version 1.1 (06/07/2021) - MIT License

    if nargin<2 || isempty(L); L = size(T,1); end
    if nargin<3 || isempty(n); n = round(log(numel(T))/log(L)); end
    if nargin<4; R = []; end
    if nargin<5; K = []; end
       
    n2 = n/2;
    d = L^n2;
    
    assert(mod(n,2)==0 && n > 0,'n is not even an even positive integer');
    
    if ~isempty(K)
        assert(~isempty(R) && mod(R,nchoosek(K+n2-1,n2))==0);
    end
    
    %% Set options here   
    try
        [maxiter, ntries, gradtol, eigtol, eigtol2, ftol] = ...
            option_parser(varargin, {'maxiter',  2000, @(x) x>0},...
                                    { 'ntries',     3, @(x) x>0},...
                                    {'gradtol', 1e-14, @(x) x>0},...
                                    { 'eigtol',  1e-8, @(x) x>0},...
                                    {'eigtol2',  1e-3, @(x) x>0},...
                                    {   'ftol',  1e-2/sqrt(L), @(x) x>0});
    catch ME
        if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
            warning('Did you add ''helper_functions'' to the MATLAB path?');
        end
        rethrow(ME)
    end
    
    timer = tic;

    %% Flatten T
    T = reshape(T,d,d);

    %%
    % In this block of code we exploit the fact that each column and row of   
    % mat(T) is a symmetric tensor of order d and therefore has a lot of 
    % repeated entries. Instead of calculating the eigen decomposition of
    % mat(T), which has L^(2d) entries, we calculate an equivalent eigen
    % decomposition of a matrix with roughly (L^d/d!)^2 entries, this way 
    % getting a speed up of approximately (d!)^3.
    
    [symind, findsym, symindscale] = symmetric_indices(L, n2);

    symindscale = sqrt(symindscale);
    findsym = reshape(findsym,[],1);
    findsymscale = 1./symindscale(findsym);
    symind = (symind-1)*(L.^(0:n2-1)')+1;
    %%
    % Eigen decomposition
    [symV, D] = eig2(symindscale.*T(symind,symind).*symindscale');

    % Determine tensor rank by the eigenvalues of mat(T)
    if isempty(R)
        R = sum(abs(D) > eigtol);
    end

    D1 = diag(1./D(1:R));

    V = findsymscale.*symV(findsym,1:R);
    

    % C_n from Lemma 4.7
    if n2<=4
      cn = sqrt(2*(n2-1)/n2);
    else
      cn = (2-sqrt(2))*sqrt(n2);
    end
    
    genrank = 1;
    
    lap = toc(timer);
    stat.extracttime = lap;
    stat.powertime = 0;
    stat.deflatetime = 0;
    stat.avgiter = 0;
    stat.nrr = 0;

    while R > 0
      
        V_ = reshape(V,[],L*R);

        for tries = 1:ntries
          
          % Initialize Ak
          Ak = randn(L,1);
          Ak = Ak/norm(Ak);
        
          for iter = 1:maxiter

            % Calculate power of Ak
            Apow = Ak;
            for i=2:n2-1
                Apow = reshape(Apow.*Ak',[],1);
            end

            % Calculate contraction of V with x^(n2-1)
            VAk = reshape(Apow'*V_,L,R);
                        
            Ak_new = VAk*(Ak'*VAk)';

            f = Ak_new'*Ak;

            % Determine optimal shift
            % Sometimes due to numerical error f can be greater than 1
            f_ = max(min(f,1),.5);
            clambda = sqrt(f_*(1-f_));
            shift = cn*clambda;

            if f < ftol
                % Ak was not a good initialization
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

        % Calculate power of Ak
        Apow = Ak;
        for i=2:n2-1
            Apow = reshape(Apow.*Ak',[],1);
        end

        VAk = reshape(Apow'*V_,L,R);

        % Using formula (D.3)
        [Ker, KerD] = eig2(eye(L)/n2+((n2-1)/n2)*(Ak*Ak') - VAk*VAk');

        if ~isempty(K)
            KerR = K;
        else
            KerR = max(sum(KerD < eigtol2),1);
        end

        Ker = Ker(:,L-KerR+1:L);

        [symind, findsym, scale] = symmetric_indices(KerR, n2);
        nR = size(symind,1);
        scale = sqrt(scale);
        Kerpow = Ker(:,symind(:,1)) .* scale';
        for j=2:n2
            Kerpow = khatri_rao_product(Kerpow, Ker(:,symind(:,j)), nR);
        end
        Kerpow = symmetrize_tensor(Kerpow, L, n2);
        nU = (Kerpow'*V)';
        

        % Kerpow has all tensor products of columns of Ker
        %Kerpow = Ker;
        %for i=2:n2
        %    Kerpow = kron(Kerpow, Ker);
        %end

        % Obtain a orthonormal basis for the projection of the columns of
        % Kerpow in V
        %[nU,nS,nV] = svd(V'*Kerpow);

        %nR = nchoosek(KerR+n2-1,n2);

        %nU = nU(:,1:nR);
        %nV = nV(:,1:nR);

        % W is given by Lemma ___
        D1nU = D1*nU;
        WV1 = nU'*D1nU;
        WV1 = (WV1 + WV1')/2;
        %WV = inv(WV1);

        A{genrank} = Ker;
        findsym = findsym(:);
        findsymscale = reshape(scale(findsym),[],1);
        iWV1 = inv(WV1);
        W_flat = iWV1(findsym,findsym)./(findsymscale.*findsymscale');
        W{genrank} = reshape(W_flat, KerR*ones(1,n));

        R = R - nR;
        
        % DEFLATION
        if R > 0
            genrank = genrank + 1;
            
            % See hh_preprocess (below) for more info
            % on Householder reflection implementation
            [X, Xt] = hh_preprocess(D1nU);
            Xt = Xt(:, nR+1:nR+R);
            
            % Using formulas in Appendix B          
            D1 = D1(:, nR+1:nR+R) - (D1 * X) * Xt;
            D1 = D1(nR+1:nR+R, :) - Xt' * (X' * D1);
            
            V = V(:, nR+1:nR+R) - (V * X) * Xt;

        end
        
        timenow = toc(timer);
        stat.deflatetime = stat.deflatetime + timenow - lap;
        lap = timenow;

    end
    
    stat.avgiter = stat.avgiter/genrank;
    stat.totaltime = toc(timer);


end

function AB = khatri_rao_product(A, B, k)

if nargin< 3 || isempty(k)
    k = size(A, 2);
end

AB = reshape(A, [], 1, k) .* reshape(B, 1, [], k);
AB = reshape(AB, [], k);

end

function [X, Xt] = hh_preprocess(Y)

    opts.UT = true;
    
    [M, N] = size(Y);
    % The lower bottom of qr gives the coefficients
    % for the householder reflection
    X = tril(qr(Y), -1);
    % Add ones in the diagonal of X
    X(1+(0:N-1)*(M+1)) = 1;
    
    T = X'*X;
    % Divide the diagonal of T by 2
    T(1+(0:N-1)*(N+1)) = diag(T)/2;
    Xt = linsolve(T, X', opts);%T \ X';
    
    % Now to build the orthogonal matrix one can use
    % the formula Q = eye(M) - X * T^(-1) * X
    %               = eye(M) - X * Xt
    % However if multiplying and M > N,
    % then is faster to use the formulas
    %      AQ = A - (A * X) * Xt
    %  or  QA = A - X * (Xt * A)
  
end