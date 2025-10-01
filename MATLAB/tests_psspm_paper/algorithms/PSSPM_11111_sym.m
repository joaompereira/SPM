function [varargout] = PSSPM_11111_sym(T, varargin)
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
%                           decomposition and generalized PCA
% https://github.com/joaompereira/SPM
% 
% version 1.1 (06/07/2021) - MIT License

    import helper_functions.*
        
    %% Set options here
    iP = inputParser;
    addOptional(iP,'rank', [], @(x) isempty(x) || x>0);
    addParameter(iP, 'maxiter', 10000, @(x) x>0);
    addParameter(iP,  'ntries', 5, @(x) x>0);
    addParameter(iP, 'gradtol', 1e-15, @(x) x>0);
    addParameter(iP, 'ranksel', 1e-4);
    addParameter(iP,    'ftol', 1e-4, @(x) x>0);
    addParameter(iP,   'flats', []);
    addParameter(iP, 'svals_only', false);
    parse(iP, varargin{:});

    opts = iP.Results;
    r = opts.rank;
                                
    timer = tic;

    dims = size(T);
    order = length(dims);

    flats = opts.flats;
    if isempty(flats)
        [flats, max_rank] = find_best_flatpair(dims);
    else
        if iscell(flats)
            cell_flats = flats;
            flats = false(2, order);
            flats(1, cell_flats{1}) = true;
            flats(2, cell_flats{2}) = true;
        end
        for i=2:-1:1
            left_dims = dims(flats(i,:));
            left_rank = prod(left_dims) - sum(left_dims);
            right_rank = prod(dims(~flats(i,:)));
            flat_ranks(i) = min(left_rank, right_rank);
        end

        max_rank = max(flat_ranks);
        assert(max_rank>1);
    end

    assert(isempty(r) || r<= max_rank);

    dim_order = [find(flats(1, :) & flats(2, :)),...
                 find(flats(1, :) & ~flats(2, :)),...
                 find(~flats(1, :) & flats(2, :)),...
                 find(~flats(1, :) & ~flats(2, :))];
    T_old = T;
    T = permute(T, dim_order);

    %i_dim_order(dim_order) = 1:order;

    %size(T)
    %dims(dim_order)
    %dims(i_dim_order)

    dims = dims(dim_order);

    m1 = sum(flats(1, :));
    m2 = sum(flats(2, :));
    mc = sum(flats(1, :) & flats(2, :));
    mu = sum(flats(1, :) | flats(2, :));

    matT_1 = reshape(T, prod(dims(1:m1)), []);
    
    [U1, S, V1] = svd(matT_1, 'econ');
    
    S = diag(S);

    % Determine tensor rank by the singular values of mat(T)
    if isempty(r)
        r = rank_selector(S, opts.rank_sel);
        if isempty(r)
            varargout = {S};
            return
        end
    end
    
    U1 = U1(:, 1:r);
    U1_copy = U1;
    S = S(1:r);
    C = diag(1./ S);
    C_copy = C;
    V1 = V1(:, 1:r);

    T = permute(T, [1:mc, m1+1:mu, mc+1:m1, mu+1:order]);
    
    matT_2 = reshape(T, prod(dims([1:mc, m1+1:mu])), []);
    [U2,S2,V2] = svd(matT_2, 'econ');
    U2 = U2(:, 1:r);
    V2 = V2(:, 1:r);
    
    % Pre-allocation of X and lambda
    factors = arrayfun(@(dim) zeros(dim, r), dims,...
            	       'UniformOutput',false);
    
    lambda = zeros(1,r);
    
    lap = toc(timer);
    stat.extracttime = lap;
    stat.powertime = 0;
    stat.deflatetime = 0;
    stat.avgiter = 0;
    stat.nrr = 0;

    Aks = arrayfun(@(dim) zeros(dim, 1), dims,...
            	       'UniformOutput',false);

    ntries = opts.ntries;
    maxiter = opts.maxiter;

    for k = r:-1:1
          
        [Aks(1:m1), f] = power_method(U1, Aks(1:m1), k);

        %log10(1+1e-14-f)

        
        %% Find right-side vectors by using the second flattening
        if mc>0 && m2 > mc
            Akpow = tensor_product(Aks(1:mc));
            U_half = Akpow' * reshape(U2, length(Akpow), []);
            U_half = reshape(U_half, [], r);
            [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);%, 1+1e-12
            %[v, ~] = svds(U_half, 1);
            %[v, ~] = svds(reshape(U_half, [], r), 1, 1+1e-12);
            [Aks(m1+1:mu), f_c] = power_method(v, Aks(m1+1:mu), 1);

            %log10(1+1e-14-f_c)

            if order > mu
                Akpow = tensor_product(Aks(m1+1:mu));
                U_half = Akpow' * reshape(V1, length(Akpow), []);
                U_half = reshape(U_half, [], r);
                [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);
                [Aks(mu+1:order), f_c] = power_method(v, Aks(mu+1:order), 1);

                %log10(1+1e-14-f_c)
            end
        else
            Akpow = tensor_product(Aks(mc+1:m1));
            U_half = Akpow' * reshape(V2, length(Akpow), []);
            U_half = reshape(U_half, [], r);
            [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);
            [Aks(mu+1:order), f_c] = power_method(v, Aks(mu+1:order), 1);            

            if m2 > mc
                Akpow = tensor_product(Aks(mu+1:order));
                U_half = pagemtimes(reshape(V1, [], length(Akpow), r), Akpow);
                U_half = reshape(U_half, [], r);
                [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);
                [Aks(m1+1:mu), f_c] = power_method(v, Aks(m1+1:mu), 1);
            end
        end

        alpha = (tensor_product(Aks(1:m1))'* U1_copy)';
        beta = (tensor_product(Aks(m1+1:order))'* V1)';

        % Solve for lambda
        Ctbeta = (beta'*C_copy)';
        lambda(k) = 1/(alpha'*Ctbeta);

        for i=1:order
            factors{i}(:,k) = Aks{i}; 
        end

        if k > 1

            % Calpha = C*alpha;

            % Calculate the new matrix D and the new subspace
            x = get_hh_reflector((beta'*C)');

            C = RHR(C,x);
            U1 = RHR(U1,x);

            % x = get_hh_reflector(Calpha);

            % C = LHR(C,x);
            % V1 = RHR(V1,x);



        end
            
        timenow = toc(timer);
        stat.deflatetime = stat.deflatetime + timenow - lap;
        lap = timenow;

    end
    
    stat.avgiter = stat.avgiter/r;
    stat.totaltime = toc(timer);

    factors(dim_order) = factors;
    
    
    factors{1} = factors{1}./vecnorm(factors{1});
    factors{2} = factors{2}./vecnorm(factors{2});
    factors{3} = factors{3}./vecnorm(factors{3});
    factors{4} = factors{4}./vecnorm(factors{4});
    factors{5} = factors{5}./vecnorm(factors{5});
    newfactor_1 = 1/4*(factors{1}*diag(sign(factors{1}(1,:)))+factors{2}*diag(sign(factors{2}(1,:)))+factors{3}*diag(sign(factors{3}(1,:)))+factors{4}*diag(sign(factors{4}(1,:))));
    newfactor_1 = newfactor_1./vecnorm(newfactor_1);
    newfactors = {newfactor_1,factors{5}};

    T_recovered = generate_lowrank_tensor(lambda.*sign(factors{1}(1,:)).*sign(factors{2}(1,:)).*sign(factors{3}(1,:)).*sign(factors{4}(1,:)), newfactors{:},[4,1]);
    err = norm(T_old-T_recovered,'fro');


    
    if nargout==1
        newfactors{1} = newfactor_1 .* lambda;
        varargout = {newfactors};
    else
        varargout = {lambda, newfactors};
    end
    if nargout==3; varargout{3} = err; end



    function [Aks, f] = power_method(U, Aks, r)

        mk = length(Aks);

        for tries = 1:ntries

            % Initialize Xk
            for j=1:mk
              Aks{j} = randn(size(Aks{j}));
              Aks{j} = Aks{j}/norm(Aks{j});  
            end
            
            for iter = 1:maxiter
            
                [Aks, f, max_shift] = power_method_iteration(U, Aks, r);
            
                if max_shift < opts.gradtol
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
              Aks_ = Aks;
            else
              stat.nrr = stat.nrr + 1;
              Aks = Aks_;
            end

        end
         
    end
    
        
        
end

function r = rank_selector(S, rank_sel)
% Select rank by looking at lues
eigenva
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


function Akpow = tensor_product(Aks)
    Akpow = Aks{1};
    for k=2:length(Aks)
        Akpow = kron(Aks{k}, Akpow);
    end
end

function [Aks, f, max_shift] = power_method_iteration(U, Aks, r)
    
    m = length(Aks);
    if m==1
        U = reshape(U, [], r);
        Ak_old = Aks{1};
        Ak_new = U * (Ak_old' * U)';
        f = Ak_new'* Ak_old;
        Ak_new = Ak_new / norm(Ak_new);
        max_shift = norm(Ak_new - Ak_old, 'inf');
        Aks{1} = Ak_new;
    else
        m2 = ceil(m/2);
        Akpow = tensor_product(Aks(1:m2));
        U_half = Akpow' * reshape(U, length(Akpow), []);
        [Aks(m2+1:m), ~, max_shift] = ...
            power_method_iteration(U_half, Aks(m2+1:m), r);
        
        Akpow = tensor_product(Aks(m2+1:m));
        U_half = pagemtimes(reshape(U, [], length(Akpow), r), Akpow);
        [Aks(1:m2), f, max_shift_] = ...
            power_method_iteration(U_half, Aks(1:m2), r);
       
        max_shift = max(max_shift, max_shift_);
    end


end

function y = get_hh_reflector(y)
% Get vector for Householder reflection
%    The last column of the corresponding Householder reflection, is a 
%    multiple of the input y.

    norm_y = norm(y);
    y(end) = y(end) + norm_y*sign(y(end));
    y = y / sqrt(abs(y(end))*norm_y);

end

function A =  RHR(A, x)
% Apply Householder reflection from the right

A = A(:, 1:end-1) - (A*x)*x(1:end-1)';

end

function A =  LHR(A, x)
% Apply Householder reflection from the right

A = A(1:end-1, :) - x(1:end-1)*(x'*A);

end
