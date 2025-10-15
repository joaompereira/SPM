function [varargout] = MSPM_asym(T, varargin)
% MSPM_asym - Decompose an asymmetric tensor using the Multi Subspace Power Method
%
%   This function calculates the CP decomposition of an asymmetric tensor 
%   using the Multi Subspace Power Method (multiSPM). It is a specialized 
%   version of multiSPM tailored for asymmetric tensors.
%
%   ** Usage **
%       [lambda, factors, stat] = MSPM_asym(T, varargin)
%
%   ** INPUT **
%        T: Asymmetric tensor of dimension d_1 x d_2 x ... x d_n
%     varargin: Optional parameters provided as name-value pairs or a struct:
%         - 'rank': Tensor rank (optional, can be estimated if not provided).
%         - 'maxiter': Maximum number of iterations for the power method
%                      (default: 5000).
%         - 'ntries': Maximum number of attempts for the power method
%                     (default: 5).
%         - 'gradtol': Gradient tolerance for convergence (default: 1e-15).
%         - 'ranksel': Tolerance for rank selection using singular values.
%                      You may provide a function that takes the singular values
%                      as input and decides the rank, or set it to 'return_sv'. 
%                      In the latter case, the function terminates early and returns 
%                      the singular values of the first flattening. (default: 1e-4). 
%         - 'ftol': Function value tolerance for restarting the power method
%                   (default: 1e-2).
%         - 'flats': Flattenings to be considered for the decomposition. This can be
%                    provided as a cell array of index vectors, or as a logical matrix
%                    where each row indicates a flattening. For example, for a tensor
%                    T(i,j,k,l), the flattening T(ij,kl) can be indicated by [1,1,0,0]
%                    or { [1,2], [3,4] }. (default: empty, which automatically selects
%                    the flattenings that allow the highest rank).
%
%   ** OUTPUT **
%     lambda: Scaling factors for the decomposition.
%    factors: Cell array containing the factor matrices of the decomposition.
%       stat: Struct containing various statistics of the decomposition process:
%             - extracttime: Time spent on extracting tensor properties.
%             - powertime: Time spent on the power method.
%             - deflatetime: Time spent on deflation.
%             - avgiter: Average number of iterations per rank.
%             - nrr: Number of restarts during the power method.
%             - totaltime: Total runtime of the decomposition.
%
%   ** NOTE **
%   - Ensure that the 'helper_functions/' directory is added to the MATLAB path.
%
%   ** Reference **
%   K. Wang, J. M. Pereira, J. Kileel, A. Seigal, "Multi-subspace power method
%   for decomposing all tensors",
% 
%   https://github.com/joaompereira/SPM
%
%   ** Version **
%   - Version 1.0 (10/15/2025) - MIT License
        
    %% Set options here
    iP = inputParser;
    addOptional(iP,'rank', [], @(x) isempty(x) || x>0);
    addParameter(iP, 'maxiter', 5000, @(x) x>0);
    addParameter(iP,  'ntries', 5, @(x) x>0);
    addParameter(iP, 'gradtol', 1e-15, @(x) x>0);
    addParameter(iP, 'ranksel', 1e-4);
    addParameter(iP,    'ftol', 1e-2, @(x) x>0);
    addParameter(iP,   'flats', []);
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
    
    T = permute(T, dim_order);

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
        
        %% Find right-side vectors by using the second flattening
        if mc>0 && m2 > mc
            Akpow = tensor_product(Aks(1:mc));
            U_half = Akpow' * reshape(U2, length(Akpow), []);
            U_half = reshape(U_half, [], r);
            [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);
            [Aks(m1+1:mu), f_c] = power_method(v, Aks(m1+1:mu), 1);

            if order > mu
                Akpow = tensor_product(Aks(m1+1:mu));
                U_half = Akpow' * reshape(V1, length(Akpow), []);
                U_half = reshape(U_half, [], r);
                [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);
                [Aks(mu+1:order), f_c] = power_method(v, Aks(mu+1:order), 1);
            end
        else
            Akpow = tensor_product(Aks(mc+1:m1));
            U_half = Akpow' * reshape(V2, length(Akpow), []);
            U_half = reshape(U_half, [], r);
            [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);
            [Aks(mu+1:order), f_c] = power_method(v, Aks(mu+1:order), 1);            

            if m2 > mc
                Akpow = tensor_product(Aks(mu+1:order));
                U_half = pagemtimes(reshape(V1, [], length(Akpow), r), Akpow);4
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

        end
            
        timenow = toc(timer);
        stat.deflatetime = stat.deflatetime + timenow - lap;
        lap = timenow;

    end
    
    stat.avgiter = stat.avgiter/r;
    stat.totaltime = toc(timer);

    factors(dim_order) = factors;
    
    if nargout==1
        factors{1} = factors{1} .* lambda;
        varargout = {factors};
    else
        varargout = {lambda, factors};
    end
    if nargout==3; varargout{3} = stat; end

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
