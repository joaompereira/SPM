function varargout = multiSPM(T, varargin)
% multiSPM - Decompose any tensor using the Multi Subspace Power Method
%
%   This function calculates the CP decomposition of a symmetric, partial
%   symmetric, or asymmetric tensor using the Multi Subspace Power Method
%   (multiSPM). It supports various options for tensor rank estimation, 
%   symmetry handling, and flattening strategies.
%
%   ** Usage **
%       [lambda, factors, symvec, stat] = multiSPM(T, varargin)
%
%   ** INPUT **
%        T: Tensor of dimension d_1 x d_2 x ... x d_n
%     varargin: Optional parameters provided as name-value pairs or a struct:
%      - 'rank': Tensor rank (optional, can be estimated if not provided).
%      - 'maxiter': Maximum number of iterations for the power method
%                  (default: 5000).
%      - 'ntries': Maximum number of attempts for the power method
%                  (default: 5).
%      - 'gradtol': Gradient tolerance for convergence (default: 1e-15).
%      - 'ranksel': Tolerance for rank selection using singular values
%                   You may provide a function that takes the flattening
%                   singular values as input and decides the rank, or set
%                   it to 'return_sv'. In the last case, multiSPM terminates 
%                   early and returns the singular values of the first
%                   flattening. (default: 1e-4). 
%      - 'ftol': Function value tolerance for restarting the power method
%                 (default: 1e-2).
%      - 'flats': Flattenings to be considered for the decomposition. This
%                 can be provided as a cell array of index vectors, or as a 
%                 logical matrix where each row indicates a flattening. For
%                 example, for a tensor T(i,j,k,l), the flattening T(ij,kl)
%                 can be indicated by [1,1,0,0] or { [1,2], [3,4] }.
%                 (default: empty, which automatically selects the 
%                           flattenings that allow the highest rank).
%       - 'symmetries': Tensor symmetries. You may inform the symmetries of
%                         the tensor in two ways:
%                       * As a cell array, where each cell contains the 
%                         indices of the modes that are symmetric. For 
%                         example, for a tensor T(i,j,k,l) with symmetries 
%                         T(i,j,k,l) = T(j,i,k,l), you would provide 
%                         {[1,2],[3],[4]}.
%                       * As a vector of length n (tensor order), where the 
%                         same number indicates symmetry. For example, for 
%                         the same tensor T(i,j,k,l) with symmetries 
%                         T(i,j,k,l) = T(j,i,k,l), you would provide [1,1,2,3].
%                       * As a vector indicating the number of symmetric modes.
%                         For example, for the same tensor T(i,j,k,l) with 
%                         symmetries T(i,j,k,l) = T(j,i,k,l), you would provide
%                         [2,1,1]. This assumes the tensor modes are ordered
%                         according to their symmetry groups. You may not use 
%                         this option to describe asymmetric tensors, as the 
%                         vector [1, 1, ...,1] would clash with the previous 
%                         option, where it refers to a fully symmetric tensor.
%                         (default value: 1:n (asymmetric tensor)).
%
%   ** OUTPUT **
%     lambda: Scaling factors for the decomposition.
%    factors: Cell array containing the factor matrices of the decomposition.
%     symvec: Vector indicating the symmetry structure of the tensor.
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
    addParameter(iP, 'symmetries', []);
    parse(iP, varargin{:});

    opts = iP.Results;
    r = opts.rank;
                                
    timer = tic;

    dims = size(T);
    order = length(dims);

    %% Processing tensor symmetries
    symmetries = opts.symmetries;
    if isempty(symmetries)
        symvec = 1:order;
        nudims = order;
    elseif iscell(symmetries)
        symvec = zeros(1, order);
        lsym = length(symmetries);
        for i=1:lsym
            assert(all(symvec(symmetries{i})==0))
            symvec(symmetries{i}) = i;
        end
        nudims = lsym + sum(~symvec);
        symvec(~symvec) = lsym+1:nudims;
    elseif length(symmetries) == order
        [vals,~,symvec] = unique(symmetries);
        symvec = symvec';
        nudims = length(vals);
    else
        nudims = length(symmetries);
        symvec = repelem(1:nudims, symmetries);
    end
    
    udims = zeros(1, nudims);
    nsyms = zeros(1, nudims);
    symorder = [];
    for i=1:nudims
        syms_i = find(symvec == i);
        symorder = [symorder, reshape(syms_i, 1, [])];
        udims(i) = dims(syms_i(1));
        assert(all(dims(syms_i)==udims(i)));
        nsyms(i) = length(syms_i);
    end

    %% Processing flattenings
    flats = opts.flats;
    if isempty(flats)
        [uflats, flat_ranks] = find_biggest_flattenings_ps(udims, nsyms, 3);
    else
        if iscell(flats)
            cell_flats = flats;
            flats = false(length(cell_flats), order);
            for i=1:length(cell_flats)
                flats(i, cell_flats{i}) = true;
            end
        end
        for i=size(flats,1):-1:1
            uflats(i, :) = accumarray(symvec, flats(i, :), udims);
            left_rank = sv_dimension(udims, uflats(i, :));
            left_rank = left_rank - sum(udims(uflats(i, :)>0));
            right_rank = sv_dimension(udims, nsyms - uflats(i, :));
            flat_ranks(i) = min(left_rank, right_rank);
        end
        
    end
    
    [uflats, max_rank] = process_uflats(nsyms, uflats, flat_ranks);
    cs_nsyms = [0, cumsum(nsyms)];
    nflats = size(uflats, 1);
    sym_breaking = nflats == 1;
    unord_flats = false(nflats, order);
    for i=1:nudims
        for j=1:nflats
            unord_flats(j, cs_nsyms(i)+1:cs_nsyms(i)+uflats(j, i)) = true;
        end
    end

    assert(isempty(r) || r<= max_rank);
    
    m1 = sum(unord_flats(1, :));
    
    if sym_breaking
        comp_flat = false(1, order);
        for i=find(uflats==0)
            comp_flat(cs_nsyms(i)+1:cs_nsyms(i+1)) = true;
        end
        dim_order = [find(unord_flats),...
                     find(~comp_flat & ~unord_flats),...
                     find(comp_flat)];
    else
        ufb = boolean(uflats) > 0;
        reord_sym = [find(ufb(1, :) & ufb(2, :)),...
                     find(ufb(1, :) & ~ufb(2, :)),...
                     find(~ufb(1, :) & ufb(2, :)),...
                     find(~ufb(1, :) & ~ufb(2, :))];
        symorder = reord_sym(symorder);
        uflats = uflats(:, reord_sym);
        ufb = ufb(:, reord_sym);
        udims = udims(reord_sym);

        m2 = sum(nsyms(ufb(1, :)));
        mc = sum(nsyms(ufb(1, :) & ufb(2, :)));
        mu = sum(nsyms(ufb(1, :) | ufb(2, :)));

        unord_flats = repelem(uflats, 1, nsyms);

        dim_order = [find(unord_flats(1, :) & unord_flats(2, :)),...
                     find(unord_flats(1, :) & ~unord_flats(2, :)),...
                     find(~unord_flats(1, :) & unord_flats(2, :)),...
                     find(~unord_flats(1, :) & ~unord_flats(2, :))];
    end

    %% Permute dimensions, to make further calculations easier
    dim_order(symorder) = dim_order;
    T = permute(T, dim_order);

    dims = dims(dim_order);

    %% Compute SVD(s) of flattening(s)
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

    if ~sym_breaking
        % If flattening does not break symmetry,
        % use another flattening to calculate all components that
        % correspond to each other
        T = permute(T, [1:mc, m1+1:mu, mc+1:m1, mu+1:order]);
        
        matT_2 = reshape(T, prod(dims([1:mc, m1+1:mu])), []);
        [U2, S2, V2] = svd(matT_2, 'econ');
        U2 = U2(:, 1:r);
        V2 = V2(:, 1:r);
    end

    % Pre-allocation of X and lambda
    factors = arrayfun(@(dim) zeros(dim, r), udims,...
            	       'UniformOutput',false);
    
    lambda = zeros(1,r);
    
    lap = toc(timer);
    stat.extracttime = lap;
    stat.powertime = 0;
    stat.deflatetime = 0;
    stat.avgiter = 0;
    stat.nrr = 0;

    Aks = arrayfun(@(dim) zeros(dim, 1), udims,...
            	       'UniformOutput',false);

    ntries = opts.ntries;
    maxiter = opts.maxiter;

    for k = r:-1:1
          
        f = power_method(U1, uflats(1,:), k);        
        %% Find right-side vectors by using the second flattening
        if sym_breaking
            if any(~uflats)
                Akpow = tensor_product(Aks, (nsyms - uflats) .* ~~uflats);
                U_half = Akpow' * reshape(V1, length(Akpow), []);
                U_half = reshape(U_half, [], r);
                [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);
                f_c = power_method(v, nsyms .* ~uflats, 1);
            end
        elseif mc>0 && m2 > mc
            mask = ufb(1,:) & ufb(2,:);
            Akpow = tensor_product(Aks, nsyms .* mask);
            U_half = Akpow' * reshape(U2, length(Akpow), []);
            U_half = reshape(U_half, [], r);
            [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);
            mask = ~ufb(1,:) & ufb(2,:);
            f_c = power_method(v, nsyms .* mask, 1);

            if order > mu
                Akpow = tensor_product(Aks, nsyms .* mask);
                U_half = Akpow' * reshape(V1, length(Akpow), []);
                U_half = reshape(U_half, [], r);
                [v, ~] = eigs(U_half * U_half', 1, 1+sqrt(eps));
                mask = ~ufb(1,:) & ~ufb(2,:);
                f_c = power_method(v, nsyms .* mask, 1);
            end
        else
            mask = ufb(1,:) & ~ufb(2,:);
            Akpow = tensor_product(Aks, nsyms .* mask);
            U_half = Akpow' * reshape(V2, length(Akpow), []);
            U_half = reshape(U_half, [], r);
            [v, ~] = eigs(U_half * U_half', 1, 1+sqrt(eps));
            mask = ~ufb(1,:) & ~ufb(2,:);
            f_c = power_method(v, nsyms .* mask, 1);          

            if m2 > mc
                Akpow = tensor_product(Aks, nsyms .* mask);
                U_half = pagemtimes(reshape(V1, [], length(Akpow), r), Akpow);
                U_half = reshape(U_half, [], r);
                [v, ~] = eigs(U_half * U_half', 1, 1+sqrt(eps));
                mask = ufb(2,:) & ~ufb(1,:);
                f_c = power_method(v, nsyms .* mask, 1);
            end
        end
        
        alpha = (tensor_product(Aks, uflats(1,:))'* U1_copy)';
        beta = (tensor_product(Aks, nsyms - uflats(1,:))'* V1)';

        % Solve for lambda
        Ctbeta = (beta'*C_copy)';
        lambda(k) = norm(alpha)*norm(beta)/(alpha'*Ctbeta);

        for i=1:nudims
            factors{i}(:,k) = Aks{i}; 
        end

        if k > 1
            %% Deflation step
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

    if ~sym_breaking
        factors(reord_sym) = factors;
        i_reord_sym(reord_sym) = 1:nudims;
        symvec = i_reord_sym(symvec);
    end

    varargout = {lambda, factors, symvec, stat};

    function f = power_method(U, nsyms, r)

        for tries = 1:ntries
            
            Ak_inds = find(nsyms);
            % Initialize Xk
            for jj=Ak_inds
              Aks{jj} = randn(size(Aks{jj}));
              Aks{jj} = Aks{jj}/norm(Aks{jj});  
            end
            
            for iter = 1:maxiter
            
                [Aks(Ak_inds), f, max_shift] = ...
                    power_method_iteration(U, Aks(Ak_inds), nsyms(Ak_inds), r);
            
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

function [uflats, max_rank, uf_ind] = process_uflats(nsyms, uflats, ranks)

if ~issorted(ranks,'descend')
    [ranks, I] = sort(ranks, 'descend');
    uflats = uflats(I, :);
end

for j=1:size(uflats, 1)
    if any(and(uflats(j, :),nsyms - uflats(j, :)))
        uflats = uflats(j, :);
        max_rank = ranks(j);
        uf_ind = j;
        return
    end
    for i = 1:j-1
        if ~all(uflats(i, :) + uflats(j, :) == nsyms)
            if ranks(i)>ranks(j)
                uf_ind = [i, j];
                uflats = uflats(uf_ind, :);
                max_rank = ranks(j);
            else
                uf_ind = [j, i];
                uflats = uflats(uf_ind, :);
                max_rank = ranks(i);
            end
            return
        end
    end
end

err_msg = "Provide either a symmetry breaking flattening,\n" + ...
          "or two flattenings that are not complementary";

error(sprintf(err_msg));
   
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


function Akpow = tensor_product(Aks, nsyms)
    if nargin < 2
        nsyms = ones(1, length(Aks));
    end
    Akpow = 1;
    for i=1:length(Aks)
        for j=1:nsyms(i)
            Akpow = reshape(Akpow * Aks{i}.', [], 1);
        end
    end
end

function [Aks, f, max_shift] = power_method_iteration(U, Aks, nsyms, r)
    
    m = length(Aks);
    if m==1
        Ak_old = Aks{1};
        if nsyms == 1
            UAk = reshape(U, [], r);
        else
            n = size(Ak_old, 1);
            UAk = reshape(tensor_product(Aks, nsyms - 1)' * ...
                                    reshape(U, [], n * r), n, r);
        end
        Ak_new = UAk * (Ak_old' * UAk)';
        f = Ak_new'* Ak_old;
        if nsyms > 1
            if f <= 2/3
                c = sqrt(1-1/nsyms) * (1-f/2);
            else
                c = sqrt(1-1/nsyms) * sqrt(2*f*max(1-f,0));
            end
            Ak_new = Ak_new + c * Ak_old;
        end
        Ak_new = Ak_new / norm(Ak_new);
        max_shift = norm(Ak_new - Ak_old, 'inf');
        Aks{1} = Ak_new;
    else
        m2 = ceil(m/2);
        Akpow = tensor_product(Aks(1:m2), nsyms(1:m2));
        U_half = Akpow' * reshape(U, length(Akpow), []);
        [Aks(m2+1:m), ~, max_shift] = ...
            power_method_iteration(U_half, Aks(m2+1:m), nsyms(m2+1:m), r);
        
        Akpow = tensor_product(Aks(m2+1:m), nsyms(m2+1:m));
        U_half = pagemtimes(reshape(U, [], length(Akpow), r), Akpow);
        [Aks(1:m2), f, max_shift_] = ...
            power_method_iteration(U_half, Aks(1:m2), nsyms(1:m2), r);
       
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
