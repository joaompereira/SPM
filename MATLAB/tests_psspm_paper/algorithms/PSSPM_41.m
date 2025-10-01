function varargout = PSSPM_41(T, varargin)
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

    %% Set options here
    iP = inputParser;
    addOptional(iP,'rank', [], @(x) isempty(x) || x>0);
    addParameter(iP, 'maxiter', 5000, @(x) x>0);
    addParameter(iP,  'ntries', 5, @(x) x>0);
    addParameter(iP, 'gradtol', 1e-15, @(x) x>0);
    addParameter(iP, 'ranksel', 1e-4);
    addParameter(iP,    'ftol', 1e-2, @(x) x>0);
    addParameter(iP,   'flats', []);
    addParameter(iP, 'svals_only', false);
    addParameter(iP, 'symmetries', []);
    parse(iP, varargin{:});

    T_input = T+0;

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
    else
        [vals,~,symvec] = unique(symmetries);
        symvec = symvec';
        nudims = length(vals);
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
        ufb = uflats > 0;
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

    dim_order(symorder) = dim_order;
    T = permute(T, dim_order);

    %i_dim_order(dim_order) = 1:order;

    %size(T)
    %dims(dim_order)
    %dims(i_dim_order)

    dims = dims(dim_order);

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
        T = permute(T, [1:mc, m1+1:mu, mc+1:m1, mu+1:order]);
        
        matT_2 = reshape(T, prod(dims([1:mc, m1+1:mu])), []);
        [U2,S2,V2] = svd(matT_2, 'econ');
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

        %log10(1+1e-14-f)
        
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
            [v, ~] = eigs(U_half * U_half', 1, 1+1e-12);%, 1+1e-12
            %[v, ~] = svds(U_half, 1);
            %[v, ~] = svds(reshape(U_half, [], r), 1, 1+1e-12);
            mask = ~ufb(1,:) & ufb(2,:);
            f_c = power_method(v, nsyms .* mask, 1);

            %log10(1+1e-14-f_c)

            if order > mu
                Akpow = tensor_product(Aks, nsyms .* mask);
                U_half = Akpow' * reshape(V1, length(Akpow), []);
                U_half = reshape(U_half, [], r);
                [v, ~] = eigs(U_half * U_half', 1, 1+sqrt(eps));
                mask = ~ufb(1,:) & ~ufb(2,:);
                f_c = power_method(v, nsyms .* mask, 1);

                %log10(1+1e-14-f_c)
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
        lambda(k) = 1/(alpha'*Ctbeta);

        for i=1:nudims
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

    if ~sym_breaking
        factors(reord_sym) = factors;
        i_reord_sym(reord_sym) = 1:nudims;
        symvec = i_reord_sym(symvec);
    end

    % varargout = {lambda, factors, symvec, stat};
    T_est = generate_lowrank_tensor(lambda,factors{:}, [4,1]);
    err = norm(T_input-T_est,'fro');
    varargout = {lambda, factors,err};
    

   

    % factors(dim_order) = factors;
    % 
    % if nargout==1
    %     factors{1} = factors{1} .* lambda;
    %     varargout = {factors};
    % else
    %     varargout = {lambda, factors};
    % end
    % if nargout==3; varargout{3} = stat; end

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
            c = sqrt(1-1/nsyms) * max(1-f/2,sqrt(2*f*max(1-f,0)));
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
