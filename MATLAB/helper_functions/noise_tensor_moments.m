function E = noise_tensor_moments(L, n, flatten, is_sparse)
% Returns the tensor moments of a Gaussian standard vector random variable
%
% function E = noise_tensor_moments(L, n)
% function E = noise_tensor_moments(L, n, tensorform)
% function E = noise_tensor_moments(L, n, tensorform, nosym )
%
% Returns E, which is the expectation of the tensor power of a Gaussian 
% standard vector random variable.
%
%          L: the length of the vector random Gaussian variable.
%          n: the power of the tensor.
% vectorform: a flag that when set to true, E is returned in vectorized
%             and sparse form.
%      nosym: a flag that when set to true, the tensor symmetrization step
%             is not performed. This is useful when there is another
%             symmetrization step afterwards.
% 

assert(~mod(n, 1) && n>0)

if nargin<3
    flatten = false;
end
if nargin<4
    is_sparse = false;
end


if mod(n,2)==1
    % if n is odd the expectation is 0
    if is_sparse
        E = is_sparse(L^n, 1);
    else
        E = zeros(L^n,1);
    end
else
    % 
    if is_sparse
        Id = reshape(speye(L),[],1);
    else
        Id = reshape(eye(L),[],1);
    end
    E = Id;

    for k=3:2:n-1
        E = kron(Id, E);

        if is_sparse
            [i, j, v] = find(E);

            [cellind{1:k+1}] = ind2sub(L*ones(1,k+1),i);
            n_ind = horzcat(cellind{:});

            indk = repmat(i, 1, k);
            
            for p=2:k
                perm = [p:k,1:p-1,k+1]'-1;
                indk(:, p) = 1+ (n_ind-1) * L.^perm;
            end

            indk = reshape(indk, [], 1);
            E = sparse(indk, repmat(j, k, 1), repmat(v, k, 1), L^(k+1), 1);
        else
            E = E * k;
        end
        
    end

    if ~is_sparse
        E = symmetrize_tensor(E, L, n);
    end
    
end    

if ~flatten && ~is_sparse
    E = reshape(E, [L*ones(1, n), 1]);
end