function [varargout] = debias_even_moments(L, sigma, varargin)

M_est = varargin;

n2 = length(M_est);

NM = cell(floor(n2),1);
for i = 1:n2
    NM{i} = noise_tensor_moments(L, 2*i);
end

for k=1:n2
    for i=1:k
        if i==k
           M_est{k} = M_est{k} - sigma^(2*k) * NM{i};
        else
           C = nchoosek(2*k,2*i)*sigma^(2*i);
           T = reshape((C*NM{i}(:))*M_est{k-i}(:)',L*ones(1,2*k));
           M_est{k} = M_est{k} - T;
        end
    end
    
    M_est{k} = symmetrize_tensor(M_est{k});
end

varargout = M_est;