function E  = noise_tensor_moments(L, n)
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

if mod(n,2)==1
    % if n is odd the expectation is 0
    E = zeros(L^n,1);
else
    % 
    Id = reshape(eye(L),[],1);
    E = Id;

    for k=3:2:n-1
        E = k * kron(Id, E);
    end
    
    E = reshape(E,L*ones(1,n));
    E = symmetrize_tensor(E, L, n);
    
end    

end