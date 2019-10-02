function [V, D] = eig2(T)
% EIG2 Return eigenvalue decomposition of symmetric matrix
% This algorithm ensures T is symmetric and returns D in vector form

% Ensure matrix is symmetric up to machine precision
% Matlab is much faster if this operation is performed,
% and the eigenvalues returned are real
if ~issymmetric(T)
    T = (T + T')/2;
end

[V, D] = eig(T,'vector');

%Order eigenvalues/eigenvectors in decreasing order
[D, ind] = sort(abs(D),'descend');
V = V(:,ind);

end

