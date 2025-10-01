function [varargout] = eig2(T, order, varargin)
% EIG2 Return eigenvalue decomposition of symmetric matrix
% This algorithm ensures T is symmetric and returns D in vector form
% Matlab is much faster in doing a symmetric matrix decomposition, and the
% eigenvalues are more accurate

if nargin<2
    order = 'descend';
end

% If not symmetric, symmetrize matrix
if ~issymmetric(T)
    T = (T + T')/2;
end

if nargout==1
    D = eig(T,'vector');
else
    [V, D] = eig(T,'vector');
end

%Order eigenvalues/eigenvectors in decreasing order
if ~issorted(abs(D), order)
    [~, ind] = sort(abs(D),order);
    D = D(ind);
    if nargout>1
        V = V(:,ind);
    end
end

if nargout==1
    varargout{1} = D;
else
    varargout{2} = D;
    varargout{1} = V;
end

end

