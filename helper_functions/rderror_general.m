function [err, perm, E] = rderror_general(A_est, A_true, W_est, W_true, n)
% Returns the error of the X_est in terms of X_true
% Code inspired by alignment algorithms in
%       https://github.com/NicolasBoumal/HeterogeneousMRA
                 
    % Calculate distance between all pairs of vectors.
    % Use abs because rank decompositions are unique up to sign flip
    n_true = numel(A_true);
    n_est = numel(A_est);
    E = zeros(n_true,n_est);
    for i=1:length(A_true)
      for j = 1:length(A_est)
        At = A_true{i};
        Rt = size(At,2);
        Wt = W_true{i};
        Ae = A_est{j};
        Re = size(Ae,2);
        We = W_est{j};
        [Ut,St,Vt] = svd(At,0);
        [Ue,Se,Ve] = svd(Ae,0);
        Wt = tucker_product(St*Vt',Wt,n);
        We = tucker_product(Se*Ve',We,n);
        W1 = Wt - tucker_product(Ut'*Ue,We,n);
        M1 = (Ut*(Ut'*Ue)-Ue)';
        W2 = tucker_product(M1*Ut,Wt,n)-...
                tucker_product(M1*Ue,We,n);
        E(i,j) = sqrt(norm(W1(:))^2 + W2(:)'*We(:));
      end
    end
    
    % Use an out-sourced implementation of Hungarian algorithm to find
    % right permutation of indices that minimize the distances above
    [perm, err] = munkres(E);
    

end