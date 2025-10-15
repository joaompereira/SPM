function [A,B, filler] = hogsvd(T,r)
    k = size(T,3);
    n = size(T,1);
    A_list = cell(1, k);
    M_list = cell(1,k);
    for i = 1:k
        A_list{i} = T(:, :, i);
        M_list{i} = T(:,:,i);
    end
    % form S
    S = zeros(n,n);
    for i = 1:k
        for j = i+1:k
            M_i = A_list{i};
            M_j = A_list{j};
            S = S+ 1/(k*(k-1))*(M_i*inv(M_j)+ M_j*inv(M_i));
        end
    end


    [eVec, eVal] = eig(S);
    [sortedVals, idx] = sort(diag(eVal), 'descend');
    sortedEval = diag(sortedVals);
    A = eVec(:,idx);
    A = A(:,1:r);
    B = zeros(k,r);
    for i = 1:k
        B(i,:)= diag(pinv(A)*M_list{i}*pinv(A'));
    end
    filler = S;