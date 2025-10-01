function [symind, findsym, ncomb] = symmetric_indices(L, n)

    % Number of rows and columns of new matrix
    dsym = nchoosek(L+n-1,n);

    % Different indices without permutations
    symind = nchoosek(1:L+n-1,n)-(0:n-1);
    symind = symind(:,end:-1:1);
    
    % Map from all indices to indices without permutations
    if nargout>2
        S = ones(dsym,1);
        ncomb = factorial(n) * S;
    
        for i=2:n
            S(symind(:,i-1)~=symind(:,i)) = 0;
            S = S + 1;
            ncomb = ncomb ./ S;
        end

    end
    if nargout>1

        if n==1
            findsym= zeros(L,1);
        else
            findsym = zeros(L * ones(1,n));
        end
        findsym((symind-1) * (L.^(0:n-1)') + 1) = 1:dsym;
        perm = 1:n;
        for i=2:n
            findsym_i = findsym;
            for k=2:i
                perm([n-i+1,n-i+k]) = [n-i+k,n-i+1];
                findsym = max(findsym, permute(findsym_i, perm));
                perm([n-i+1,n-i+k]) = [n-i+1,n-i+k];
            end
        end        
        
    end


    
end
