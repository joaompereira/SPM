function [symind, findsym, ncomb] = symmetric_indices(L,n)

    % Number of rows and columns of new matrix
    dsym = nchoosek(L+n-1,n);

    % Different indices without permutations
    symind = nchoosek(1:L+n-1,n)-(0:n-1);
    symind = symind(:,end:-1:1);
    
    % Map from all indices to indices without permutations
    if nargout>2
        ncomb = factorial(n)*ones(dsym,1);
        
        for i=2:n
            ncomb = ncomb ./ sum(symind(:,1:i)==symind(:,i),2);
        end
    end
    if nargout>1
        findsym = zeros(L * ones(1,n));
        for perm = perms(1:n)'
            findsym((symind(:,perm)-1) * (L.^(0:n-1)') + 1) = 1:dsym;
        end
    end
    
end
