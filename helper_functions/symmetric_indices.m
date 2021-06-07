function [symind, findsym, ncomb] = symmetric_indices(L,n)

% Number of rows and columns of new matrix
    dsym = nchoosek(L+n-1,n);

    % Different indices without permutations
    symind = nchoosek(1:L+n-1,n)-(0:n-1);
    symind = symind(:,end:-1:1);
    

    % Map from all indices to indices without permutations
    if nargout>1
        findsym = zeros(ones(1,n)*L);
        ncomb = zeros(dsym,1);

        for i=1:dsym
            % Scaling to keep have the same dot product in the new space
            if nargout>2
                ncomb(i) = nperm(symind(i,:));
            end
            for perm = perms(1:n)'
                commalist = num2cell(symind(i,perm));
                findsym(commalist{:}) = i;
            end
        end
    end