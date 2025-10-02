
function score = norm_reorder_cosine_sim(A,B)
    % make columns of A,B have norm 1
    A_normed = A./ sqrt(sum(A.^2, 1));
    B_normed = B./sqrt(sum(B.^2,1));
    numcolA = size(A_normed,2);
    B_reordered = reorder_columns(A_normed,B_normed);
    summatrix = sum(A_normed.*B_reordered(:,1:size(A,2)));
    score = mean(summatrix(1:end,1:numcolA));
end



function B_reordered = reorder_columns(A, B)
    % Get the number of columns
    numCols = size(A, 2);

    % Initialize the reordered matrix B_reordered
    B_reordered = zeros(size(B));
    
    % Track columns in B that have been assigned
    assignedCols = false(1, size(B,2));

    for i = 1:numCols
        % Calculate the inner product with each column of B
        innerProducts = A(:, i)' * B;
        
        % Take the absolute values to find the maximum in magnitude
        [~, maxIdx] = max(abs(innerProducts) .* ~assignedCols);
        
        % Assign the column with the maximum inner product to B_reordered
        B_reordered(:, i) = sign(innerProducts(maxIdx)) * B(:, maxIdx);
        
        % Mark the selected column as assigned
        assignedCols(maxIdx) = true;
    end
end