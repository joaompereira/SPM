function   [A,B,filler] = pca_c(T,my_sample_list,pca_comp_num,r)
    % n = size(my_sample_list{1},2);
    k = size(my_sample_list,1);
    
 
    pca_comp_cell = cell(1,k);
    for i =1:k
        coeff = pca(my_sample_list{i},'NumComponents',pca_comp_num(i));
        pca_comp_cell{i}=coeff;   
    end

    pre_cluster_vectors=horzcat(cell2mat(pca_comp_cell));
    
    numClusters = r;  
    [~,C] = kmeans(pre_cluster_vectors', numClusters,'Distance','cosine');
    A = C'./vecnorm(C');

    krp = reshape(A, [], 1, r) .* reshape(A, 1, [], r);
    size(krp)
    Tmatrix = reshape(T,[],k);
    krp = reshape(krp, [], r);
    Btranspose = pinv(krp)*Tmatrix;
    B = Btranspose';
    filler = [];
    
end
