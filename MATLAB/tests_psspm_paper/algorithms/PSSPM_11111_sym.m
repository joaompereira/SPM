function [lambda_new, newfactors] = PSSPM_11111_sym(T, R)
    
    [lambda, factors] = asym_SPM(T,'rank',R);
    newfactor_1 = 1/4*(factors{1}*diag(sign(factors{1}(1,:)))+factors{2}*diag(sign(factors{2}(1,:)))+factors{3}*diag(sign(factors{3}(1,:)))+factors{4}*diag(sign(factors{4}(1,:))));
    newfactor_1 = newfactor_1./vecnorm(newfactor_1);
    newfactors = {newfactor_1,factors{5}};
    lambda_new = lambda.*sign(factors{1}(1,:)).*sign(factors{2}(1,:)).*sign(factors{3}(1,:)).*sign(factors{4}(1,:));

    


    
