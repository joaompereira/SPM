function [A,B,filler] = tensorlab_nls_better_initialization(T, r)
    model=struct;
    [U,V,~] = svd_ortho(T,r);
    model.variables.u=U;
    model.variables.v=V;
    model.factors.U='u';
    model.factors.V='v';
    model.factorizations.myfac.data = T;
    model.factorizations.myfac.cpd={'U','U','V'};
    sol = ccpd_nls(model);
    A = sol{1};
    B = sol{2}.* vecnorm(A).^2;
    A = A./vecnorm(A);
    filler = [];
end