function [lambda,newfactors,err] = tensorlab_minf_11111_avg(T, r)
    n  = size(T,1);
    k = size(T,5);
    model=struct;
    model.variables.u=randn(n,r);
    model.variables.w=randn(n,r);
    model.variables.c= randn(n,r);
    model.variables.d = randn(n,r);
    model.variables.v=randn(k,r);
    model.factors.U='u';
    model.factors.V='v';
    model.factors.W = 'w';
    model.factors.C = 'c';
    model.factors.D = 'd';

    model.factorizations.myfac.data = T;
    model.factorizations.myfac.cpd={'U','W','C','D','V'};
    sol_old = ccpd_minf(model);

    sol = sol_old;

    
    sol{1} = sol{1}./vecnorm(sol{1});
    sol{2} = sol{2}./vecnorm(sol{2});
    sol{3} = sol{3}./vecnorm(sol{3});
    sol{4} = sol{4}./vecnorm(sol{4});
    sol{5} = sol{5}./vecnorm(sol{5});
    newfactor_1 = 1/4*(sol{1}*diag(sign(sol{1}(1,:)))+sol{3}*diag(sign(sol{3}(1,:)))+sol{4}*diag(sign(sol{4}(1,:)))+sol{5}*diag(sign(sol{5}(1,:))));
    newfactor_1 = newfactor_1./vecnorm(newfactor_1);
    newfactors = {newfactor_1,sol{2}};
    
    lambda = vecnorm(sol_old{1}).*vecnorm(sol_old{2}).*vecnorm(sol_old{3}).*vecnorm(sol_old{4}).*vecnorm(sol_old{5});
    T_recovered = generate_lowrank_tensor(lambda.*sign(sol{1}(1,:)).*sign(sol{5}(1,:)).*sign(sol{3}(1,:)).*sign(sol{4}(1,:)), newfactors{:},[4,1]);
    err = norm(T-T_recovered,'fro');
end