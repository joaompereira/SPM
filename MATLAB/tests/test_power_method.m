clearvars
clc
addpath '../'
addpath '../helper_functions/'

L = 20;
n2 = 2;
vecR = 120:1:200;
atries = 10000;
maxiter = 1e4;

freq = zeros(size(vecR));

if n2<=4
  cn = sqrt(2*(n2-1)/n2);
else
  cn = (2-sqrt(2))*sqrt(n2);
end
                        
for Riter=1:length(vecR)
    
R = vecR(Riter)

freqa = 0;
                        
parfor aiter=1:atries

    a_true = randn(L,R);
    a_true = a_true./vecnorm(a_true);

    decG = decomposition((a_true'*a_true).^n2);

    Ak = randn(L,1);
    Ak = Ak/norm(Ak);

    fval = 0;

    for tries=1:maxiter

        Ax = Ak'*a_true;

        Ak_new = ((Ax.^(n2-1)).*a_true)*(decG \ (Ax.^n2)');

        f = Ak_new'*Ak;

        % Determine optimal shift
        % Sometimes due to numerical error f can be greater than 1
        f_ = max(min(f,1),.5);
        clambda = sqrt(f_*(1-f_));
        shift = cn*clambda;

        % Shifted power method
        Ak_new = Ak_new + shift*Ak;
        Ak_new = Ak_new/norm(Ak_new);

        if norm(Ak - Ak_new) < 1e-14
            % Algorithm converged
            Ak = Ak_new;
            break
        else
            Ak = Ak_new;
        end
    end

    err = 2*min(1-abs(Ak'*a_true));

    if err>1e-3
       freqa = freqa + 1;
    end

end

freq(Riter) = 1-freqa/atries;

end

clear decG
save results/test_power_method.mat

% Plot results
if ispc
    hf = figure(1);
    ha = plot(vecR, freq, 'LineWidth',2);
    hold on
    plot([190 190], [0 1], 'LineStyle','--');
    xlabel('$R$','Interpreter','latex');
    ylabel('Relative frequency',...
           'Interpreter','latex');
    pdfprint('results/test_power_method', hf);
end
