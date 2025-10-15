clearvars
close(figure(1))
clf
addpath '../'
addpath '../helper_functions/'
load results/test_power_method.mat
colors = [0 0.3470 0.6410
          0.9790 0.8040 0.2250
          0.7500 0.2750 0.0980
          0.4660 0.6740 0.1880];

colororder(colors);

legends = {'Tensor component';
           'Spurious local maxima';
           'Other global maxima';
           'Did not converge'};

close(figure(1));
hf1 = figure(1);
hf1.Position = [100 100 500 200];
area(vecR, freq(:, [1, 2, 3, 4])); %'LineWidth',2
hold on
plot([190 190], [0 1], 'LineStyle','--', 'Color','k');
hold off
xlabel('$r$','Interpreter','latex');
ylabel('Relative frequency',...
       'Interpreter','latex');
set(gca,'FontSize', 10, 'TickLabelInterpreter','latex')
legend(legends, Interpreter="latex", Location="southwest")

% close(figure(2));
% hf2 = figure(2);
% hf2.Position = [100 100 500 200];
% error_shaded(vecR', f_hit_it_avgs, f_hit_it_stds,'LineWidth',2)
% hold on
% plot([190 190], [0 40], 'LineStyle','--', 'Color','k');
% hold off
% xlabel('$r$','Interpreter','latex');
% ylabel('$\arg\min\{k:f(x_k)\ge 0.95\}$',...
%        'Interpreter','latex');
% set(gca,'FontSize', 10, 'TickLabelInterpreter','latex')

close(figure(2));
hf2 = figure(2);
hf2.Position = [100 100 500 200];
error_shaded(vecR', xn_hit_it_avgs(:,2), xn_hit_it_stds(:,2),'LineWidth',2)
hold on
plot([190 190], [1e2 1e4], 'LineStyle','--', 'Color','k');
hold off
xlabel('$r$','Interpreter','latex');
ylabel('$I$',...
       'Interpreter','latex','FontSize',8);
set(gca,'FontSize', 10, 'TickLabelInterpreter','latex', 'YScale','log')
%legend({'$\lambda=10^{-4}$','$\lambda=10^{-10}$'},Interpreter="latex", Location="southeast")

pdfprint('results/test_power_method', hf1);
pdfprint('results/test_power_method_niter', hf2);