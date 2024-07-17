clearvars
close(figure(1))
clf
addpath '../'
addpath '../helper_functions/'
load results/test_power_method_5000.mat
colors = [0 0.3470 0.6410
          0.9790 0.8040 0.2250
          0.7500 0.2750 0.0980
          0.4660 0.6740 0.1880];

colororder(colors);

legends = {'Tensor component';
           'Spurious local maxima';
           'Other global maxima';
           'Did not converge'};


hf = figure(1);
hf.Position = [100 100 500 200];
ha = area(vecR, freq(:, [1, 2, 3, 4])); %'LineWidth',2
hold on
plot([190 190], [0 1], 'LineStyle','--', 'Color','k');
xlabel('$r$','Interpreter','latex');
ylabel('Relative frequency',...
       'Interpreter','latex');
set(gca,'FontSize', 10, 'TickLabelInterpreter','latex')
legend(legends, Interpreter="latex", Location="southwest")

hold off
pdfprint('results/test_power_method', hf);