load results/testrd_time_order3_300.mat

addpath '../helper_functions/'

close(figure(1));

% Average over all runs
avg_time = zeros(nalgs, nvals);
low_quantile = zeros(nalgs, nvals);
high_quantile = zeros(nalgs, nvals);
freq_converged = zeros(nalgs, nvals);
for i=1:nalgs
    convergedi = error{i} < 1e-8;
    ntriesi = size(convergedi, 2);
    freq_convergedi = mean(convergedi, 2);
    freq_convergedi(isnan(time{i}(:, 1))) = NaN;
    freq_converged(i, :) = freq_convergedi;
    avg_time(i, :) = mean(time{i}, 2);
    low_quantile(i, :) = quantile(time{i}, .2, 2);
    high_quantile(i, :) = quantile(time{i}, .8, 2);
end

hf = figure(1);
hf.Position = [100 100 450 300];
ax = gca;

[hpl, hpa] = error_shaded(Lval, avg_time', low_quantile', high_quantile', 'LineWidth', 1);
hpl(2).LineStyle = '--';
hpl(3).LineStyle = '-.';

ax.XScale = 'log';
ax.YScale = 'log';

ylim([1e-3 1e2]);
xlim([Lval(1), Lval(end)]);
ylabel('Time (sec)','Interpreter','latex');

legends = cell(nalgs, 1);
for i=1:nalgs
    legends{i} = sprintf('%s (%.0f\\%%)', Algs{i,1}, 100*mean(freq_converged(i, :), "omitmissing"));
end

legend(ax,hpl, legends, 'Interpreter','latex','Location','southeast');
xlabel('$L$','Interpreter','latex');

ax.FontSize = 10;
ax.TickLabelInterpreter = "latex";
pdfprint('results/testrd_time_order3', hf);
