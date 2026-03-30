%% standard_plotting_2nd_order.m
% Plotting script for the second-order within-host model solver.
% Uses infection_level (accumulated 1D array) instead of the full I matrix.
% Called from run_full_model.m after within_host_model_2nd_order returns.

%% iRBC abundance and gametocyte dynamics
figure(1);
hold on;
plot(x/24, infection_level, '-', 'LineWidth', 3);
xlabel('infection age x (days)');
ylabel('abundance');
yscale log;
xlim([0 800]);
ylim([1 (10.5)^6]);
set(gca, 'TickDir', 'out');

length_infection_out = find(G > G_threshold, 1, 'last');
if isempty(length_infection_out)
    length_infection_out = length(x);
end

plot(x/24, G, '-', 'LineWidth', 3);
title('parasite dynamics', 'FontWeight', 'normal');
legend('iRBC', 'gametocyte');
legend('boxoff');
set(gca, 'TickDir', 'out');

disp(['Length of infection within-host: ', ...
      num2str(h*length_infection_out/24), ' days']);

%% Strategy
figure(15);
hold on;
plot(x(1:length_infection_out)/24, 100*CC(1:length_infection_out), 'LineWidth', 3);
scatter(x(length_infection_out)/24, 100*CC(length_infection_out), 100, 'k', 'diamond', 'filled');
axis tight;
xlabel('infection age x (days)');
ylabel('transmission investment');
ytickformat('percentage');
ylim([0 50]);
xlim([0 600]);
LimitsX = xlim; LimitsY = ylim;
title('A. optimal strategies', 'FontWeight', 'Normal', ...
      'HorizontalAlignment', 'left', 'position', [LimitsX(1), LimitsY(2)]);
set(gca, 'TickDir', 'out');
box off;

%% Cumulative infectiousness
cum_inf1_time = h * cumtrapz(betaHV(G), 1) / 24;
figure(4);
hold on;
plot(x(1:length_infection_out)/24, cum_inf1_time(1:length_infection_out), '-', 'LineWidth', 3);
temp_int = h * cumtrapz(betaHV(1e10)*ones(length(x),1), 1) / 24;
plot(x/24, temp_int, '--', 'LineWidth', 3);
scatter(x(length_infection_out)/24, cum_inf1_time(length_infection_out), 100, 'k', 'diamond', 'filled');
xlim([0 600]);
ylim([0 300]);
xlabel('infection age x (days)');
ylabel('cumulative infectiousness');
LimitsX = xlim; LimitsY = ylim;
title('B. host infectiousness', 'FontWeight', 'Normal', ...
      'HorizontalAlignment', 'left', 'position', [LimitsX(1), LimitsY(2)]);
set(gca, 'TickDir', 'out');

disp(['Cumulative infectiousness (pre-threshold): ', num2str(cum_inf1_time(length_infection_out))]);

%% Expected merozoites invading (using B and M only)
immune_activation = h * find(infection_level > P.IT, 1, 'first') / 24;
if isempty(immune_activation)
    immune_activation = X_max;
end

E_merozoites = (1 - P.c) .* P.beta * (P.p*B ./ (P.p*B + P.muM));

figure(12);
hold on;
plot(x(1:length_infection_out)/24, E_merozoites(1:length_infection_out), '.-', 'LineWidth', 3);
xlabel('infection age x (days)');
ylabel('E[merozoites invading]', 'FontWeight', 'Normal');
set(gca, 'TickDir', 'out');
ylim([0.6 2.8]);
xlim([0 600]);
LimitsX = xlim; LimitsY = ylim;
title('A. expected merozoites invading', 'FontWeight', 'Normal', ...
      'HorizontalAlignment', 'left', 'position', [LimitsX(1), LimitsY(2)]);
