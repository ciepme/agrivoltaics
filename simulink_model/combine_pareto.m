%% Load Data
clear;
close;
clc;

load("agrivoltaic_multiobjective_GA_main_data_berry_pop50.mat", ...
    "social_cost_set", "social_cost_pop", "yearly_biomass_set", "yearly_biomass_pop");
social_benefit_50 = [-social_cost_set; -social_cost_pop];
yearly_biomass_50 = [yearly_biomass_set; yearly_biomass_pop];
load("agrivoltaic_multiobjective_GA_main_data_berry_pop80.mat", ...
    "social_cost_set", "social_cost_pop", "yearly_biomass_set", "yearly_biomass_pop");
social_benefit_80 = [-social_cost_set; -social_cost_pop];
yearly_biomass_80 = [yearly_biomass_set; yearly_biomass_pop];


%% Combine Pareto Fronts

social_benefit = [social_benefit_50; social_benefit_80];
yearly_biomass = [yearly_biomass_50; yearly_biomass_80];

pareto_data = [social_benefit yearly_biomass];
pareto_data_inverse = -1.*pareto_data;
pareto_indices = getNonDominated(pareto_data_inverse);
pareto_front = pareto_data(pareto_indices,:);
pareto_front = sortrows(pareto_front, 1);

x_limit = [min(pareto_data(:,1)) max(pareto_data(:,1))];
y_limit = [min(pareto_data(:,2)) max(pareto_data(:,2))];

front_size = 200;
val_size = 50;
utopia_size = 800;
fig3 = figure;
theme(fig3,"light");
hold on;
title("Pareto Front");
xlabel("Social Profit ($M)");
ylabel("Annual Raspberry Production (g/m^2)");
ylim(y_limit);
xlim(x_limit);
scatter(pareto_data(:,1), pareto_data(:,2), val_size, 'black', 'filled');
scatter(pareto_front(:,1), pareto_front(:,2), front_size, 'green', 'filled');
scatter(x_limit(2), y_limit(2), utopia_size, 'cyan', 'filled', "pentagram");
plot(pareto_front(:,1), pareto_front(:,2), 'green', 'LineWidth', 3);
legend("Values from GA Population", "Pareto Front", "Utopia Point", 'Location','southwest');
figure_name = "graphs/official_pareto" + "_combined" + ".png";
saveas(fig3,figure_name);
hold off;