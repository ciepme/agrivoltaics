%% Load Data
clear;
close;
clc;

addpath(genpath(pwd));

% prep model

agrivoltaics_variable_definition;

%modify parameters
agriParams.crop_elec_price = 0.054; % changed from 0.5; 0.054 is new nominal

load("agrivoltaic_multiobjective_GA_main_data_berry_pop50.mat", ...
    "social_cost_set", "social_cost_pop", "yearly_biomass_set", ...
    "ga_final_pop_set", "ga_set", "yearly_biomass_pop");
social_cost_set_50 = social_cost_set;
social_cost_pop_50 = social_cost_pop;
ga_set_50 = ga_set;
ga_final_pop_set_50 = ga_final_pop_set;
social_benefit_50 = [-social_cost_set; -social_cost_pop];
yearly_biomass_50 = [yearly_biomass_set; yearly_biomass_pop];

load("agrivoltaic_multiobjective_GA_main_data_berry_pop80.mat", ...
    "social_cost_set", "social_cost_pop", "yearly_biomass_set", ...
    "ga_final_pop_set", "ga_set", "yearly_biomass_pop");
social_cost_set_80 = social_cost_set;
social_cost_pop_80 = social_cost_pop;
ga_set_80 = ga_set;
ga_final_pop_set_80 = ga_final_pop_set;
social_benefit_80 = [-social_cost_set; -social_cost_pop];
yearly_biomass_80 = [yearly_biomass_set; yearly_biomass_pop];

load("agrivoltaic_multiobjective_GA_main_data_berry_pop100.mat", ...
    "social_cost_set", "social_cost_pop", "yearly_biomass_set", ...
    "ga_final_pop_set", "ga_set", "yearly_biomass_pop");
social_cost_set_100 = social_cost_set;
social_cost_pop_100 = social_cost_pop;
ga_set_100 = ga_set;
ga_final_pop_set_100 = ga_final_pop_set;
social_benefit_100 = [-social_cost_set; -social_cost_pop];
yearly_biomass_100 = [yearly_biomass_set; yearly_biomass_pop];


%% Combine Pareto Fronts

social_benefit = [social_benefit_50; social_benefit_80; social_benefit_100];
yearly_biomass = [yearly_biomass_50; yearly_biomass_80; yearly_biomass_100];

pareto_data = [social_benefit yearly_biomass];
pareto_data_inverse = -1.*pareto_data;
pareto_indices = getNonDominated(pareto_data_inverse);
pareto_front = pareto_data(pareto_indices,:);
pareto_front = sortrows(pareto_front, 1);

x_limit = [min(pareto_data(:,1)) max(pareto_data(:,1))];
y_limit = [min(pareto_data(:,2)) max(pareto_data(:,2))];

% find optimal design
social_design = [1.05215 678.568];
biomass_design = [0.381511 1319.81];
balanced_design = [0.832745 840.469];

tol = 1e-3;

social_index = find(all(abs(pareto_front - social_design) < 1e-3, 2), 1);
biomass_index = find(all(abs(pareto_front - biomass_design) < 1e-2, 2), 1);
balanced_index = find(all(abs(pareto_front - balanced_design) < 1e-3, 2), 1);
intesting_designs = [pareto_front(social_index,:); pareto_front(biomass_index,:); pareto_front(balanced_index,:)];

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
scatter(intesting_designs(:,1), intesting_designs(:,2), front_size, 'blue', 'd', 'filled');
scatter(x_limit(2), y_limit(2), utopia_size, 'cyan', 'filled', "pentagram");
plot(pareto_front(:,1), pareto_front(:,2), 'green', 'LineWidth', 3);
scatter(intesting_designs(:,1), intesting_designs(:,2), front_size, 'blue', 'd', 'filled');
legend("Values from GA Population", "Pareto Front", "Designs of Interest", "Utopia Point", 'Location','southwest');
figure_name = "graphs/official_pareto" + "_combined" + ".png";
saveas(fig3,figure_name);
hold off;

%% Find Indices

social_index_full = find(all(abs(pareto_data - social_design) < 1e-3, 2), 1);
biomass_index_full = find(all(abs(pareto_data - biomass_design) < 1e-2, 2), 1);
balanced_index_full = find(all(abs(pareto_data - balanced_design) < 1e-3, 2), 1);

% balanced
if balanced_index_full < length(social_benefit_50)
    variable_balanced_index = balanced_index_full - length(ga_set_50);
    variables_balanced = ga_final_pop_set_50(variable_balanced_index,:);
    agrivoltaic_wrapper(variables_balanced, agriParams)
end

% social
if social_index_full <= length(social_cost_set_50)
    variable_social_index = social_index_full;
    variables_social = ga_set_50(social_index_full,:);
    agrivoltaic_wrapper(variables_social, agriParams)
end

% biomass
if biomass_index_full > (length(social_benefit_50) + length(social_benefit_80))
    variable_biomass_index = biomass_index_full - length(social_benefit_50) - length(social_benefit_80) - length(ga_set_100);
    variables_biomass = ga_final_pop_set_100(variable_biomass_index,:);
    agrivoltaic_wrapper(variables_biomass, agriParams)
end