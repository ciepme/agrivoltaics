%% Graph Profit
clear;
close;
clc;

%load data
load('agrivoltaic_multiobjective_GA_main_data.mat', ...
    'E_set', 'P_set', 'E_pop', 'P_pop');
E_default = [E_set; E_pop];
P_default = [P_set; P_pop];

load('agrivoltaic_multiobjective_GA_main_data_low_pv_deep.mat', ...
    'E_set', 'P_set', 'E_pop', 'P_pop');
E_5p4 = [E_set; E_pop];
P_5p4 = [P_set; P_pop];

load('agrivoltaic_multiobjective_GA_main_data_2p7_deep.mat', ...
    'E_set', 'P_set', 'E_pop', 'P_pop');
E_2p7 = [E_set; E_pop];
P_2p7 = [P_set; P_pop];

load('agrivoltaic_multiobjective_GA_main_data_p1_deep.mat', ...
    'E_set', 'P_set', 'E_pop', 'P_pop');
E_p1 = [E_set; E_pop];
P_p1 = [P_set; P_pop];

%% Start Graphing

dot_size = 50;
fig1 = figure;
set(fig1,'position',[200,200,1100,800]);
scatter(P_default, E_default, dot_size, 'g');
hold on;
scatter(P_5p4, E_5p4, dot_size, 'b', 'LineWidth', 3);
scatter(P_2p7, E_2p7, dot_size, 'b', 'LineWidth', 3);
scatter(P_p1, E_p1, dot_size, 'b', 'LineWidth', 3);
scatter(star_x_position_emission, 3.5, utopia_size, 'cyan', 'filled', "pentagram");
theme(fig1, "light");
xlim([min(hand_per_year) max(hand_per_year)]);
ylabel("Net Present Value (USD)");
xlabel("Pairs of HotHands Per Year");
legend("50 cent per kWh", "5.4 cent per kWh", ...
    "2.7 cent per kWh", "0.1 cent per kWh", "Utopia Point", ...
    'Location', 'northwest')
title("Impact of HotHand Usage Frequency on Net Present Value");
saveas(fig1, "hothand.jpeg");