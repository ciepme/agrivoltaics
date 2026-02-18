%% Clear
clear;
clc;

%% Add Path

addpath(genpath(pwd));

%%

%open data dictionary busses
data_dict = Simulink.data.dictionary.open('agrivoltaics_v1_data_dict.sldd');
open_system("agrivoltaics_v1");

%% create variables

variable_a = [0 5; 1 5];
variable_b = [0 6; 1 6];

%%
out = sim("agrivoltaics_v1.slx");