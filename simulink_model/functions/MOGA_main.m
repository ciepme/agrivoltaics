%% Clear and Setup

clc; clear; close all; % Clear command window, workspace, and close figures

addpath(genpath(pwd));
agrivoltaics_variable_definition;

%User define statements
GA_SELECTOR = 1;
file_suffix = "_berry_pop100";