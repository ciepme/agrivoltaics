%% Clear
clear;
clc;

%% Add Path

addpath(genpath(pwd));

%%

%open data dictionary busses
data_dict = Simulink.data.dictionary.open('agrivoltaics_v1_data_dict.sldd');
open_system("agrivoltaics_v1");


%% 1. Parameter Definition (Fixed values)

% Land parameters
params.land.x = 50;       % length of base (m)
params.land.y = 50;       % length of height (m)
params.land.angle = 0;    % rotation (rad)

% 2. Load the processed weather data from the .mat file
weather_data = load('pv_inputs.mat'); 

% 3. Assign the loaded data into parameter structure
params.weather.DNI = weather_data.DNI;       
params.weather.DHI = weather_data.DHI;       
params.weather.beta_s = weather_data.beta_s; 
params.weather.phi_s = weather_data.phi_s;

% PV parameters
params.PV.n_p = 0.2;      % panel efficiency

% Crop & Econ parameters
params.crop.elec_price = 0.5; 
params.crop.crop_price = 14.46; %USD/kg
params.crop.HI = 0.3; %this is harvest index, for raspberries it is roughly .3, so 30% of the plant weight is berries- changes based on crop choice
params.crop.MC = .85; %this is moisture content, raspberries are about 85% water
params.crop.RUE = 2.0; %radiation use efficiency, g of biomass per MJ of light, crop dependent
params.crop.k =  0.65; %light extinction coefficient, crop dependent
params.crop.LAI = 3.0; %lead area index (sq meters of leaves per sq meter of ground, crop dependent
params.crop.capital_cost = 1500; %USD per kW in net costs
params.crop.payback_period = 20;

% environmental parameters
base_dir = fileparts(mfilename('fullpath'));

if isempty(base_dir) % in case there are issues with referencing the main dir
    base_dir = pwd;
end
data_dir = fullfile(base_dir,'parameterData');
ci_avg_hourly = get_july1_hourly_carbon_intensity(data_dir);
ci_July_1 = ci_avg_hourly; 
params.env.ci_marginal_hourly_miso = ci_July_1;


%% 2. Design Variables

% Panel layout variables
var.PV.z_p = 2;           % panel height (m) 
var.PV.l_p = 1;           % panel length (m)
var.PV.w_p = 1;           % panel width (m)
var.PV.phi = 0;           % azimuth (rad)
var.PV.sigma = pi/4;      % tilt (rad)
var.PV.psi = 0;           % field layout angle (rad)
var.PV.y_p = 1;           % row distance (m)
var.PV.x_p = 0.1;         % panel distance (m)

%% 3.Simulink Bus Objects
Simulink.Bus.createObject(params);
Simulink.Bus.createObject(var);

% Rename the auto-generated buses to match your model ports
% createObject makes 'slBus1', 'slBus2' etc. We need to rename them.

params_bus = slBus1; 
var_bus = slBus2; 

%out = sim("agrivoltaics_v1.slx");

function ci_avg_hourly = get_july1_hourly_carbon_intensity(data_dir)
emissions_2024 = readtable(fullfile(data_dir, ...
    'generated_emissions_MISO_20240701T0400_20240702T0359_co2_for_electricity_5m.csv'));
emissions_2025 = readtable(fullfile(data_dir, ...
    'generated_emissions_MISO_20250701T0400_20250702T0359_co2_for_electricity_5m.csv'));

ci_2024_5min = emissions_2024.total_co2_marginal_intensity_lbs_per_mwh;
ci_2025_5min = emissions_2025.total_co2_marginal_intensity_lbs_per_mwh;
ci_avg_5min = (ci_2024_5min + ci_2025_5min) / 2;

% 12 five-minute points per hour --> 24 hourly vals
ci_avg_hourly = mean(reshape(ci_avg_5min, 12, []), 1).'; % making this a matrix, transposing, and taking average of each col corresponding to each hour
% disp(ci_avg_hourly); % lbs/MWh
end