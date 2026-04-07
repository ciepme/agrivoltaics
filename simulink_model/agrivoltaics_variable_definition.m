%% Add Path
clear
addpath(genpath(pwd));

%%

%open data dictionary busses
% data_dict = Simulink.data.dictionary.open('agrivoltaics_v1_data_dict.sldd');
% open_system("agrivoltaics_v1");

%% Variable Bound Definition
% (1) = PV_z_p; panel height above the ground (m)
% (2) = PV_l_p; panel length (m)
% (3) = PV_w_p; panel width (m)
% (4) = PV_phi; azimuth angle (radians) - relative to true South, going ccw e.g. pi/2 rad is East
% (5) = PV_sigma; tilt angle (radians) - fixed sloping angle of PV relative to horizontal (xy) plane
% (6) = PV_y_p; Distance between rows (m) of panels
% (7) = PV_x_p; Distance between panels in a pair/row (m)
% Order: [Height, Length, Width, Azimuth, Tilt, Row Gap, Panel Gap]
lb = [2.5, 1.0, 1.0, -pi/2, 0,    2.5,  0.1]; 
ub = [4.5, 2.5, 1.5,  pi/2, pi/2, 10.0, 1.0];

%% Tracking System Definition (Fixed vs Single Axis)
% 0 = Fixed Axis, 1 = Single-Axis Tracking 
agriParams.tracking_mode = 0; 

% Max rotation angle for the single-axis tracker (mechanical limit before
% hitting motor housing or sturcture
agriParams.PV_max_tilt = 60 * (pi/180); % Convert degrees to radians

% tracking angles for optimizer matrix definition
% 4x24 matrix (4 seasons, 24 hours). 
agriVar.tracking_angles = zeros(4, 24);
%% Parameter Definition (Fixed values)

% Land parameters
agriParams.land_x = 50;       % length of base (m)
agriParams.land_y = 50;       % length of height (m)
agriParams.land_angle = 0;    % rotation (rad)

% 2. Load the processed weather data from the .mat file
load('pv_weather_4seasons.mat', 'weather_struct');
agriParams.weather = weather_struct;

% 3. Assign the loaded data into parameter structure (changed to seasonal
% input with the weather_struct
% agriParams.weather_DNI = weather_data.DNI;       
% agriParams.weather_DHI = weather_data.DHI;       
% agriParams.weather_beta_s = weather_data.beta_s; 
% agriParams.weather_phi_s = weather_data.phi_s;

% PV parameters
agriParams.PV_n_p = 0.2;      % panel efficiency
agriParams.PV_psi = 0;        % field layout angle (rad)
agriParams.PV_startup_period=1; %years spent creating the solar array

% Crop & Econ parameters
agriParams.crop_elec_price = 0.5; 
agriParams.crop_price = 16.03; %USD/kg
agriParams.crop_HI = 0.3; %this is harvest index, for raspberries it is roughly .3, so 30% of the plant weight is berries- changes based on crop choice
agriParams.crop_MC = .85; %this is moisture content, raspberries are about 85% water
agriParams.crop_RUE = 1.0; %can't find specific data for rasperries--> 1.0 is a typical RUE for tomatoes and lettuce
agriParams.crop_k =  0.65; %light extinction coefficient, crop dependent
%LAI for spring, summer, fall, winter
agriParams.crop_LAI = [1.0, 3.0, 1.0, 0.0]; %leaf area index (sq meters of leaves per sq meter of ground, crop dependent

agriParams.agrivoltaic_capital_cost_rate = 1500; %USD per kW in net costs
agriParams.investigation_period = 20;
agriParams.discount_rate = 0.07;

agriParams.year0_cost=402.89;%year 0 planning costs per acre
agriParams.year1_cost=12229.19;%year 1 planting costs per acre
agriParams.year2_cost=2943.4;%year 2 low yield establishment costs per acre
agriParams.startup_years=3;%number of startup years for raspberries prior to steady state production

agriParams.ongoing_OMcost = 18543.62; %ongoing cost of operations and maintenance (post startup years) for raspberry growing per acre

% additions for microclimate
agriParams.crop_PAR_frac = 0.48; % portion of sunlight usable for photosynthesis--> changed to 0.48 from 0.45 (from more recent study)
agriParams.crop_T_base = 5; % all temps in celcius --> below this temp raspberries can't grow
agriParams.crop_T_opt = 20; % this is the optimal temperature for raspberry growth
agriParams.crop_T_max = 30; % about the max temp that raspberries can  grow at
agriParams.crop_c_T = 2.5; % degrees celcius per unit SF--> conflicting data for this--> could be anywhere from 1.5 to 3.7 deg celcius
agriParams.crop_GCF = 0.3; % Ground Cover Fraction --> the portion of the plot that is actually covered my raspberries (calc. from Penn State metrics for open field red raspberries)
% reasonable GCF range --> [0.15-0.5]

% environmental parameters
base_dir = fileparts(mfilename('fullpath'));

if isempty(base_dir) % in case there are issues with referencing the main dir
    base_dir = pwd;
end
data_dir = fullfile(base_dir, 'parameterData');
%for each of the solstice days--> timezone diff relative to UTC--> also
%for each of the solstice days--> timezone diff relative to UTC--> also
%Created representative days for each season --> averaged values for each
%seasonal day across the 90ish days of data gathered for season, spanning
%from December 2024 to November 2025
ci_jan21 = get_season_hourly_ci(data_dir,'representativeDaysMisoMarginalEmissions/representative_winter_et_5min.csv'); % winter
ci_mar21 = get_season_hourly_ci(data_dir,'representativeDaysMisoMarginalEmissions/representative_spring_et_5min.csv'); % spring
ci_jun21 = get_season_hourly_ci(data_dir,'representativeDaysMisoMarginalEmissions/representative_summer_et_5min.csv'); % summer
ci_sep21 = get_season_hourly_ci(data_dir,'representativeDaysMisoMarginalEmissions/representative_fall_et_5min.csv'); % fall
% Combine all seasons in order, starting from the first file passed
% [96 x 1]
ci_all_seasons = [ci_jan21; ci_mar21; ci_jun21; ci_sep21];
% in lbs /MWh --> marginal CO2 emissions which should be about equal to (but a little less than) marginal CO2eq emissions
agriParams.env_ci_marginal_hourly_miso = ci_all_seasons;
clear ci_jan21 ci_sep21 ci_jun21 ci_mar21 ci_all_seasons;



%% Design Variables

% Panel layout variables
agriVar.PV_z_p = 2.5;           % panel height (m)
agriVar.PV_l_p = 2.46;           % panel length (m)
agriVar.PV_w_p = 1.37;           % panel width (m)
agriVar.PV_phi = .28;           % azimuth (rad)
agriVar.PV_sigma = 1.12;      % tilt (rad)
agriVar.PV_y_p = 2.5;           % row distance (m)
agriVar.PV_x_p = 0.1;         % panel distance (m)

%override variables
% agriVar.PV_z_p = 1;
% agriVar.PV_sigma = 0;

%% if single-axis tracking, set bounds + physics-based initialization
if agriParams.tracking_mode == 1
    
    % Create bounds for 96 tracking variables
    tracking_lb = zeros(1, 96); 
    tracking_ub = ones(1, 96) * agriParams.PV_max_tilt; 
    
    % Append to existing bounds
    lb = [lb, tracking_lb]; 
    ub = [ub, tracking_ub]; 

    % Physics-based initialization
    agriVar.tracking_angles = generate_physics_tracking(agriParams, agriVar);

else
    % Required for Simulink bus consistency
    agriVar.tracking_angles = zeros(4,24);
end
%%  Simulink Bus Objects
agriParams = orderfields(agriParams);
agriVar    = orderfields(agriVar);

% Create the bus for parameters
info_1 = Simulink.Bus.createObject(agriParams);
params_bus = eval(info_1.busName); 
clear(info_1.busName); 

% Create the bus for variables
info_2 = Simulink.Bus.createObject(agriVar);
var_bus = eval(info_2.busName);
clear(info_2.busName);

% Clean up
clear info_1 info_2;


%% helper functions
function ci_avg_hourly = get_season_hourly_ci(data_dir, file_season)
    ci_data_season = readtable(fullfile(data_dir, file_season));
    ci_season_day_5min = ci_data_season.mean_marginal_co2_lbs_per_mwh;
    % 12 five-minute points per hour --> 24 hourly vals--> in lbs/MWh
    ci_avg_hourly = mean(reshape(ci_season_day_5min, 12, []), 1).';
end


