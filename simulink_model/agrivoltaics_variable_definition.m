%% Add Path
clear
addpath(genpath(pwd));

%% Persona Mode Definition
%Options: 
% None (Default model)
% Fred (Large-scale commodity crops) 
% Janice (small-scale specialtycrop)
persona_mode = 'Fred';

%% Tracking System Definition (Fixed vs Single Axis)
% 0 = Fixed Axis, 1 = Single-Axis Tracking 
% run: bdclose('all'); clear; clear classes; clear functions;
% clc; when switching between fixed axis and single axis modes
agriParams.tracking_mode = 0; 

% Max rotation angle for the single-axis tracker (mechanical limit before
% hitting motor housing or sturcture
agriParams.PV_max_tilt = 60 * (pi/180); % Convert degrees to radians

% Max slew rate for single-axis tracker (radians per hour)
agriParams.max_slew_per_hour = deg2rad(20);

% tracking angles for optimizer matrix definition
% 4x24 matrix (4 seasons, 24 hours). 
agriVar.tracking_angles = zeros(4, 24);


%% Parameter Definition (fixed values)
agriParams.machinery_width=0;
% Persona specific parameters
if strcmp(persona_mode, 'Fred')
    agriParams.persona = 1;
    agriParams.land_x = 2000;       % length of base (m)
    agriParams.land_y = 2000;       % length of height (m)
    % crop parameters for corn
    agriParams.crop_price = .25; % USD/kg, selling price
    agriParams.crop_HI = 0.45; % harvest index, % of the plant weight that is harvestable
    agriParams.crop_MC = .15; % moisture content
    agriParams.crop_RUE = 3.5; % g/MJ PAR, radiation use efficiency
    agriParams.crop_k =  0.65; % light extinction coefficient
    agriParams.crop_T_base = 10; % all temps in C --> below this temp raspberries can't grow
    agriParams.crop_T_opt = 30; % optimal temperature for raspberry growth
    agriParams.crop_T_max = 40; % max temp that raspberries can  grow at
    crop_height = 3.5; % (m)
    min_height = crop_height; % lower bound panel height
    agriParams.crop_canopy_cover = .9; % ground cover of corn
    % LAI for spring, summer, fall, winter
    agriParams.crop_LAI = [1.0, 5.0, 1.5, 0.0]; %leaf area index (sq meters of leaves per sq meter of ground)
    % economic parameters 
    agriParams.agrivoltaic_capital_cost_rate = 1400; %USD per kW in net costs
    agriParams.year0_cost=40; %year 0 planning costs per acre
    agriParams.year1_cost=0; %year 1 planting costs per acre
    agriParams.year2_cost=0; %year 2 low yield establishment costs per acre
    agriParams.startup_years=0; % number of startup years for steady state production
    agriParams.ongoing_OMcost = 800; % ongoing cost of operations and maintenance (post startup years) for raspberry growing per acre
    % farming parameters
    agriParams.machinery_width = 7.5; % (m) for combine header
    agriParams.harvest_clearance = 1; % (m) space needed between edge of crops and panel
    min_row_width = 9.5; % 2*harvest clearance + machinery width
    max_row_width = agriParams.machinery_width*3;
    agriParams.foundation_flag = 2; % for ground screw foundation
    % for yearly electricity demand
    land_area_m2 = agriParams.land_x * agriParams.land_y;
    % baseload + per-square-meter variable load
    base_facility_kWh = 14000; % baseline cold storage + shop/packing overhead
    per_m2_kWh        = 0.42;  % combined irrigation pumping & variable cooling per m²
    agriParams.farm_elec_demand_kWh = base_facility_kWh + (land_area_m2 * per_m2_kWh);
    agriParams.net_metering_mode =  0; %sell back to the grid = 0, self-use = 1
elseif strcmp(persona_mode, 'Janice')
    agriParams.persona = 2;
    agriParams.land_x = 250;       % length of base (m)
    agriParams.land_y = 250;       % length of height (m)
    % crop parameters for raspberries
    agriParams.crop_price = 16.03; % USD/kg
    agriParams.crop_HI = 0.3; % harvest index, % of the plant weight that is harvestable
    agriParams.crop_MC = .85; % moisture content
    agriParams.crop_RUE = 1.0; % can't find specific data for rasperries--> 1.0 is a typical RUE for tomatoes and lettuce
    agriParams.crop_k =  0.65; % light extinction coefficient
    agriParams.crop_T_base = 5; % all temps in C --> below this temp raspberries can't grow
    agriParams.crop_T_opt = 20; % optimal temperature for raspberry growth
    agriParams.crop_T_max = 30; % max temp that raspberries can  grow at
    agriParams.crop_canopy_cover = .6; % ground cover of raspberries, need gaps between rows for sunlight and airflow
    % LAI for spring, summer, fall, winter
    agriParams.crop_LAI = [1.0, 3.0, 1.0, 0.0]; %leaf area index (sq meters of leaves per sq meter of ground)
    % economic parameters 
    agriParams.agrivoltaic_capital_cost_rate = 2000; %USD per kW in net costs
    agriParams.year0_cost=400;%year 0 planning costs per acre
    agriParams.year1_cost=0;% year 1 planting costs per acre
    agriParams.year2_cost=0;%year 2 low yield establishment costs per acre
    agriParams.startup_years=0;%number of startup years for raspberries prior to steady state production
    agriParams.ongoing_OMcost = 18000; %ongoing cost of operations and maintenance (post startup years) for raspberry growing per acre
    % farming parameters
    agriParams.harvest_clearance = .5;  
    min_height = 2; % lower bound panel height, based on human height clearance
    plant_width = 1;
    min_row_width = 2; % 2*harvest clearance + plant width
    max_row_width = 10;
    agriParams.foundation_flag = 1; % for driven pier/concrete foundation
    % for yearly electricity demand
    land_area_m2 = agriParams.land_x * agriParams.land_y;
    % baseload + per-square-meter variable load
    base_facility_kWh = 14000; % baseline cold storage + shop/packing overhead
    per_m2_kWh        = 0.42;  % combined irrigation pumping & variable cooling per m²
    agriParams.farm_elec_demand_kWh = base_facility_kWh + (land_area_m2 * per_m2_kWh);
    agriParams.net_metering_mode =  1; %sell back to the grid = 0, self-use = 1
elseif strcmp(persona_mode, 'None')
    agriParams.persona = 3;
    agriParams.land_x = 50;       % length of base (m)
    agriParams.land_y = 50;       % length of height (m)
    % crop parameters for raspberries
    agriParams.crop_price = 16.03; %USD/kg
    agriParams.crop_HI = 0.3; %this is harvest index, for raspberries it is roughly .3, so 30% of the plant weight is berries- changes based on crop choice
    agriParams.crop_MC = .85; %this is moisture content, raspberries are about 85% water
    agriParams.crop_RUE = 1.0; %can't find specific data for rasperries--> 1.0 is a typical RUE for tomatoes and lettuce
    agriParams.crop_k =  0.65; %light extinction coefficient, crop dependent
    agriParams.crop_LAI = [1.0, 3.0, 1.0, 0.0]; %leaf area index (sq meters of leaves per sq meter of ground)
    agriParams.crop_T_base = 5; % all temps in C --> below this temp raspberries can't grow
    agriParams.crop_T_opt = 20; % optimal temperature for raspberry growth
    agriParams.crop_T_max = 30; % max temp that raspberries can  grow at
    agriParams.crop_canopy_cover = .6;
    % economic parameters
    agriParams.agrivoltaic_capital_cost_rate = 1600; %USD per kW in net costs
    agriParams.year0_cost=400;%year 0 planning costs per acre
    agriParams.year1_cost=12000;% year 1 planting costs per acre
    agriParams.year2_cost=3000;%year 2 low yield establishment costs per acre
    agriParams.startup_years=3;%number of startup years for raspberries prior to steady state production
    agriParams.ongoing_OMcost = 18000; %ongoing cost of operations and maintenance (post startup years) for raspberry growing per acre
    % farming parameters
    agriParams.harvest_clearance = .5;  
    min_height = 2; % lower bound panel height, based on human height clearance
    plant_width = 1;
    min_row_width = 2; % 2*harvest clearance + plant width
    max_row_width = 10;
    agriParams.foundation_flag = 1; % for driven pier/concrete foundation
    % for yearly electricity demand
    land_area_m2 = agriParams.land_x * agriParams.land_y;
    % baseload + per-square-meter variable load
    base_facility_kWh = 14000; % baseline cold storage + shop/packing overhead
    per_m2_kWh        = 0.42;  % combined irrigation pumping & variable cooling per m²
    agriParams.farm_elec_demand_kWh = base_facility_kWh + (land_area_m2 * per_m2_kWh);
    agriParams.net_metering_mode =  0; %sell back to the grid = 0, self-use = 1
else 
    disp('error:input persona mode')
    agriParams.land_x = 0;       % length of base (m)
    agriParams.land_y = 0;       % length of height (m)
    agriParams.land_angle = 0;    % rotation (rad)
end

% Weather parameter loading (currently set for Ann Arbor, MI)
load('pv_weather_4seasons.mat', 'weather_struct');
agriParams.weather = weather_struct;

% PV parameters
agriParams.PV_n_p = 0.2;      % panel efficiency
agriParams.PV_psi = 0;        % field layout angle (rad)
agriParams.PV_startup_period=1; % years spent building the solar array

% Economic parameters
agriParams.elec_sell_price = 0.054; %selling price for electricity to the grid
agriParams.elec_buy_price = 0.21; %buying price from utility
agriParams.investigation_period = 30;
agriParams.discount_rate = 0.07;

% Microclimate parameters
agriParams.crop_PAR_frac = 0.48; % portion of sunlight usable for photosynthesis
agriParams.crop_c_T = 2.5; % degrees C per unit SF

% Environmental parameters
% specific for MISO area grid emissions
base_dir = fileparts(mfilename('fullpath'));
if isempty(base_dir) % in case there are issues with referencing the main dir
    base_dir = pwd;
end
%data_dir = fullfile(base_dir, 'parameterData');
%for each of the solstice days--> timezone diff relative to UTC
%Created representative days for each season --> averaged values for each
%seasonal day across the 90ish days of data gathered for season, spanning
%from December 2024 to November 2025
ci_jan21 = get_season_hourly_ci('representative_winter_et_5min.csv'); % winter
ci_mar21 = get_season_hourly_ci('representative_spring_et_5min.csv'); % spring
ci_jun21 = get_season_hourly_ci('representative_summer_et_5min.csv'); % summer
ci_sep21 = get_season_hourly_ci('representative_fall_et_5min.csv'); % fall
% Combine all seasons in order, starting from the first file passed
% [96 x 1]
ci_all_seasons = [ci_jan21; ci_mar21; ci_jun21; ci_sep21];
% in lbs /MWh --> marginal CO2 emissions which should be about equal to (but a little less than) marginal CO2eq emissions
agriParams.env_ci_marginal_hourly_miso = ci_all_seasons;
clear ci_jan21 ci_sep21 ci_jun21 ci_mar21 ci_all_seasons;

% Hardware parameters 
agriParams.support_width = 0.5;      % Width of the steel solar mount base

%% Design Variables

% Variable Bound Definition
% (1) = PV_z_p; panel height above the ground (m)
% (2) = PV_l_p; panel length, goes in y_p direction (m)
% (3) = PV_w_p; panel width, goes in x_p direction along th (m)
% (4) = PV_phi; azimuth angle (radians) - farm/row orientation relative to true South, going ccw e.g. pi/2 rad is East
% (5) = PV_sigma; tilt angle (radians) - fixed sloping angle of PV relative to horizontal (xy) plane
% (6) = PV_y_p; Distance between rows (m) of panels
% (7) = PV_x_p; Distance between panels in a pair/row (m)
% Order: [Height, Length, Width, Azimuth, Tilt, Row Gap, Panel Gap]
if agriParams.tracking_mode == 1; 
    min_phi = pi/4;
    max_phi = 3*pi/4;
    start_phi = pi/2;
else 
    min_phi = -pi/2;
    max_phi = pi/2; 
    start_phi = 0;
end

lb = [min_height, 1.0, 1.0, min_phi, 0, min_row_width,  0]; 
ub = [5.5, 2.5, 1.0,  max_phi, pi/2, max_row_width, 1.0];

% Panel layout variables
agriVar.PV_z_p = 5;           % panel height (m)
agriVar.PV_l_p = 2;           % panel length (m)
agriVar.PV_w_p = 1.45;        % panel width (m)
agriVar.PV_phi = start_phi;   % azimuth (rad)
agriVar.PV_sigma = 0;         % tilt (rad)
agriVar.PV_y_p = min_row_width;         % row distance (m)
agriVar.PV_x_p = 0.1;         % panel distance (m)


%% Single-axis tracking bounds and initialization
if agriParams.tracking_mode == 1

    % Physics-based initial tracking curve
    agriVar.tracking_angles = generate_physics_tracking(agriParams, agriVar);

    % Clamp initial physics curve to mechanical tilt limit
    agriVar.tracking_angles = max(agriVar.tracking_angles, -agriParams.PV_max_tilt);
    agriVar.tracking_angles = min(agriVar.tracking_angles,  agriParams.PV_max_tilt);

    % Initialize 96 tracking bounds
    tracking_lb = zeros(1, 96);
    tracking_ub = zeros(1, 96);

    seasons = {'spring', 'summer', 'fall', 'winter'};

    max_tilt = agriParams.PV_max_tilt;
    max_slew = agriParams.max_slew_per_hour;

    % Number of hours needed to move from 0 to max tilt
    ramp_hours = ceil(max_tilt / max_slew);

    % Store masks for debugging/plotting
    agriParams.tracking_daytime_mask = false(4, 24);
    agriParams.tracking_move_mask    = false(4, 24);

    for s = 1:4

        beta_s = agriParams.weather.(seasons{s}).beta_s;
        is_daytime = beta_s > 0;

        day_idx = find(is_daytime);

        % Default: tracker must be stowed at 0
        can_move = false(1, 24);

        if ~isempty(day_idx)

            first_day_hour = day_idx(1);
            last_day_hour  = day_idx(end);

            % Allow motion before sunrise and after sunset
            move_start = max(1,  first_day_hour - ramp_hours);
            move_end   = min(24, last_day_hour  + ramp_hours);

            can_move(move_start:move_end) = true;

        end

        local_start = (s-1)*24 + 1;
        local_end   = s*24;

        tracking_lb(local_start:local_end) = -max_tilt .* can_move;
        tracking_ub(local_start:local_end) =  max_tilt .* can_move;

        agriParams.tracking_daytime_mask(s, :) = is_daytime;
        agriParams.tracking_move_mask(s, :)    = can_move;

    end

    % Append tracking bounds to seven layout-variable bounds
    lb = [lb(1:7), tracking_lb];
    ub = [ub(1:7), tracking_ub];

else

    % Fixed-axis case
    agriVar.tracking_angles = zeros(4, 24);

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
function ci_avg_hourly = get_season_hourly_ci(file_season)
    ci_data_season = readtable(fullfile(file_season));
    ci_season_day_5min = ci_data_season.mean_marginal_co2_lbs_per_mwh;
    % 12 five-minute points per hour --> 24 hourly vals--> in lbs/MWh
    ci_avg_hourly = mean(reshape(ci_season_day_5min, 12, []), 1).';
end


