model_name = 'agrivoltaics_v1'; 

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
params.crop.crop_price = 0.5;

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
