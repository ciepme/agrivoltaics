model_name = 'agrivoltaics_v1'; 

%% Parameter Definition
% land parameters
% land area parameters
    %  assuming rectangular plot
land_x = 50; % length of base of plot (m)
land_y = 50; % length of height of plot (m)
land_angle = 0; % rotation of plot along z axis (m) - ccw, 0 rad would assume x runs parallel along latitude and y along longitude

params.land = [land_x; land_y; land_angle];
% weather parameters

params.weather = [];
% PV parameters
n_p = .2; % panel efficiency

params.PV = [n_p];

% crop parameters
 electricity_price = .5; % $/kWh
 crop_price = .5; % $/kg 

 params.crop = [electricity_price; crop_price];

 % sustainability parameters

 params.sus = [];

 % economic parameters

 params.econ = [];

 % all parameters
parameters = [params.land;
   params.weather;
   params.PV;
   params.crop;
   params.sus;
   params.econ];

 %% Design Variables
% includes initial guess for design variable
% panel layout variables
z_p = 2; % panel height (m) 
l_p = 1; % panel length (m)
w_p = 1; % panel width (m)
phi = 0; % azimuth angle (radians) - relative to true South, going ccw e.g. pi/2 rad is East
theta = pi/4; % tilt angle (radians) - fixed sloping angle of PV relative to horizontal (xy) plane
psi = 0; % field layout angle (radians) -  direction of main axis of PV array rows measured from west-east line 
y_p = 1; %row distance (m) - distance between parallel rows
x_p = .1; %panel distance (m) - distance between panels within a row

var.PV = [z_p; l_p; w_p; phi; theta; psi; y_p;x_p];

% crop layout variables

var.crop = [];

% all design variables

variables = [var.PV; var.crop];

