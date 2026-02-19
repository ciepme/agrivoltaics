%% Import Excel
excelData = readtable('DHIDNI_Data.xls');

DNI = excelData.DNI;
DHI = excelData.DHI;

Topocentric_zenith_angle = [128.98; 129.30; 126.01; 119.74; 111.34; 101.57; 90.97; ...
                            79.88; 68.87; 58.16; 48.29; 40.16; 35.18; 34.80; 39.17; ...
                            46.92; 56.60; 67.24; 78.26; 88.98; 100.10; 110.07; 118.78; 125.49];

Topocentric_azimuth_angle = [168.86; 188.04; 206.37; 222.44; 236.04; 247.73; 258.22; ...
                             268.25; 278.57; 290.08; 303.99; 321.92; 344.96; 11.00; 34.75; ...
                             53.42; 67.83; 79.60; 90.02; 100.03; 110.39; 121.82; 135.08; 150.77];

% 2. Calculate beta_s (Altitude) in Radians
% Logic: 90 degrees - Zenith
beta_s = deg2rad(90 - Topocentric_zenith_angle);

% 3. Calculate phi_s (Azimuth) in Radians
% Logic: Negate the NREL value to switch from CW (Westward) to CCW (Eastward)
phi_s = deg2rad(-Topocentric_azimuth_angle);

%set hourly vector
time_vector = [1:24]';



%% Save .mat file
save('pv_inputs.mat', ...
     'DNI','DHI','beta_s','phi_s');
