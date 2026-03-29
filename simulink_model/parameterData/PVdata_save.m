
weather_data = readtable('Weather_Data_AA.csv', 'HeaderLines', 2);
spa_data = readtable('Spa_Data_AA.csv'); 

%parse SPA dats
spa_dates = datetime(spa_data{:, 1});
spa_months = month(spa_dates);
spa_days = day(spa_dates);

%target days of equinox and solstice
target_days = [3, 21;   % Spring
               6, 21;   % Summer
               9, 21;   % Fall
              12, 21];  % Winter

season_names = {'spring', 'summer', 'fall', 'winter'};
weather_struct = struct(); % Initialize empty structure

%4 season loop
for i = 1:4
    m = target_days(i, 1);
    d = target_days(i, 2);
    curr_season = season_names{i};
    
    % --- MATCH TMY DATA ---
    tmy_idx = (weather_data.Month == m) & (weather_data.Day == d);
    day_weather = weather_data(tmy_idx, :);
    
    % --- MATCH SPA DATA ---
    spa_idx = (spa_months == m) & (spa_days == d);
    day_spa = spa_data(spa_idx, :); 
    
    % Pack the Weather Data
    weather_struct.(curr_season).DNI = day_weather.DNI;
    weather_struct.(curr_season).DHI = day_weather.DHI;
    weather_struct.(curr_season).T_air = day_weather.Temperature;
    weather_struct.(curr_season).Wind = day_weather.WindSpeed;
    
    % Pack the Geometry Data (Col 3 = Zenith, Col 4 = Azimuth)
    zenith_deg = day_spa{:, 3}; 
    spa_azimuth_deg = day_spa{:, 4}; 
    
    % Convert Zenith to Altitude (beta_s) in radians
    weather_struct.(curr_season).beta_s = deg2rad(90 - zenith_deg);
    
    % Flipped by -1 to make East positive (Counter-Clockwise)
    weather_struct.(curr_season).phi_s = deg2rad(-1 * spa_azimuth_deg);
end

%.mat file
save('pv_weather_4seasons.mat', 'weather_struct');
