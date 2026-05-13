weather_data = readtable('Weather_Data_AA.csv', 'HeaderLines', 2);
spa_data = readtable('Spa_Data_AA.csv'); 

% Parse SPA dates
spa_dates = datetime(spa_data{:, 1});
spa_months = month(spa_dates);
spa_days = day(spa_dates);

% Define seasons by their constituent months
seasons.spring = [3, 4, 5];
seasons.summer = [6, 7, 8];
seasons.fall   = [9, 10, 11];
seasons.winter = [12, 1, 2];

season_names = {'spring', 'summer', 'fall', 'winter'};

% Target days of equinox and solstice (used ONLY for representative geometry)
target_days = [3, 21;   % Spring
               6, 21;   % Summer
               9, 21;   % Fall
              12, 21];  % Winter

weather_struct = struct(); % Initialize empty structure

% 4 season loop
for i = 1:4
    curr_season = season_names{i};
    season_months = seasons.(curr_season);
    
    % --- average TMY weather data ---
    % rows that fall within the current season's months
    season_idx = ismember(weather_data.Month, season_months);
    season_weather = weather_data(season_idx, :);
    
    % Group by Hour (and Minute) and calculate the mean for the weather variables
    % This creates an "average day" profile for the season
    avg_weather = groupsummary(season_weather, {'Hour', 'Minute'}, 'mean', ...
        {'DNI', 'DHI', 'Temperature', 'WindSpeed'});
        
    % Sort to ensure the data is strictly chronological from 0:00 to 23:00
    avg_weather = sortrows(avg_weather, {'Hour', 'Minute'});
    
    % Pack the Averaged Weather Data
    weather_struct.(curr_season).DNI = avg_weather.mean_DNI;
    weather_struct.(curr_season).DHI = avg_weather.mean_DHI;
    weather_struct.(curr_season).T_air = avg_weather.mean_Temperature;
    weather_struct.(curr_season).Wind = avg_weather.mean_WindSpeed;
    
    % --- 2. MATCH SPA DATA (Using the 21st as the representative day) ---
    m = target_days(i, 1);
    d = target_days(i, 2);
    
    spa_idx = (spa_months == m) & (spa_days == d);
    day_spa = spa_data(spa_idx, :); 
    
    % Pack the Geometry Data (Col 3 = Zenith, Col 4 = Azimuth)
    zenith_deg = day_spa{:, 3}; 
    spa_azimuth_deg = day_spa{:, 4}; 
    
    % Convert Zenith to Altitude (beta_s) in radians
    weather_struct.(curr_season).beta_s = deg2rad(90 - zenith_deg);
    
    % Flipped by -1 to make East positive (Counter-Clockwise)
    weather_struct.(curr_season).phi_s = deg2rad(-1 * spa_azimuth_deg);
end

% Save to .mat file
save('pv_weather_4seasons.mat', 'weather_struct');