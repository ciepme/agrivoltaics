function biomass_cropgrowth = calculate_crop_growth(SF_spring, SF_summer, SF_fall, SF_winter, agriParams)
    
    %  General Crop Parameters
    RUE = 2.0;       % Radiation Use Efficiency 
    k = 0.65;        % Light extinction coefficient 
    days_per_season = 91.25; 
    
    %  Seasonal Leaf Area Index (LAI)
    LAI_spring = 1.0; 
    LAI_summer = 3.0; 
    LAI_fall   = 1.0; 
    LAI_winter = 0.0; % Dormant in winter
    
    %  Calculate Daily Biomass for each representative day
    % calls daily function from the bottom
    daily_spring = calc_daily_growth(agriParams.weather.spring, SF_spring, LAI_spring, RUE, k);
    daily_summer = calc_daily_growth(agriParams.weather.summer, SF_summer, LAI_summer, RUE, k);
    daily_fall   = calc_daily_growth(agriParams.weather.fall, SF_fall, LAI_fall, RUE, k);
    daily_winter = calc_daily_growth(agriParams.weather.winter, SF_winter, LAI_winter, RUE, k);
    
    % 4. Scale up to the full year
    biomass_cropgrowth = (daily_spring * days_per_season) + ...
                     (daily_summer * days_per_season) + ...
                     (daily_fall * days_per_season) + ...
                     (daily_winter * days_per_season);
end

% original daily growth calculation
function daily_biomass = calc_daily_growth(weather_season, SF_season, LAI, RUE, k)
    
    total_daily_radiation_W = 0; 
    
    % Loop through the 24 hours of this specific season
    for t = 1:24
        DNI = weather_season.DNI(t);       
        DHI = weather_season.DHI(t);       
        beta_s = weather_season.beta_s(t); 
        SF = SF_season(t);                 
        
        if beta_s > 0 % If the sun is up
            R_unshaded = DNI * sin(beta_s) + DHI;
            R_shaded = DHI;
            
            % Blend based on Shading Factor
            R_actual = (R_unshaded * (1 - SF)) + (R_shaded * SF);
            total_daily_radiation_W = total_daily_radiation_W + R_actual;
        end
    end
    
    % Convert to Megajoules, then to PAR
    total_daily_radiation_MJ = total_daily_radiation_W * 0.0036;
    PAR = total_daily_radiation_MJ * 0.45;
    
    % The EPIC Crop Equation
    daily_biomass = RUE * PAR * (1 - exp(-k * LAI));
    
    % Failsafe for mathematical edge cases
    if isnan(daily_biomass)
        daily_biomass = 0;
    end
end