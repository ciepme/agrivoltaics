function annual_biomass = calculate_crop_growth(SF_spring, SF_summer, SF_fall, SF_winter, agriParams)
   %  General Crop Parameters
  RUE = agriParams.crop_RUE;
  k = agriParams.crop_k;
  T_base = agriParams.crop_T_base;
  T_opt = agriParams.crop_T_opt;
  T_max = agriParams.crop_T_max;
  c_T  = agriParams.crop_c_T;
  PAR_frac = agriParams.crop_PAR_frac;
  days_per_season = 91.25;
   %  Seasonal Leaf Area Index (LAI)
  LAI_spring = agriParams.crop_LAI(1);
  LAI_summer = agriParams.crop_LAI(2);
  LAI_fall   = agriParams.crop_LAI(3);
  LAI_winter = agriParams.crop_LAI(4);
   %  Calculate Daily Biomass for each representative day
  % calls daily function from the bottom
  daily_spring = calc_daily_growth(agriParams.weather.spring, SF_spring, LAI_spring, RUE, k, c_T, T_base, T_opt, T_max, PAR_frac);
  daily_summer = calc_daily_growth(agriParams.weather.summer, SF_summer, LAI_summer, RUE, k, c_T, T_base, T_opt, T_max, PAR_frac);
  daily_fall   = calc_daily_growth(agriParams.weather.fall, SF_fall, LAI_fall, RUE, k, c_T, T_base, T_opt, T_max, PAR_frac);
  daily_winter = calc_daily_growth(agriParams.weather.winter, SF_winter, LAI_winter, RUE, k, c_T, T_base, T_opt, T_max, PAR_frac);
   % 4. Scale up to the full year
  annual_biomass = (daily_spring * days_per_season) + ...
                   (daily_summer * days_per_season) + ...
                   (daily_fall * days_per_season) + ...
                   (daily_winter * days_per_season);
end
% original daily growth calculation
function daily_biomass = calc_daily_growth(weather_season, SF_season, LAI, RUE, k, c_T, T_base, T_opt, T_max, PAR_frac)
   daily_biomass = 0;
   % Loop through the 24 hours of this specific season
  for t = 1:24
      DNI = weather_season.DNI(t);     
      DHI = weather_season.DHI(t);     
      beta_s = weather_season.beta_s(t);
      T_air  = weather_season.T_air(t);
      SF = SF_season(t);               
    
      if beta_s > 0 % If the sun is up
          R_unshaded = DNI * sin(beta_s) + DHI;
          R_shaded = DHI;
        
          % Blend based on Shading Factor
          R_actual = (R_unshaded * (1 - SF)) + (R_shaded * SF);
        
          % Changing: convert to PAR in MJ/m²/hr
          PAR_MJ = R_actual * 0.0036 * PAR_frac;
          % Beer-Lambert interception eq
          APAR = PAR_MJ * (1 - exp(-k * LAI)); % absorbed PAR
        
          % adding microclimate temp--> depends on the shading factor
          % T_crop = T_air - c_T * 0; % testing shading having no impact
          % on microclimate temp
          T_crop = T_air - c_T * SF;
          % ADDING TEMP RESPONSE TO RUE: temperature response on RUE
          % f_T is basically photosynthesizing efficiency (to be multiplied with RUE)--> max (1.0),
          % when T_crop = 20
          f_T = compute_fT(T_crop, T_base, T_opt, T_max);
          % Accumulate hourly biomass
          daily_biomass = daily_biomass + (RUE * f_T * APAR);
          %this replaces the previous version of this equation:
          % % % The EPIC Crop Equation
          % % daily_biomass = RUE * PAR * (1 - exp(-k * LAI));
      end
  end
    % Failsafe for mathematical edge cases
  if isnan(daily_biomass)
      daily_biomass = 0;
  end
end
function f = compute_fT(T, T_base, T_opt, T_max)
  if T < T_base || T > T_max
      f = 0;
  elseif T <= T_opt
      f = (T - T_base) / (T_opt - T_base);
  else
      f = (T_max - T) / (T_max - T_opt);
  end
end
