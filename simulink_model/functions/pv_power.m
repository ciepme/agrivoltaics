function [P_annual, P_spring_hourly, P_summer_hourly, P_fall_hourly, P_winter_hourly, n_cols, n_rows, total_panels] = pv_power(var, params)
    %PV model calculating incidence angle, total radiation, and estimated power
    %output
        % DNI ; %direct normal irradiance
        % DHI ; % diffuse horizontal irradiance
        % beta_s ; %sun altitude angle (angle of sun above the local horizon)
        % phi_s ; %sun azimuth angle (angle between the sun and true south)
        l_p = var.PV_l_p; % panel length (m)
        w_p = var.PV_w_p; % panel width (m)
        phi = var.PV_phi; % azimuth angle (radians) - relative to true South, going ccw e.g. pi/2 rad is East
        sigma =var.PV_sigma; % tilt angle (radians) - fixed sloping angle of PV relative to horizontal (xy) plane 
        n_p = params.PV_n_p;%panel efficiency
        x_p = var.PV_x_p; % Distance between panels in a pair/row (m)
        y_p = var.PV_y_p; % Distance between rows (m)
        x = params.land_x; % Total plot width (m)
        y = params.land_y; % Total plot length (m)
    
    %number of panels that fit on the plot
if isfield(params, 'geometry_mode') && params.geometry_mode == 1
    n_cols = params.slice_count;
    n_rows = params.row_count;
elseif params.tracking_mode == 1
    % TRACKING LAYOUT: Panels roll East/West. 
    % They need max clearance when flat (width = w_p). 
    % The length (l_p) does not tilt, so it is just l_p.
    n_cols = max(1, floor(x / (w_p + x_p)));
    n_rows = max(1, floor(y / (l_p + y_p))); 
else
    % FIXED TILT LAYOUT: Panels pitch South. 
    % The length footprint is reduced by the cosine of the tilt.
    n_cols = max(1, floor(x / (w_p + x_p)));
    n_rows = max(1, floor(y / (l_p*cos(sigma) + y_p)));
end
    A_p = l_p*w_p ; %panel area m^2
    
    if n_cols == Inf
        n_cols = 0;
    end
    if n_rows == Inf
        n_rows = 0;
    end
    total_panels = n_rows * n_cols;
    
    [P_spring_daily, P_spring_hourly] = calc_daily_pv(params.weather.spring, phi, sigma, n_p, A_p, total_panels, params, var.tracking_angles(1,:));
    [P_summer_daily, P_summer_hourly] = calc_daily_pv(params.weather.summer, phi, sigma, n_p, A_p, total_panels, params, var.tracking_angles(2,:));
    [P_fall_daily, P_fall_hourly]     = calc_daily_pv(params.weather.fall, phi, sigma, n_p, A_p, total_panels, params, var.tracking_angles(3,:));
    [P_winter_daily, P_winter_hourly] = calc_daily_pv(params.weather.winter, phi, sigma, n_p, A_p, total_panels, params, var.tracking_angles(4,:));
    days_per_season = 365 / 4;
    P_annual = (P_spring_daily + P_summer_daily + P_fall_daily + P_winter_daily) * days_per_season;
    
end
%
function [P_daily, P_hourly] = calc_daily_pv(season_weather, phi, sigma, n_p, A_p, total_panels, params, hourly_angles)
 
    p_tot=0;
    P_hourly = zeros(24,1);
    DNI = season_weather.DNI;
    DHI = season_weather.DHI;
    beta_s = season_weather.beta_s;
    phi_s = season_weather.phi_s;
    for i = 1:24
        if beta_s(i) <= 0
            continue; % Skip to next hour (P_inst remains 0 effectively)
        end
%Toggle logic for single-axis versus fixed axis tracking
        if params.tracking_mode == 1
            % SINGLE-AXIS TRACKING (Optimizer Driven)
            raw_angle = hourly_angles(i); 
            
            % The physical tilt of the panel is the absolute magnitude of the roll
            sigma_current = abs(raw_angle);
            
            % The tracker axis faces South, so as it rolls, the panel surface 
            % faces perfectly East (+pi/2) or perfectly West (-pi/2).
            if raw_angle < 0
                phi_current = pi/2;  % Rolled East
            elseif raw_angle > 0
                phi_current = -pi/2; % Rolled West
            else
                phi_current = phi;   % Flat at solar noon
            end
            
        else
            % FIXED AXIS
            % Ignore the matrix, just use the static variables
            sigma_current = sigma;
            phi_current = phi;
        end
    % Incidence angle
    % FIX: Now uses sigma_current and phi_current to calculate actual intercepted light
    cos_theta = cos(beta_s(i))*cos(phi_s(i) - phi_current)*sin(sigma_current) + sin(beta_s(i))*cos(sigma_current);
    
    % Prevent negative values if panel is facing away from sun
        cos_theta = max(0, cos_theta);
    
    % Direct beam
    I_db = DNI(i) .* cos_theta;
    
    % Diffuse irradiance
    I_diff = DHI(i) .* (1 + cos(sigma_current)) / 2;
    
    % Reflected (Need to double check equation)
    %I_ref = (DNI*sin(beta_s) + DHI) .* albedo .* (1 - cos(tilt)) / 2;
    
    % Total Intercepted Radiation
    R = I_db + I_diff; %+ I_ref;
    
    % Simple power model
    
    p_hour = n_p * A_p * R;          % hourly power per panel
    P_hourly(i) = p_hour*total_panels; % hourly power for entire plot
    p_tot = p_tot + P_hourly(i); % in kWh
    
    end
    P_daily=p_tot;
    if isnan(P_daily)
        P_daily = 0;
    end
    P_hourly(isnan(P_hourly)) = 0; 
end
