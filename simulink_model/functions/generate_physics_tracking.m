function tracking_angles = generate_physics_tracking(agriParams, agriVar)
    max_tilt = agriParams.PV_max_tilt;
    weather  = agriParams.weather;
    tracking_angles = zeros(4,24);
    
    seasons = {'spring','summer','fall','winter'};
    
    for d = 1:4
        season_weather = weather.(seasons{d});
        beta_s = season_weather.beta_s; % [1x24]
        phi_s  = season_weather.phi_s;  % [1x24]
        theta_day = zeros(1,24);
        
        for h = 1:24
            if beta_s(h) <= 0
                theta_day(h) = 0;
                continue;
            end
            
            % Bulletproof Morning/Afternoon logic based strictly on the HOUR
            if h < 13
                % Morning: Sun is in the East
                phi_current = pi/2;
            else
                % Afternoon: Sun is in the West
                phi_current = -pi/2;
            end
            
            % Solve for absolute optimal tilt magnitude (sigma) 
            sigma_opt = atan( ...
                cos(beta_s(h)) * cos(phi_s(h) - phi_current) / ...
                sin(beta_s(h)) ...
            );
            
            % --- Smooth saturation (gradient-friendly) ---
            sigma_opt = max_tilt * tanh(sigma_opt / max_tilt);
            
            % FORCE THE SIGN based on the hour!
            % Negative = East (Morning), Positive = West (Afternoon)
            if h < 13
                theta_day(h) = -abs(sigma_opt); % Tip negative to roll East
            else
                theta_day(h) =  abs(sigma_opt); % Tip positive to roll West
            end
        end
        
        % --- Smooth across time ---
        theta_day = smoothdata(theta_day, 'gaussian', 3);
        tracking_angles(d,:) = theta_day;
    end
end