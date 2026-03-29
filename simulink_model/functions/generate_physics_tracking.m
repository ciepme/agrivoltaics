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

            % --- MATCH SIMULINK LOGIC ---
            % East in morning, West in afternoon
            if phi_s(h) > 0
                phi_current = pi/2;
            else
                phi_current = -pi/2;
            end

            % Solve for optimal tilt (sigma) for this fixed azimuth
            % Derived from maximizing cos(theta)
            sigma_opt = atan( ...
                cos(beta_s(h)) * cos(phi_s(h) - phi_current) / ...
                sin(beta_s(h)) ...
            );

            % --- Smooth saturation (gradient-friendly) ---
            sigma_opt = max_tilt * tanh(sigma_opt / max_tilt);

            theta_day(h) = sigma_opt;
        end

        % --- Smooth across time ---
        theta_day = smoothdata(theta_day, 'gaussian', 3);

        tracking_angles(d,:) = theta_day;
    end
end