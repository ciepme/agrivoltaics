function [SF_spring, SF_summer, SF_fall, SF_winter] = calculate_shading(var, params)
    
    % pv parameters
    PV_length = var.PV_l_p;
    PV_width = var.PV_w_p;
    PV_tilt = rad2deg(var.PV_sigma);    
    PV_azimuth = rad2deg(var.PV_phi);   
    PV_height = var.PV_z_p;
    
    row_distance = var.PV_y_p;
    pair_distance = var.PV_x_p;
    crop_azimuth = PV_azimuth;
    

    unit_width  = PV_width  + pair_distance;  % w_p + x_p
    unit_length = PV_length + row_distance;   % l_p + y_p

    
    n_pv = 1;
    PV_pairs = [unit_width/2, unit_length/2, PV_height];
    
    % 24 hour shading for each representative day
    % FIX 1: Pass params and the specific row of var.tracking_angles down!
    SF_spring = calc_seasonal_shading(params.weather.spring, PV_azimuth, PV_tilt, n_pv, PV_pairs, PV_length, PV_width, unit_width, unit_length, crop_azimuth, params, var.tracking_angles(1,:));
    SF_summer = calc_seasonal_shading(params.weather.summer, PV_azimuth, PV_tilt, n_pv, PV_pairs, PV_length, PV_width, unit_width, unit_length, crop_azimuth, params, var.tracking_angles(2,:));
    SF_fall   = calc_seasonal_shading(params.weather.fall,   PV_azimuth, PV_tilt, n_pv, PV_pairs, PV_length, PV_width, unit_width, unit_length, crop_azimuth, params, var.tracking_angles(3,:));
    SF_winter = calc_seasonal_shading(params.weather.winter, PV_azimuth, PV_tilt, n_pv, PV_pairs, PV_length, PV_width, unit_width, unit_length, crop_azimuth, params, var.tracking_angles(4,:));
end
% 24 hour loops for each season
% FIX 2: Update inputs to catch params and hourly_angles
function SF_season = calc_seasonal_shading(weather_season, PV_azimuth, PV_tilt, n_pv, PV_pairs, PV_length, PV_width, unit_width, unit_length, crop_azimuth, params, hourly_angles)
    
    SF_season = zeros(24, 1);
    
    for t = 1:24
        sun_alt_deg = rad2deg(weather_season.beta_s(t));
        sun_az_deg = rad2deg(weather_season.phi_s(t));
        
        if sun_alt_deg > 0
            
%Toggle logic for single-axis versus fixed axis tracking
if params.tracking_mode == 1
                % SINGLE-AXIS TRACKING
                % Because we fixed the GA bounds, hourly_angles ALREADY contains 
                % negative numbers for morning (East) and positive for afternoon (West). 
                tilt_current = rad2deg(hourly_angles(t)); 
                
                % Lock the farm's torque tube orientation (0 = South)
                phi_current = PV_azimuth; 
            else
                % FIXED AXIS
                tilt_current = PV_tilt;
                phi_current = PV_azimuth;
end

% Pass the tilt_current, phi_current, and params.tracking_mode down:
SF_temp = PV_shading_factor_calc(sun_az_deg, sun_alt_deg, ...
                            phi_current, tilt_current, ...
                            n_pv, PV_pairs, ...
                            PV_length, PV_width, ...
                            unit_width, unit_length, crop_azimuth, params.tracking_mode);
            SF_season(t) = SF_temp;
        else
            SF_season(t) = 0; 
        end
    end
end
%shading math
function SF = PV_shading_factor_calc(sun_azimuth, sun_altitude, PV_azimuth, PV_tilt, n_pv, PV_pairs, PV_length, PV_width, crop_width, crop_length, crop_azimuth, tracking_mode)
   shading_area = repmat(polyshape(), 1, n_pv);
    for pv = 1:n_pv
        % Pass tracking_mode down one more level to the projection
        shading_area(pv) = shade_projection(PV_pairs(pv,:)', PV_width, PV_length, PV_azimuth, PV_tilt, sun_azimuth, sun_altitude, crop_width, crop_length, crop_azimuth, tracking_mode);
    end     
    total_shadow = union(shading_area);
    SF = area(total_shadow) / (crop_width * crop_length);
end
% 3D polygon projection
function shading_box = shade_projection(center_point, PV_width, PV_length, PV_azimuth, PV_tilt, sun_azimuth, sun_altitude, crop_width, crop_length, crop_azimuth, tracking_mode)
    panel_2D = [PV_width/2, PV_length/2, 0; -PV_width/2, PV_length/2, 0; -PV_width/2, -PV_length/2, 0; PV_width/2, -PV_length/2,  0];
    
    % Rz keeps your farm's global orientation locked
    Rz_pv = [-sind(PV_azimuth), -cosd(PV_azimuth), 0; cosd(PV_azimuth), -sind(PV_azimuth), 0; 0, 0, 1];
    Rz_crop = [-sind(crop_azimuth), -cosd(crop_azimuth), 0; cosd(crop_azimuth), -sind(crop_azimuth), 0; 0, 0, 1]; 
    
    if tracking_mode == 1
        % BARREL ROLL: Rotates around the Y-axis (North-South torque tube)
        Ry_pv = [cosd(PV_tilt), 0, sind(PV_tilt); 0, 1, 0; -sind(PV_tilt), 0, cosd(PV_tilt)];
        panel_global_pv = repmat(Rz_crop*center_point, 1, 4) + Rz_pv*Ry_pv*panel_2D';
    else
        % FIXED TILT: Pitches forward around the X-axis
        Rx_pv = [1,0,0; 0, cosd(PV_tilt), -sind(PV_tilt); 0, sind(PV_tilt), cosd(PV_tilt)];
        panel_global_pv = repmat(Rz_crop*center_point, 1, 4) + Rz_pv*Rx_pv*panel_2D';
    end
    ncrop = [0,0,1]; 
    sun_vector = [cosd(sun_azimuth)*cosd(sun_altitude), sind(sun_azimuth)*cosd(sun_altitude), sind(sun_altitude)]; 
    PV_vertex_projected = zeros(4,3);
    for vertex = 1:4
        t = -(ncrop(1)*panel_global_pv(1,vertex) + ncrop(2)*panel_global_pv(2,vertex) + ncrop(3)*panel_global_pv(3,vertex)) / ...
             (ncrop(1)*sun_vector(1) + ncrop(2)*sun_vector(2) + ncrop(3)*sun_vector(3)); 
        PV_vertex_projected(vertex,:) = [panel_global_pv(1,vertex) + sun_vector(1)*t, panel_global_pv(2,vertex) + sun_vector(2)*t, panel_global_pv(3,vertex) + sun_vector(3)*t]; 
    end
    PV_vertex_projected_ENZ = transpose(Rz_crop) * PV_vertex_projected'; 
    crop_shape = polyshape([0, crop_width, crop_width, 0], [0, 0, crop_length, crop_length]);
    shading_box = intersect(polyshape(PV_vertex_projected_ENZ(1,:), PV_vertex_projected_ENZ(2,:)), crop_shape);
end