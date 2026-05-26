function [SF_spring, SF_summer, SF_fall, SF_winter] = calculate_shading(var, params)
    if isfield(params, 'geometry_mode') && params.geometry_mode == 1
        [SF_spring, SF_summer, SF_fall, SF_winter] = calculate_row_centered_shading(var, params);
        return;
    end

    
    % pv parameters
    PV_length = var.PV_l_p;
    PV_width = var.PV_w_p;
    PV_tilt = rad2deg(var.PV_sigma);    
    PV_azimuth = rad2deg(var.PV_phi);   
    PV_height = var.PV_z_p;
    
    row_distance = var.PV_y_p;
    pair_distance = var.PV_x_p;
    crop_azimuth = rad2deg(params.land_angle);
    
    % unit cell definition (speeds up calc)
    if params.tracking_mode == 1
        % Legacy single-axis rows run north-south along PV_l_p. The
        % cross-row gap is PV_y_p; the along-row panel gap is PV_x_p.
        unit_width = PV_width + row_distance;
        unit_length = PV_length + pair_distance;
    else
        unit_width = PV_width + pair_distance;
        unit_length = PV_length + row_distance;
    end
    
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

function [SF_spring, SF_summer, SF_fall, SF_winter] = calculate_row_centered_shading(var, params)
    geom = build_row_centered_geometry(var, params);
    PV_tilt = rad2deg(var.PV_sigma);

    SF_spring = calc_row_centered_season(params.weather.spring, PV_tilt, geom, params, var.tracking_angles(1, :));
    SF_summer = calc_row_centered_season(params.weather.summer, PV_tilt, geom, params, var.tracking_angles(2, :));
    SF_fall   = calc_row_centered_season(params.weather.fall,   PV_tilt, geom, params, var.tracking_angles(3, :));
    SF_winter = calc_row_centered_season(params.weather.winter, PV_tilt, geom, params, var.tracking_angles(4, :));
end

function SF_season = calc_row_centered_season(weather_season, PV_tilt, geom, params, hourly_angles)
    sample_heights = get_shading_sample_heights(params);
    SF_by_height = zeros(24, numel(sample_heights));

    for t = 1:24
        sun_alt_deg = rad2deg(weather_season.beta_s(t));
        sun_az_deg = wrap_to_180(rad2deg(weather_season.phi_s(t)));

        if sun_alt_deg <= 0
            continue;
        end

        if params.tracking_mode == 1
            max_tilt = params.PV_max_tilt;
            tilt_current = rad2deg(max(-max_tilt, min(hourly_angles(t), max_tilt)));
        else
            tilt_current = PV_tilt;
        end

        for h = 1:numel(sample_heights)
            SF_by_height(t, h) = row_centered_shading_factor( ...
                sun_az_deg, sun_alt_deg, tilt_current, geom, sample_heights(h));
        end
    end

    SF_season = mean(SF_by_height, 2);
    SF_season = max(0, min(SF_season, 1));
end

function SF = row_centered_shading_factor(sun_azimuth, sun_altitude, PV_tilt, geom, plane_height)
    central_idx = find(abs(geom.panel_offsets(:, 1)) < 1e-12 & ...
        abs(geom.panel_offsets(:, 2)) < 1e-12, 1);
    if isempty(central_idx)
        error('Row-centered geometry is missing a central source panel.');
    end

    base_shadow = row_centered_shadow_projection(geom.panel_centers(central_idx, :)', ...
        geom.panel_span, geom.panel_length, PV_tilt, sun_azimuth, sun_altitude, ...
        geom.crop_azimuth_deg, plane_height);
    [base_x, base_y] = boundary(base_shadow);
    if isempty(base_x)
        SF = 0;
        return;
    end

    hedge_box = poly_bounds(geom.hedge_shape);
    hedge_shadows = repmat(polyshape(), 1, size(geom.panel_centers, 1));
    active_count = 0;

    for pv = 1:size(geom.panel_centers, 1)
        dx = geom.panel_centers(pv, 1) - geom.panel_centers(central_idx, 1);
        dy = geom.panel_centers(pv, 2) - geom.panel_centers(central_idx, 2);
        translated_box = [
            min(base_x) + dx, max(base_x) + dx, ...
            min(base_y) + dy, max(base_y) + dy
        ];

        if ~boxes_overlap(translated_box, hedge_box)
            continue;
        end

        shadow = intersect(polyshape(base_x + dx, base_y + dy), geom.hedge_shape);
        if area(shadow) <= 1e-12
            continue;
        end

        active_count = active_count + 1;
        hedge_shadows(active_count) = shadow;
    end

    if active_count == 0 || geom.hedge_area <= 0
        SF = 0;
        return;
    end

    total_shadow = union(hedge_shadows(1:active_count));
    SF = area(total_shadow) / geom.hedge_area;
    SF = max(0, min(SF, 1));
end

function shadow_poly = row_centered_shadow_projection(center_point, panel_span, panel_length, PV_tilt, sun_azimuth, sun_altitude, crop_azimuth, plane_height)
    if sind(sun_altitude) <= 0
        shadow_poly = polyshape();
        return;
    end

    panel_2D = [
        panel_span / 2,  panel_length / 2, 0
       -panel_span / 2,  panel_length / 2, 0
       -panel_span / 2, -panel_length / 2, 0
        panel_span / 2, -panel_length / 2, 0
    ];

    Rz_crop = [-sind(crop_azimuth), -cosd(crop_azimuth), 0; cosd(crop_azimuth), -sind(crop_azimuth), 0; 0, 0, 1];
    Ry_pv = [cosd(PV_tilt), 0, sind(PV_tilt); 0, 1, 0; -sind(PV_tilt), 0, cosd(PV_tilt)];
    panel_row = repmat(center_point, 1, 4) + Ry_pv * panel_2D';
    panel_global = Rz_crop * panel_row;
    sun_vector = [cosd(sun_azimuth) * cosd(sun_altitude), sind(sun_azimuth) * cosd(sun_altitude), sind(sun_altitude)];
    projected = zeros(4, 3);

    for vertex = 1:4
        t = (plane_height - panel_global(3, vertex)) / sun_vector(3);
        projected(vertex, :) = [
            panel_global(1, vertex) + sun_vector(1) * t, ...
            panel_global(2, vertex) + sun_vector(2) * t, ...
            plane_height
        ];
    end

    projected_ENZ = transpose(Rz_crop) * projected';
    shadow_poly = polyshape(projected_ENZ(1, :), projected_ENZ(2, :));
end

function sample_heights = get_shading_sample_heights(params)
    if isfield(params, 'shading_sample_heights_m') && ~isempty(params.shading_sample_heights_m)
        sample_heights = params.shading_sample_heights_m;
    else
        sample_heights = 0;
    end

    sample_heights = sample_heights(:).';
    if any(~isfinite(sample_heights)) || any(sample_heights < 0)
        error('shading_sample_heights_m must contain finite nonnegative heights.');
    end
end

function bbox = poly_bounds(poly)
    [x, y] = boundary(poly);
    bbox = [min(x), max(x), min(y), max(y)];
end

function tf = boxes_overlap(a, b)
    tf = a(1) <= b(2) && a(2) >= b(1) && a(3) <= b(4) && a(4) >= b(3);
end

function wrapped = wrap_to_180(angle_deg)
    wrapped = mod(angle_deg + 180, 360) - 180;
end
