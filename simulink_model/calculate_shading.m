function SF_hourly = calculate_shading(var, params)
    % 1. Initialize output array (24 hours)
    SF_hourly = zeros(24, 1);
    
    % 2. Extract parameters
    PV_length = var.PV.l_p;
    PV_width = var.PV.w_p;
    PV_tilt = rad2deg(var.PV.sigma);    
    PV_azimuth = rad2deg(var.PV.phi);   
    PV_height = var.PV.z_p;
    
    row_distance = var.PV.y_p;
    pair_distance = var.PV.x_p;
    crop_azimuth = rad2deg(params.land.angle);

    % 3. THE SPEED HACK: Define the "Unit Cell"
    % Instead of the whole farm, we look at the land for just ONE panel
    unit_width = PV_width + pair_distance;
    unit_length = PV_length + row_distance;
    
    % Place our single panel exactly in the center of this unit cell
    n_pv = 1;
    PV_pairs = [unit_width/2, unit_length/2, PV_height];
    
    % 4. Loop through the 24-hour sun data
    for t = 1:24
        sun_alt_deg = rad2deg(params.weather.beta_s(t));
        sun_az_deg = rad2deg(params.weather.phi_s(t));
        
        if sun_alt_deg > 0
            % Calculate shading for just this single unit cell!
            SF_temp = PV_shading_factor_calc(sun_az_deg, sun_alt_deg, ...
                                        PV_azimuth, PV_tilt, ...
                                        n_pv, PV_pairs, ...
                                        PV_length, PV_width, ...
                                        unit_width, unit_length, crop_azimuth);
            SF_hourly(t) = SF_temp;
        else
            SF_hourly(t) = 0; 
        end
    end
end

%% --- Helper Functions ---
function SF = PV_shading_factor_calc(sun_azimuth, sun_altitude, PV_azimuth, PV_tilt, n_pv, PV_pairs, PV_length, PV_width, crop_width, crop_length, crop_azimuth)
    shading_area = repmat(polyshape(), 1, n_pv);
    for pv = 1:n_pv
        shading_area(pv) = shade_projection(PV_pairs(pv,:)', PV_width, PV_length, PV_azimuth, PV_tilt, sun_azimuth, sun_altitude, crop_width, crop_length, crop_azimuth);
    end      
    total_shadow = union(shading_area);
    SF = area(total_shadow) / (crop_width * crop_length);
end

function shading_box = shade_projection(center_point, PV_width, PV_length, PV_azimuth, PV_tilt, sun_azimuth, sun_altitude, crop_width, crop_length, crop_azimuth)
    panel_2D = [PV_width/2, PV_length/2, 0; -PV_width/2, PV_length/2, 0; -PV_width/2, -PV_length/2, 0; PV_width/2, -PV_length/2,  0];
    Rx_pv = [1,0,0; 0, cosd(PV_tilt), -sind(PV_tilt); 0, sind(PV_tilt), cosd(PV_tilt)];
    Rz_pv = [-sind(PV_azimuth), -cosd(PV_azimuth), 0; cosd(PV_azimuth), -sind(PV_azimuth), 0; 0, 0, 1];
    Rz_crop = [-sind(crop_azimuth), -cosd(crop_azimuth), 0; cosd(crop_azimuth), -sind(crop_azimuth), 0; 0, 0, 1]; 
    
    panel_global_pv = repmat(Rz_crop*center_point, 1, 4) + Rz_pv*Rx_pv*panel_2D';
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