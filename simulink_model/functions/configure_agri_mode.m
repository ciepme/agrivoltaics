function [agriParams, agriVar, lb, ub] = configure_agri_mode(agriParams, tracking_mode, geometry_mode)
    if ischar(geometry_mode) || isstring(geometry_mode)
        geometry_is_row_centered = strcmpi(string(geometry_mode), "row_centered");
    else
        geometry_is_row_centered = geometry_mode == 1;
    end

    agriParams.tracking_mode = tracking_mode;
    agriParams.geometry_mode = double(geometry_is_row_centered);
    agriParams.shading_sample_heights_m = get_param_default(agriParams, 'shading_sample_heights_m', 0);

    lb = [2.5, 1.0, 1.0, -pi/2, 0, 2, 0.1];
    ub = [4.5, 2.5, 1.5,  pi/2, pi/2, 10.0, 1.0];

    agriVar = struct();
    agriVar.PV_z_p = 3.79;
    agriVar.PV_l_p = 2.43;
    agriVar.PV_w_p = 1.45;
    agriVar.PV_phi = 0;
    agriVar.PV_sigma = 0;
    agriVar.PV_y_p = 2.0;
    agriVar.PV_x_p = 0.1;
    agriVar.tracking_angles = zeros(4, 24);

    if ~geometry_is_row_centered && tracking_mode == 1
        agriVar.PV_l_p = 1.45;
        agriVar.PV_w_p = 2.43;
        ub(2:3) = [1.5, 2.5];
    end

    if geometry_is_row_centered
        agriParams.land_x = 50;
        agriParams.land_y = 50;
        agriParams.row_count = 11;
        agriParams.row_pitch = agriParams.land_y / agriParams.row_count;
        agriParams.row_length = agriParams.land_x;
        agriParams.slice_count = 50;
        agriParams.fixed_panel_length = agriParams.row_length / agriParams.slice_count;
        agriParams.hedge_width = get_param_default(agriParams, 'hedge_width', 2 * 0.3048);
        agriParams.PV_d_norm_min = get_param_default(agriParams, 'PV_d_norm_min', -1);
        agriParams.PV_d_norm_max = get_param_default(agriParams, 'PV_d_norm_max', 1);
        agriParams.row_centered_fixed_phi = get_param_default(agriParams, 'row_centered_fixed_phi', pi / 2);
        agriParams.shading_optimization_row_offsets = get_param_default(agriParams, 'shading_optimization_row_offsets', -2:2);
        agriParams.shading_optimization_slice_offsets = get_param_default(agriParams, 'shading_optimization_slice_offsets', -10:10);

        agriVar.PV_l_p = agriParams.fixed_panel_length;
        agriVar.PV_w_p = 1.75;
        agriVar.PV_x_p = 0;
        agriVar.PV_y_p = agriParams.row_pitch - agriVar.PV_w_p;
        agriVar.PV_d_norm = 0;

        if tracking_mode == 1
            agriVar.PV_phi = agriParams.land_angle;
            agriVar.PV_sigma = 0;
            lb = [2.5, 1.0, agriParams.PV_d_norm_min];
            ub = [4.5, 2.5, agriParams.PV_d_norm_max];
        else
            agriVar.PV_phi = agriParams.row_centered_fixed_phi;
            agriVar.PV_sigma = 0;
            lb = [2.5, 1.0, 0, agriParams.PV_d_norm_min];
            ub = [4.5, 2.5, pi/2, agriParams.PV_d_norm_max];
        end
    end

    if tracking_mode == 1
        agriVar.tracking_angles = generate_physics_tracking(agriParams, agriVar);
        agriVar.tracking_angles = max(agriVar.tracking_angles, -agriParams.PV_max_tilt);
        agriVar.tracking_angles = min(agriVar.tracking_angles,  agriParams.PV_max_tilt);

        [tracking_lb, tracking_ub, agriParams] = build_tracking_bounds(agriParams);
        lb = [lb, tracking_lb];
        ub = [ub, tracking_ub];
    end

    agriParams = orderfields(agriParams);
    agriVar = orderfields(agriVar);
end

function [tracking_lb, tracking_ub, agriParams] = build_tracking_bounds(agriParams)
    tracking_lb = zeros(1, 96);
    tracking_ub = zeros(1, 96);

    seasons = {'spring', 'summer', 'fall', 'winter'};
    max_tilt = agriParams.PV_max_tilt;
    max_slew = agriParams.max_slew_per_hour;
    ramp_hours = ceil(max_tilt / max_slew);

    agriParams.tracking_daytime_mask = false(4, 24);
    agriParams.tracking_move_mask = false(4, 24);

    for s = 1:4
        beta_s = agriParams.weather.(seasons{s}).beta_s;
        is_daytime = beta_s > 0;
        day_idx = find(is_daytime);
        can_move = false(1, 24);

        if ~isempty(day_idx)
            first_day_hour = day_idx(1);
            last_day_hour = day_idx(end);
            move_start = max(1, first_day_hour - ramp_hours);
            move_end = min(24, last_day_hour + ramp_hours);
            can_move(move_start:move_end) = true;
        end

        local_start = (s-1)*24 + 1;
        local_end = s*24;
        tracking_lb(local_start:local_end) = -max_tilt .* can_move;
        tracking_ub(local_start:local_end) =  max_tilt .* can_move;
        agriParams.tracking_daytime_mask(s, :) = is_daytime;
        agriParams.tracking_move_mask(s, :) = can_move;
    end
end

function value = get_param_default(params, fieldname, default_value)
    if isfield(params, fieldname)
        value = params.(fieldname);
    else
        value = default_value;
    end
end
