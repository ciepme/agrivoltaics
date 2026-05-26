function print_design_summary(var, params)
    if isfield(params, 'geometry_mode') && params.geometry_mode == 1
        print_row_centered_design_summary(var, params);
        return;
    end

    fprintf('\n=== OPTIMIZED DESIGN VARIABLES ===\n');
    fprintf('Geometry mode              : legacy\n');
    fprintf('Tracking mode              : %s\n', tracking_mode_label(params));
    fprintf('Panel height               : %.2f m\n', var.PV_z_p);
    fprintf('Panel length               : %.2f m\n', var.PV_l_p);
    fprintf('Panel width                : %.2f m\n', var.PV_w_p);
    fprintf('Azimuth                    : %.2f rad (%.1f deg)\n', var.PV_phi, rad2deg(var.PV_phi));
    fprintf('Tilt                       : %.2f rad (%.1f deg)\n', var.PV_sigma, rad2deg(var.PV_sigma));
    if isfield(params, 'tracking_mode') && params.tracking_mode == 1
        fprintf('Cross-row gap (PV_y_p)     : %.2f m\n', var.PV_y_p);
        fprintf('Along-row gap (PV_x_p)     : %.2f m\n', var.PV_x_p);
    else
        fprintf('Panel row gap (PV_y_p)     : %.2f m\n', var.PV_y_p);
        fprintf('Panel pair gap (PV_x_p)    : %.2f m\n', var.PV_x_p);
    end
end

function print_row_centered_design_summary(var, params)
    row_pitch = params.row_pitch;
    hedge_width = params.hedge_width;
    available_half_gap = (row_pitch - var.PV_w_p) / 2;
    lateral_offset_m = var.PV_d_norm * available_half_gap;
    offset_label = offset_direction_label(lateral_offset_m);

    fprintf('\n=== OPTIMIZED DESIGN VARIABLES ===\n');
    fprintf('Geometry mode              : row_centered\n');
    fprintf('Tracking mode              : %s\n', tracking_mode_label(params));
    fprintf('Panel height               : %.2f m\n', var.PV_z_p);
    fprintf('Panel length along row     : %.2f m (fixed by slice length)\n', var.PV_l_p);
    fprintf('Panel span across rows     : %.2f m (tilting dimension)\n', var.PV_w_p);
    fprintf('Raspberry row count        : %.0f rows\n', params.row_count);
    fprintf('Raspberry row pitch        : %.2f m (fixed = land_y / row_count)\n', row_pitch);
    fprintf('Raspberry hedge width      : %.2f m (%.2f ft)\n', hedge_width, hedge_width / 0.3048);
    fprintf('Crop ground cover fraction : %.3f (hedge_width / row_pitch)\n', hedge_width / row_pitch);
    fprintf('Panel center offset d_norm : %.2f (unitless)\n', var.PV_d_norm);
    if abs(lateral_offset_m) < 1e-9
        fprintf('Panel center offset        : 0.00 m (centered on hedge center)\n');
    else
        fprintf('Panel center offset        : %.2f m %s\n', abs(lateral_offset_m), offset_label);
    end
    fprintf('Offset scale               : d_norm * ((row_pitch - panel_span) / 2)\n');
    fprintf('Maximum center offset      : %.2f m when |d_norm| = 1\n', available_half_gap);

    if params.tracking_mode == 1
        fprintf('Tracker azimuth reference  : %.2f rad (%.1f deg)\n', var.PV_phi, rad2deg(var.PV_phi));
        fprintf('Hourly tracking angles     : optimized separately for 96 representative hours\n');
    else
        fprintf('Fixed panel azimuth        : %.2f rad (%.1f deg)\n', var.PV_phi, rad2deg(var.PV_phi));
        fprintf('Fixed panel tilt           : %.2f rad (%.1f deg)\n', var.PV_sigma, rad2deg(var.PV_sigma));
    end

    fprintf('Derived PV_y_p             : %.2f m (open cross-row gap, not row pitch)\n', var.PV_y_p);
    fprintf('Derived PV_x_p             : %.2f m\n', var.PV_x_p);
end

function label = tracking_mode_label(params)
    if isfield(params, 'tracking_mode') && params.tracking_mode == 1
        label = 'single-axis';
    else
        label = 'fixed-axis';
    end
end

function label = offset_direction_label(offset_m)
    tol = 1e-9;
    if offset_m > tol
        label = 'east of hedge center (local +x)';
    elseif offset_m < -tol
        label = 'west of hedge center (local -x)';
    else
        label = 'at';
    end
end
