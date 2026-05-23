function geom = build_row_centered_geometry(var, params, slice_length)
    if nargin < 3 || isempty(slice_length)
        slice_length = params.row_length / params.slice_count;
    end

    required_params = {'row_pitch', 'row_length', 'row_count', 'hedge_width', 'fixed_panel_length'};
    for i = 1:numel(required_params)
        if ~isfield(params, required_params{i})
            error('Missing row-centered geometry parameter: %s', required_params{i});
        end
    end

    required_vars = {'PV_z_p', 'PV_l_p', 'PV_w_p', 'PV_d_norm'};
    for i = 1:numel(required_vars)
        if ~isfield(var, required_vars{i})
            error('Missing row-centered geometry variable: %s', required_vars{i});
        end
    end

    row_pitch = params.row_pitch;
    hedge_width = params.hedge_width;
    panel_span = var.PV_w_p;
    panel_length = var.PV_l_p;
    panel_height = var.PV_z_p;

    if row_pitch <= 0 || hedge_width <= 0 || slice_length <= 0
        error('Row pitch, hedge width, and slice length must be positive.');
    end
    if panel_span <= 0 || panel_length <= 0 || panel_height < 0
        error('Panel span and length must be positive, and panel height must be nonnegative.');
    end
    if panel_span > row_pitch
        error('PV_w_p cannot exceed row_pitch in row-centered geometry.');
    end
    if ~isfinite(var.PV_d_norm) || abs(var.PV_d_norm) > 1 + 1e-9
        error('PV_d_norm must be finite and within [-1, 1].');
    end

    hedge_center = row_pitch / 2;
    lateral_offset = var.PV_d_norm * (row_pitch - panel_span) / 2;
    panel_center_x = hedge_center + lateral_offset;

    row_offsets = get_offsets(params, 'shading_optimization_row_offsets', -2:2);
    slice_offsets = get_offsets(params, 'shading_optimization_slice_offsets', -10:10);

    if ~any(row_offsets == 0) || ~any(slice_offsets == 0)
        error('Row-centered shading source offsets must include zero in both dimensions.');
    end

    x_shifts = row_offsets(:) * row_pitch;
    y_shifts = slice_offsets(:) * slice_length;
    [x_grid, y_grid] = ndgrid(x_shifts, y_shifts);

    panel_centers = [
        panel_center_x + x_grid(:), ...
        (slice_length / 2) + y_grid(:), ...
        repmat(panel_height, numel(x_grid), 1)
    ];
    panel_offsets = [
        x_grid(:) ./ row_pitch, ...
        y_grid(:) ./ slice_length
    ];

    geom = struct();
    geom.row_pitch = row_pitch;
    geom.row_length = params.row_length;
    geom.row_count = params.row_count;
    geom.hedge_width = hedge_width;
    geom.hedge_center = hedge_center;
    geom.slice_length = slice_length;
    geom.panel_span = panel_span;
    geom.panel_length = panel_length;
    geom.panel_height = panel_height;
    geom.panel_center_x = panel_center_x;
    geom.lateral_offset = lateral_offset;
    geom.source_row_offsets = row_offsets;
    geom.source_slice_offsets = slice_offsets;
    geom.panel_centers = panel_centers;
    geom.panel_offsets = panel_offsets;
    geom.unit_cell = polyshape([0, row_pitch, row_pitch, 0], [0, 0, slice_length, slice_length]);
    geom.hedge_shape = polyshape( ...
        [hedge_center - hedge_width / 2, hedge_center + hedge_width / 2, hedge_center + hedge_width / 2, hedge_center - hedge_width / 2], ...
        [0, 0, slice_length, slice_length]);
    geom.unit_area = row_pitch * slice_length;
    geom.hedge_area = hedge_width * slice_length;
    geom.crop_azimuth_deg = rad2deg(params.land_angle);
end

function offsets = get_offsets(params, field_name, default_offsets)
    if isfield(params, field_name)
        offsets = params.(field_name);
    else
        offsets = default_offsets;
    end

    if isempty(offsets) || ~isnumeric(offsets) || any(~isfinite(offsets(:)))
        error('%s must be a finite numeric vector.', field_name);
    end

    offsets = unique(round(offsets(:).'), 'stable');
end
