function struck = agriVarArray2Struct(AgriVar_array, agriParams)

    if isfield(agriParams, 'geometry_mode') && agriParams.geometry_mode == 1
        n = numel(AgriVar_array);
        if agriParams.tracking_mode == 1
            expected_n = 99;
        else
            expected_n = 4;
        end

        if n ~= expected_n
            if n == 7 || n == 103
                error('Stale design vector length %d. Row-centered vectors use 4 variables in fixed mode or 99 variables in single-axis mode.', n);
            end
            error('Unexpected design vector length %d. Expected %d for the current row-centered tracking mode.', n, expected_n);
        end

        struck.PV_z_p = AgriVar_array(1);
        struck.PV_l_p = agriParams.fixed_panel_length;
        struck.PV_w_p = AgriVar_array(2);

        if agriParams.tracking_mode == 1
            struck.PV_phi = agriParams.land_angle;
            struck.PV_sigma = 0;
            struck.PV_d_norm = AgriVar_array(3);
            struck.tracking_angles = reshape(AgriVar_array(4:end), [24, 4])';
        else
            struck.PV_phi = agriParams.row_centered_fixed_phi;
            struck.PV_sigma = AgriVar_array(3);
            struck.PV_d_norm = AgriVar_array(4);
            struck.tracking_angles = zeros(4, 24);
        end

        struck.PV_x_p = 0;
        struck.PV_y_p = agriParams.row_pitch - struck.PV_w_p;

        struck = orderfields(struck);
        return;
    end

    % Base 7 variables
    struck.PV_z_p = AgriVar_array(1);
    struck.PV_l_p = AgriVar_array(2);
    struck.PV_w_p = AgriVar_array(3);
    struck.PV_phi = AgriVar_array(4);
    struck.PV_sigma = AgriVar_array(5);
    struck.PV_y_p = AgriVar_array(6);
    struck.PV_x_p = AgriVar_array(7);

    % --- Tracking handling ---
    if agriParams.tracking_mode == 1
        tracking_flat = AgriVar_array(8:end);

        % reshape into [4 x 24]
        struck.tracking_angles = reshape(tracking_flat, [24,4])';
    else
        % MUST exist for bus compatibility
        struck.tracking_angles = zeros(4,24);
    end

    % Enforce field order (VERY IMPORTANT for Simulink)
    struck = orderfields(struck);
end
