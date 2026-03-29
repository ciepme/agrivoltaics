function struck = agriVarArray2Struct(AgriVar_array, agriParams)

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