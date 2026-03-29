function array = agriVarStruct2Array(AgriVar, agriParams)

    % Base variables
    base_array = [
        AgriVar.PV_z_p, ...
        AgriVar.PV_l_p, ...
        AgriVar.PV_w_p, ...
        AgriVar.PV_phi, ...
        AgriVar.PV_sigma, ...
        AgriVar.PV_y_p, ...
        AgriVar.PV_x_p
    ];

    if agriParams.tracking_mode == 1
        % Flatten [4 x 24] → [1 x 96]
        tracking_flat = reshape(AgriVar.tracking_angles', [1,96]);

        array = [base_array, tracking_flat];
    else
        array = base_array;
    end
end