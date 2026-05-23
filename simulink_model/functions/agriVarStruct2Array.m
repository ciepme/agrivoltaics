function array = agriVarStruct2Array(AgriVar, agriParams)

    if isfield(agriParams, 'geometry_mode') && agriParams.geometry_mode == 1
        if ~isfield(AgriVar, 'PV_d_norm')
            error('agriVar is missing PV_d_norm. Row-centered geometry requires hedge-relative offset.');
        end

        if agriParams.tracking_mode == 1
            if ~isfield(AgriVar, 'tracking_angles') || ~isequal(size(AgriVar.tracking_angles), [4, 24])
                error('Single-axis mode requires tracking_angles with size 4 x 24.');
            end

            tracking_flat = reshape(AgriVar.tracking_angles', [1, 96]);
            array = [
                AgriVar.PV_z_p, ...
                AgriVar.PV_w_p, ...
                AgriVar.PV_d_norm, ...
                tracking_flat
            ];
        else
            array = [
                AgriVar.PV_z_p, ...
                AgriVar.PV_w_p, ...
                AgriVar.PV_sigma, ...
                AgriVar.PV_d_norm
            ];
        end
        return;
    end

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
