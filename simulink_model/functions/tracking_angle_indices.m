function idx = tracking_angle_indices(agriParams)
    if agriParams.tracking_mode ~= 1
        idx = [];
        return;
    end

    if isfield(agriParams, 'geometry_mode') && agriParams.geometry_mode == 1
        first_idx = 4;
    else
        first_idx = 8;
    end

    idx = first_idx:(first_idx + 95);
end
