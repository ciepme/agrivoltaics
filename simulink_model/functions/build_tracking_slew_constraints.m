function [A, B, Aeq, Beq] = build_tracking_slew_constraints(num_vars, agriParams)
    A = [];
    B = [];
    Aeq = [];
    Beq = [];

    if agriParams.tracking_mode ~= 1
        return;
    end

    tracking_idx = tracking_angle_indices(agriParams);
    if numel(tracking_idx) ~= 96
        error('Expected 96 tracking variables, found %d.', numel(tracking_idx));
    end

    max_slew_per_hour = agriParams.max_slew_per_hour;
    num_steps = 23;
    total_constraints = 4 * num_steps * 2;

    A = zeros(total_constraints, num_vars);
    B = ones(total_constraints, 1) * max_slew_per_hour;

    row = 1;
    for s = 1:4
        season_idx = tracking_idx((s-1)*24 + (1:24));
        for h = 1:num_steps
            v1 = season_idx(h);
            v2 = season_idx(h + 1);

            A(row, v1) = -1;
            A(row, v2) =  1;
            row = row + 1;

            A(row, v1) =  1;
            A(row, v2) = -1;
            row = row + 1;
        end
    end
end
