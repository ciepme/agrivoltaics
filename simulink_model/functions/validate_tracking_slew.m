function [max_jump, max_violation] = validate_tracking_slew(x, agriParams)
    max_jump = 0;
    max_violation = 0;

    if agriParams.tracking_mode ~= 1
        return;
    end

    tracking_idx = tracking_angle_indices(agriParams);
    for s = 1:4
        season_idx = tracking_idx((s-1)*24 + (1:24));
        curve = x(season_idx);
        max_jump = max(max_jump, max(abs(diff(curve))));
    end

    [A, B] = build_tracking_slew_constraints(numel(x), agriParams);
    if ~isempty(A)
        max_violation = max(A*x(:) - B);
    end
end
