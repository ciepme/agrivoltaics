function pop = build_initial_population(agriVar, agriParams, lb, ub, pop_size, seed)
    if nargin >= 6 && ~isempty(seed)
        rng(seed);
    end

    x0 = agriVarStruct2Array(agriVar, agriParams);
    num_vars = numel(lb);
    pop = zeros(pop_size, num_vars);
    pop(1, :) = max(min(x0(:).', ub(:).'), lb(:).');

    tracking_idx = tracking_angle_indices(agriParams);
    if agriParams.tracking_mode == 1
        pop(1, :) = enforce_tracking_slew(pop(1, :), lb, ub, agriParams);
    end

    for i = 2:pop_size
        candidate = pop(1, :);

        if agriParams.tracking_mode == 1
            span = ub(tracking_idx) - lb(tracking_idx);
            random_angles = lb(tracking_idx) + rand(1, numel(tracking_idx)) .* span;
            smoothed_angles = smoothdata(random_angles, 'gaussian', 5);

            candidate(tracking_idx) = max(smoothed_angles, lb(tracking_idx));
            candidate(tracking_idx) = min(candidate(tracking_idx), ub(tracking_idx));
            candidate = enforce_tracking_slew(candidate, lb, ub, agriParams);
        else
            span = ub(1:num_vars) - lb(1:num_vars);
            candidate(1:num_vars) = lb(1:num_vars) + rand(1, num_vars) .* span;
            candidate = max(candidate, lb);
            candidate = min(candidate, ub);
        end

        pop(i, :) = candidate;
    end
end

function candidate = enforce_tracking_slew(candidate, lb, ub, agriParams)
    tracking_idx = tracking_angle_indices(agriParams);
    max_slew = agriParams.max_slew_per_hour;

    for s = 1:4
        season_idx = tracking_idx((s-1)*24 + (1:24));
        curve = candidate(season_idx);
        curve = enforce_curve_slew(curve, lb(season_idx), ub(season_idx), max_slew);
        candidate(season_idx) = curve;
    end
end

function curve = enforce_curve_slew(curve, lb_curve, ub_curve, max_slew)
    curve = max(curve, lb_curve);
    curve = min(curve, ub_curve);

    for iter = 1:10
        for h = 2:numel(curve)
            curve(h) = min(max(curve(h), curve(h-1) - max_slew), curve(h-1) + max_slew);
            curve(h) = max(curve(h), lb_curve(h));
            curve(h) = min(curve(h), ub_curve(h));
        end

        for h = numel(curve)-1:-1:1
            curve(h) = min(max(curve(h), curve(h+1) - max_slew), curve(h+1) + max_slew);
            curve(h) = max(curve(h), lb_curve(h));
            curve(h) = min(curve(h), ub_curve(h));
        end
    end
end
