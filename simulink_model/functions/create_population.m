function x0 = create_population(pop_size)
    agrivoltaics_variable_definition;

    % Tell GA exactly how many variables we are optimizing (7 or 103, dependent
    %create a smart guess to seed GA population for if tracking mode
    x0_base = agriVarStruct2Array(agriVar, agriParams);

    % on fixed or single-axis)
    num_vars = length(lb);
    %adds an InitialPopulationMatrix for better initial guess
    
    x0 = agriVarStruct2Array(agriVar, agriParams);
    
    % build a better initial population
    % Make sure consistent
    num_vars = length(lb);
    pop = zeros(pop_size, num_vars);
    % 1First member = physics-based smart guess
    pop(1,:) = x0;
    
    % population = small perturbations
    for i = 2:pop_size
        candidate = x0;
    
        if agriParams.tracking_mode == 1
            % Only perturb tracking angles (more stable)
            idx = tracking_angle_indices(agriParams);
    
            noise = 0.1 * agriParams.PV_max_tilt * randn(size(idx));
            candidate(idx) = candidate(idx) + noise;
    
            % Clamp
            candidate(idx) = max(candidate(idx), lb(idx));
            candidate(idx) = min(candidate(idx), ub(idx));
        else
            % Fixed-axis: perturb all vars slightly
            noise = 0.05 * (ub - lb) .* randn(1, num_vars);
            candidate = candidate + noise;
    
            candidate = max(candidate, lb);
            candidate = min(candidate, ub);
        end
    
        pop(i,:) = candidate;
    end
end
