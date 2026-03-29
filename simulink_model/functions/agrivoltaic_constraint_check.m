function isValid = agrivoltaic_constraint_check(custom_var)
    if custom_var.PV_z_p > 2
        isValid = false; % Set isValid to true if the condition is met
    else
        isValid = true; % Set isValid to false otherwise
    end
    