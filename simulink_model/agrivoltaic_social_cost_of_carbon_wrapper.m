function social_cost = agrivoltaic_social_cost_of_carbon_wrapper(custom_var)
    results = agrivoltaic_wrapper(custom_var);
    E = results(1);
    P = results(2);
    social_cost = -1.*(P + 190 .* (E ./ 1000));
    %fprintf("Social Cost of %f\n", social_cost);
end