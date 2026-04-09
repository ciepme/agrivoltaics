function cost = agrivoltaic_weighted_sum_emissions_wrapper(custom_var, agriParams, delta)
    results = agrivoltaic_wrapper(custom_var, agriParams);
    E = results(1);
    P = results(2);
    cost = -1.*(delta.*P + (1-delta).*E);
    closeExcel;
end