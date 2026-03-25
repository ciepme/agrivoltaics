function print_value_breakdown(x_opt, agriParams)
    results = agrivoltaic_wrapper(x_opt, agriParams);

    E        = results.E;
    P        = results.P;
    r_cr     = results.crop_revenue;
    r_el     = results.pv_revenue;
    c_pv     = results.pv_cost;
    c_cr     = results.crop_cost;
    inc_cr   = results.crop_income;
    inc_pv   = results.pv_income;
    yield_kg = results.raspberry_yield_kgs_annual;

    land_m2    = agriParams.land_x * agriParams.land_y;
    land_ha    = land_m2 / 10000;
    land_acres = land_m2 / 4046.86;
    inv_yrs    = agriParams.investigation_period;
    disc       = agriParams.discount_rate * 100;

    E_tons       = E / 1000;
    scc          = 190 * E_tons;
    social_value = P + scc;


    fprintf('\n=== AGRIVOLTAICS VALUE BREAKDOWN ===\n');
    fprintf('NPV over %d years | Discount rate: %.1f%%\n', inv_yrs, disc);
    fprintf('Land: %.0f m2 (%.3f acres)\n\n', land_m2, land_acres);

    fprintf('REVENUES (20-yr NPV)\n');
    fprintf('  Raspberry crop revenue   : $%12.2f\n', r_cr);
    fprintf('  PV electricity revenue   : $%12.2f\n', r_el);
    fprintf('  Social cost of carbon    : $%12.2f  ($190/ton x %.1f tons)\n', scc, E_tons);
    fprintf('  -------------------------------------------\n');
    fprintf('  Total revenue            : $%12.2f\n\n', r_cr + r_el + scc);

    fprintf('COSTS (20-yr NPV)\n');
    fprintf('  PV costs (CAPEX + OPEX)  : $%12.2f\n', c_pv);
    fprintf('  Crop costs (all)         : $%12.2f\n', c_cr);
    fprintf('  ---------------------------------------------\n');
    fprintf('  Total costs              : $%12.2f\n\n', c_pv + c_cr);

    fprintf('NET\n');
    fprintf('  Total profit (P)         : $%12.2f\n', P);
    fprintf('  Total social value       : $%12.2f\n\n', social_value);

    fprintf('ANNUAL (undiscounted)\n');
    fprintf('  Annual crop income       : $%12.2f /yr\n', inc_cr);
    fprintf('  Annual PV income         : $%12.2f /yr\n', inc_pv);
    fprintf('  Annual CO2 displaced     : %12.1f tons/yr\n', E_tons / inv_yrs);
    fprintf('  Raspberry yield          : %12.1f kg/yr (%.0f kg/ha/yr)\n', ...
        yield_kg, yield_kg / land_ha);
end