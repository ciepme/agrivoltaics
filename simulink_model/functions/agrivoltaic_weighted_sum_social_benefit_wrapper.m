function cost = agrivoltaic_weighted_sum_social_benefit_wrapper(custom_var, agriParams, delta)
    results = agrivoltaic_wrapper(custom_var, agriParams);
    social_benefit = -results(3);
    yearly_biomass = results(6);% changed from 7
    cost = -1.*(delta.*social_benefit + (1-delta).*yearly_biomass);
    closeExcel;
end