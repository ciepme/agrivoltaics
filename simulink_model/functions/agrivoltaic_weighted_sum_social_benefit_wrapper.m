function cost = agrivoltaic_weighted_sum_social_benefit_wrapper(custom_var, agriParams, delta)
    results = agrivoltaic_wrapper(custom_var, agriParams);
    result_struct = agriObjArray2Struct(results);
    cost = -1.*(delta.*(result_struct.social_profit) + (1-delta).*result_struct.yearly_biomass);
    closeExcel;
end