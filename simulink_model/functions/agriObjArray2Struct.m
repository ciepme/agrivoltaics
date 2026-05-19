function struck = agriObjArray2Struct(obj_array)

    %results = [E, P, social_cost, pv_revenue, crop_revenue, yearly_biomass, total_panels];
    % 7 Outputs for each simulation
    struck.emission_reduction = obj_array(1);
    struck.fiscal_profit = obj_array(2);
    struck.social_profit = -obj_array(3);
    struck.PV_revenue = obj_array(4);
    struck.crop_revenue = obj_array(5);
    struck.yearly_biomass = obj_array(6);
    struck.n_panels = obj_array(7);
    struck.yearly_energy = obj_array(8);

    % Enforce field order (VERY IMPORTANT for Simulink)
    struck = orderfields(struck);
end