function results = agrivoltaic_wrapper(custom_var, agriParams)

    %set custom variables
    agriVar = agriVarArray2Struct(custom_var);

    %create Simulation objects
    agriSim = Simulink.SimulationInput('agrivoltaics_v1');
    %set local as default variable (normal looks in base workspace)
    agriSim = agriSim.setVariable('agriVar', agriVar);
    %ensuring params get updated in parallel workspaces
    agriSim = agriSim.setVariable('agriParams', agriParams);

    %run sim
    out = sim(agriSim);

    % changing to a struct for results
    % results = [out.e.Data, out.p.Data];
    results.E                          = out.e.Data;
    results.P                          = out.p.Data;
    results.crop_revenue               = out.crop_revenue.Data;
    results.pv_revenue                 = out.pv_revenue.Data;
    results.pv_cost                    = out.pv_cost.Data;
    results.crop_cost                  = out.crop_cost.Data;
    results.crop_income                = out.crop_income.Data;
    results.pv_income                  = out.pv_income.Data;
    results.raspberry_yield_kgs_annual = out.raspberry_yield_kgs_annual.Data;
end