function results = agrivoltaic_wrapper(custom_var, agriParams)

    %set custom variables
    agriVar = agriVarArray2Struct(custom_var, agriParams);
    agriVar = orderfields(agriVar);

    %create Simulation objects
    agriSim = Simulink.SimulationInput('agrivoltaics_v1');

    %set local as default variable (normal looks in base workspace)
    agriSim = agriSim.setVariable('agriVar', agriVar);
    agriSim = agriSim.setVariable('agriParams', agriParams);

    %set in rapid accelerator
    %agriSim = agriSim.setModelParameter('SimulationMode','Rapid');
    %agriSim = agriSim.setModelParameter('RapidAcceleratorUpToDateCheck','off'); 

    %fprintf('z_p: %.2f\n', agriVar.PV_z_p);

    %run sim
    out = sim(agriSim);

    %fprintf('P: %.2f\n', out.p.Data);

    E = out.e.Data;
    P = out.p.Data;

    social_cost = -1.*(P + 190 .* (E ./ 1000));

    crop_revenue = out.crop_revenue.Data;
    total_panels = out.total_panels.Data;
    pv_revenue = out.pv_revenue.Data;
    yearly_biomass = out.yearly_biomass.Data;

    results = [E, P, social_cost, pv_revenue, crop_revenue, yearly_biomass, total_panels];
    closeExcel;
end