function results = agrivoltaic_wrapper(custom_var, agriParams)

    %set custom variables
    agriVar = agriVarArray2Struct(custom_var, agriParams);
    agriVar = orderfields(agriVar);

    %create Simulation objects
    agriSim = Simulink.SimulationInput('agrivoltaics_v1');

    %set local as default variable (normal looks in base workspace)
    agriSim = agriSim.setVariable('agriVar', agriVar);
    agriSim = agriSim.setVariable('agriParams', agriParams);

    %fprintf('z_p: %.2f\n', agriVar.PV_z_p);

    %run sim
    out = sim(agriSim);

    %fprintf('P: %.2f\n', out.p.Data);

    results = [out.e.Data, out.p.Data];
end