function results = agrivoltaic_wrapper(custom_var, agriParams)

    %set custom variables
agriVar = agriVarArray2Struct(custom_var, agriParams);
agriVar = orderfields(agriVar);

    %create Simulation objects
    agriSim = Simulink.SimulationInput('agrivoltaics_v1');

    %set local as default variable (normal looks in base workspace)
agriSim = agriSim.setVariable('var', agriVar);
agriSim = agriSim.setVariable('params', agriParams); 

    %run sim
    out = sim(agriSim);

    results = [out.e.Data, out.p.Data];
end