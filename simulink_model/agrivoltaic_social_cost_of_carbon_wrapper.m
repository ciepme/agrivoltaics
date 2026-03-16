function social_profit = agrivoltaic_social_cost_of_carbon_wrapper(custom_var)
    %parse custom_var
    custom_var_z_p = custom_var(1);
    custom_var_l_p = custom_var(2);
    custom_var_w_p = custom_var(3);
    custom_var_phi = custom_var(4);
    custom_var_sigma = custom_var(5);
    custom_var_psi = custom_var(6);
    custom_var_y_p = custom_var(7);
    custom_var_x_p = custom_var(8);

    %set standard variables
    load("agrivoltaics_variable_definition_data.mat");

    %set custom variables
    agriVar.PV_z_p =    custom_var_z_p;         % panel height (m)
    agriVar.PV_l_p =    custom_var_l_p;         % panel length (m)
    agriVar.PV_w_p =    custom_var_w_p;         % panel width (m)
    agriVar.PV_phi =    custom_var_phi;         % azimuth (rad)
    agriVar.PV_sigma =  custom_var_sigma;       % tilt (rad)
    agriVar.PV_psi =    custom_var_psi;         % field layout angle (rad)
    agriVar.PV_y_p =    custom_var_y_p;         % row distance (m)
    agriVar.PV_x_p =    custom_var_x_p;         % panel distance (m)

    %create Simulation objects
    agriSim = Simulink.SimulationInput('agrivoltaics_v1');

    %set local as default variable (normal looks in base workspace)
    agriSim = agriSim.setVariable('agriVar', agriVar);
    agriSim = agriSim.setVariable('agriVar', agriVar);

    %run sim
    out = sim(agriSim);

    E = out.e.Data;
    P = out.p.Data;

    social_profit = P + 190 .* (E./1000);
end