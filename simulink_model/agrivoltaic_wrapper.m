function results = agrivoltaic_wrapper(custom_var)

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
    var.PV.z_p =    custom_var_z_p;         % panel height (m) 
    var.PV.l_p =    custom_var_l_p;         % panel length (m)
    var.PV.w_p =    custom_var_w_p;         % panel width (m)
    var.PV.phi =    custom_var_phi;         % azimuth (rad)
    var.PV.sigma =  custom_var_sigma;       % tilt (rad)
    var.PV.psi =    custom_var_psi;         % field layout angle (rad)
    var.PV.y_p =    custom_var_y_p;         % row distance (m)
    var.PV.x_p =    custom_var_x_p;         % panel distance (m)

    out = sim("agrivoltaics_v1.slx");
    results = [out.e.Data, out.p.Data];
end