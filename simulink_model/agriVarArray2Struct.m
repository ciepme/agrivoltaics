function struck = agriVarArray2Struct(AgriVar_array)

    struck.PV_z_p = AgriVar_array(1);
    struck.PV_l_p = AgriVar_array(2);
    struck.PV_w_p = AgriVar_array(3);
    struck.PV_phi = AgriVar_array(4);
    struck.PV_sigma = AgriVar_array(5);
    struck.PV_y_p = AgriVar_array(6);
    struck.PV_x_p = AgriVar_array(7);
end