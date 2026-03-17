function array = agriVarStruct2Array(AgriVar)
    array = ones(1,7);
    
    array(1) = AgriVar.PV_z_p;
    array(2) = AgriVar.PV_l_p;
    array(3) = AgriVar.PV_w_p;
    array(4) = AgriVar.PV_phi;
    array(5) = AgriVar.PV_sigma;
    array(6) = AgriVar.PV_y_p;
    array(7) = AgriVar.PV_x_p;
end