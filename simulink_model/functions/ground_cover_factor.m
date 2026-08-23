function crop_GCF = ground_cover_factor(var, agriParams)
    
    if agriParams.tracking_mode == 1

        unit_row_pitch = var.PV_l_p + var.PV_y_p; %total width of a row, support to support
    else

        unit_row_pitch = var.PV_l_p * cos(var.PV_sigma) + var.PV_y_p; 
    end
    
    % Subtract the unplantable areas
    if agriParams.persona == 1 %this is Fred's persona code
    unplantable_width = 2*agriParams.harvest_clearance + var.PV_l_p; %buffer on each side + solar panel length
    else
    unplantable_width = 2*agriParams.harvest_clearance; %buffer on each side, solar above the person
    end
    
    plantable_width = unit_row_pitch - unplantable_width;
    
    land_utilization_fraction = plantable_width / unit_row_pitch;
    % Update the GCF dynamically
    crop_GCF = land_utilization_fraction * agriParams.crop_canopy_cover;
    
end