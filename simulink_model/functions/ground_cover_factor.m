function crop_GCF = ground_cover_factor(var, agriParams)
    if isfield(agriParams, 'geometry_mode') && agriParams.geometry_mode == 1
        crop_GCF = agriParams.hedge_width / agriParams.row_pitch;
        return;
    end
    
    % Total distance from one solar row to the next
    unit_row_width = var.PV_w_p + var.PV_y_p; 
    
    % Subtract the unplantable areas
    unplantable_width = agriParams.support_width + agriParams.harvest_clearance + agriParams.safety_buffer;
    plantable_width = unit_row_width - unplantable_width;
    
    % Failsafe: if panels are too close together, a tractor can't fit at all!
    if plantable_width < 0
        plantable_width = 0; 
        warning('Panels are too close together! Tractors cannot fit. Yield will be 0.');
    end
    
    % Update the GCF dynamically inside the struct
    crop_GCF = plantable_width / unit_row_width;
    
end
