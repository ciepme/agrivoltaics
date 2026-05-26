function crop_GCF = ground_cover_factor(var, agriParams)
    
    if agriParams.tracking_mode == 1
        % TRACKING LAYOUT: Tractors drive N-S between the rows. 
        % The aisle width is measured E-W.
        unit_row_pitch = var.PV_w_p + var.PV_y_p; 
    else
        % FIXED TILT LAYOUT: Tractors drive E-W between the rows. 
        % The aisle width is measured N-S (accounting for panel tilt).
        unit_row_pitch = var.PV_l_p * cos(var.PV_sigma) + var.PV_y_p; 
    end
    
    % Subtract the unplantable areas
    unplantable_width = agriParams.support_width + agriParams.harvest_clearance + agriParams.safety_buffer;
    plantable_width = unit_row_pitch - unplantable_width;
    
    % Failsafe: if panels are too close together, a tractor can't fit at all!
    if plantable_width < 0
        plantable_width = 0; 
        warning('Panels are too close together! Tractors cannot fit. Yield will be 0.');
    end
    
    % Update the GCF dynamically
    crop_GCF = plantable_width / unit_row_pitch;
    
end