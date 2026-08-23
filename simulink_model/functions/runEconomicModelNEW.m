function [Cost, total_capex, discounted_opex] = runEconomicModelNEW(system_size_kW, height, tracking_mode, foundation_flag, investigation_period)
    % if no solar
    if system_size_kW <= 0
        Cost = 0; total_capex = 0; discounted_opex = 0;
        return;
    end
    
    % Scale and Baseline Hardware Costs
    module_cost_per_kW = 320; 
    size_min = 500;  
    size_max = 5000; 
    
    scale_factor = (system_size_kW - size_min) / (size_max - size_min);
    scale_factor = max(0, min(1, scale_factor));
    
    inverter_cost_per_kW = 120 - scale_factor * (120 - 60); 
    ebos_cost_per_kW     = 250 - scale_factor * (250 - 150);
    annual_om_per_kW     = 25  - scale_factor * (25 - 15);
    soft_cost_multiplier = 1.25 - scale_factor * (1.25 - 1.15);
    
    % Structural Installation Base Costs
    if tracking_mode == 1
        base_sbos_mat_per_kW = 180; 
        base_sbos_lab_per_kW = 120;
    else
        base_sbos_mat_per_kW = 110; 
        base_sbos_lab_per_kW = 90;
    end
    
    % Agrivoltaic Height Multipliers
    standard_height = 1.5;
    height_delta = max(0, height - standard_height);
    steel_multiplier = 1 + (0.15 * height_delta); 
    labor_multiplier = 1 + (0.20 * height_delta);
    
    % Apply Foundation Multipliers
    if foundation_flag == 1
        foundation_mat_mult = 1.0; 
        foundation_lab_mult = 1.0; 
    else
        foundation_mat_mult = 1.25; 
        foundation_lab_mult = 0.80; 
    end
    
    % Final CapEx Rates
    final_sbos_mat = base_sbos_mat_per_kW * steel_multiplier * foundation_mat_mult;
    final_sbos_lab = base_sbos_lab_per_kW * labor_multiplier * foundation_lab_mult;
    
    capex_rate_per_kW = module_cost_per_kW + inverter_cost_per_kW + ...
                        ebos_cost_per_kW + final_sbos_mat + final_sbos_lab;
                        
    % Calculate Base CapEx
    gross_capex = system_size_kW * capex_rate_per_kW * soft_cost_multiplier;
    
    % --- APPLY 30% INVESTMENT TAX CREDIT (ITC) ---
    % The federal government effectively covers 30% of the capital cost.
    total_capex = gross_capex * (1 - 0.30); 
    
    % OpEx (Net Present Value)
    discount_rate = 0.07;
    annual_opex = system_size_kW * annual_om_per_kW;
    
    years = 1:investigation_period;
    discounted_opex = sum(annual_opex ./ ((1 + discount_rate) .^ years));
    
    % Final Cost
    Cost = total_capex + discounted_opex;
end