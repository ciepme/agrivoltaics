function [profit, pv_revenue_npv, crop_revenue_npv, total_capex] = value_model( ...
    P_annual, params, crop_income_annual, ...
    pv_capex, pv_discounted_opex, crop_capex, crop_discounted_opex)

    % Calculates the Net Present Value (NPV) and overall profit of the system

    %converting to kWh from Wh
    pv_annual_kWh = P_annual* 0.001;

    elec_sell_price      = params.elec_sell_price;
    discount_rate        = params.discount_rate;
    investigation_period = params.investigation_period;
    pv_startup_year      = params.PV_startup_period;
    crop_startup_year    = params.startup_years;
    elec_buy_price       = params.elec_buy_price;
    farm_elec_demand_kWh = params.farm_elec_demand_kWh;
    net_metering_mode    = params.net_metering_mode; 

    %Annual Electricity Value
    if net_metering_mode == 1
        % self-use, using the electricity rather than just selling back
        if pv_annual_kWh >= farm_elec_demand_kWh
            % solar covers all demand, plus excess is sold back to grid
            savings = farm_elec_demand_kWh * elec_buy_price;
            excess_sold = (pv_annual_kWh - farm_elec_demand_kWh) * elec_sell_price;
            pv_annual_value = savings + excess_sold;
        else
            % solar only partially covers demand
            pv_annual_value = pv_annual_kWh * elec_buy_price;
        end
    else
        % sell all generation straight to the grid
        pv_annual_value = pv_annual_kWh * elec_sell_price;
    end

    % 3. Calculate Discounted Cash Flows (NPV of Revenues)
    pv_years = pv_startup_year : investigation_period ;
    crop_years = crop_startup_year : (investigation_period - 1);

    pv_revenue_npv = sum(pv_annual_value ./ ((1 + discount_rate) .^ pv_years));
    crop_revenue_npv = sum(crop_income_annual ./ ((1 + discount_rate) .^ crop_years));

    % 4. Aggregate Costs
    total_capex = pv_capex + crop_capex;
    total_opex_npv = pv_discounted_opex + crop_discounted_opex;

    % 5. Calculate Final Project Profit
    profit = (pv_revenue_npv + crop_revenue_npv) - (total_capex + total_opex_npv);

end