function Cost = runEconomicModelNEW(params, height, width, length, system_size)

    % 1. Determine Tracking Mode
    if params.tracking_mode == 0
        mode = 'fixed';
    else
        mode = 'single_axis';
    end

    % 2. Load baseline NREL parameters
    p = load_params(mode);

    % 3. Overwrite NREL defaults with YOUR simulation parameters
    p.ModuleEfficiency = params.PV_n_p;
    p.ModuleHeight     = length; % Length of panel in y_p direction
    p.ModuleWidth      = width; % Width of panel in x_p direction
    p.DiscountRate     = params.discount_rate;
    p.LifetimeYears    = params.investigation_period;
    p.ElectricityRate  = params.crop_elec_price;
    
    system_size_kwdc = system_size; 
    module_area_m2 = system_size_kwdc / p.ModuleEfficiency;
    
    % If 0 panels (GA failed layout), return 0 costs to avoid NaN errors
    if system_size_kwdc == 0
        Cost = 0;
        return;
    end

    %% --- Agrivoltaic Cost Multipliers (Using agriVar.PV_z_p/height) ---
    baseline_h = 1.5; % Standard utility height
    steel_multiplier = max(1, height / baseline_h);
    labor_multiplier = max(1, 1 + (height - baseline_h) * 0.15); 

    %% --- Derived sizing ---
    kwac         = system_size_kwdc / p.ILR;
    land_area_ha = module_area_m2   / p.LandArea;

    %% --- Component-level MSP ---
    msp_module   = p.ModuleMarketPrice;  
    msp_inverter = p.InverterMarketPrice;

    % SBOS — dynamically scaled by the steel_multiplier
    if strcmp(mode, 'fixed')
        sbos_torque_tubes  = p.TorqueTubeWeight * 2.5 * steel_multiplier;
        sbos_driven_piers  = (p.PiersPerRow / (p.RowLength * p.ModuleHeight)) * 200 * steel_multiplier;
        sbos_module_rails  = (1 / p.ModuleWidth) * 2.1 * (1 + p.Tariff301);
        sbos_fasteners     = p.FastenerWeight    * 2.2 * (1 + p.Tariff301);
        msp_sbos_per_m2    = sbos_torque_tubes + sbos_driven_piers + sbos_module_rails + sbos_fasteners;
    else
        sbos_torque_tubes  = (p.TorqueTubeWeight * 2.7 - 0.5 * p.TorqueTube45X * p.TorqueTubeWeight) * steel_multiplier;
        sbos_driven_piers  = (p.PiersPerRow / (p.RowLength * p.ModuleHeight)) * 100 * steel_multiplier;
        sbos_module_rails  = (1 / p.ModuleWidth) * 1.8 * (1 + p.Tariff301);
        sbos_fasteners     = p.FastenerWeight    * 3.3 - 0.5 * p.Fastener45X * p.FastenerWeight;
        sbos_slew_drive    = (1 / (p.ModuleHeight * p.RowLength)) * 300 * (1 + p.Tariff301);
        sbos_dampers       = (p.PiersPerRow / (p.RowLength * p.ModuleHeight)) * 6.4 * (1 + p.Tariff301);
        sbos_motor         = (1 / (p.ModuleHeight * p.RowLength)) * 440 * (1 + p.Tariff301);
        sbos_controls      = (1 / (p.ModuleHeight * p.RowLength)) * 110 * (1 + p.Tariff301);
        msp_sbos_per_m2    = sbos_torque_tubes + sbos_driven_piers + sbos_module_rails + ...
                             sbos_fasteners    + sbos_slew_drive   + sbos_dampers      + ...
                             sbos_motor        + sbos_controls;
    end

    % EBOS 
    if strcmp(mode, 'fixed')
        msp_ebos_per_kwac = 90 + 30*(1+p.Tariff301) + 60 + 25 + 30*(1+p.Tariff301);
    else
        msp_ebos_per_kwac = 40 + 6.5*(1+p.Tariff301) + 15*(1+p.Tariff301) + 21 + ...
                            9.6*(1+p.Tariff301) + 9.5 + 33*(1+p.Tariff301) + 17 + 66;
    end

    % Fieldwork — scaled by labor_multiplier
    ebos_labor_rate  = p.PVElectricalLabor  * (1 + p.LaborBurdenRate);
    sbos_labor_rate  = p.PVConstructionLabor * (1 + p.LaborBurdenRate) * labor_multiplier;

    msp_fieldwork_per_m2 = ebos_labor_rate * p.ElectricalWage + ...
                           sbos_labor_rate * p.ConstructionWage + ...
                           p.SitePrep_per_m2 + p.EquipmentRental_per_m2;

    % Officework
    msp_officework_per_kwdc = p.Officework_per_kwdc;

    %% --- System-level CapEx ---
    capex_hardware       = (system_size_kwdc * msp_module) + (kwac * msp_inverter) + (module_area_m2 * msp_sbos_per_m2) + (kwac * msp_ebos_per_kwac);
    capex_all_excl_other = capex_hardware + (module_area_m2 * msp_fieldwork_per_m2) + (system_size_kwdc * msp_officework_per_kwdc);
    
    capex_sales_tax   = capex_hardware       * p.SalesTaxRate;
    capex_contingency = capex_all_excl_other * p.ContingencyRate;
    capex_overhead    = capex_all_excl_other * p.OverheadRate;
    capex_profit      = (capex_all_excl_other + capex_sales_tax + capex_contingency + capex_overhead) * p.DeveloperProfitMSP;

    total_capex = capex_all_excl_other + capex_sales_tax + capex_contingency + capex_overhead + capex_profit;

    %% --- Annual OpEx ---
    om_cleaning      = module_area_m2   * (1/p.ModuleEfficiency) * p.Cleaning_per_m2;
    om_inspection    = module_area_m2   * (1/p.ModuleEfficiency) * p.Inspection_per_m2;
    om_new_bos       = module_area_m2   * msp_sbos_per_m2  * p.PartsLossRate;
    om_new_modules   = system_size_kwdc * msp_module       * p.ModuleLossRate;
    om_new_inverters = kwac             * msp_inverter     * p.InverterLossRate;
    
    annual_opex = om_cleaning + om_inspection + om_new_bos + om_new_modules + ...
                  om_new_inverters + (land_area_ha * p.LandLease_per_ha) + ...
                  (total_capex * p.PropertyTaxRate) + (total_capex * p.InsuranceRate) + p.OM_management_fixed;

exponent_value = params.PV_startup_period:p.LifetimeYears-1;

discounted_opex = sum(annual_opex ./ ...
    (1 + p.DiscountRate).^exponent_value);

Cost = total_capex + discounted_opex;

end
%%
% 
function p = load_params(mode)
% load_params  Returns factor struct for the selected racking mode.
% Values pulled directly from the PVSCM .txt parameter files.

% --- Shared factors (identical in both files) ---
p.ModuleEfficiency    = 0.200;
p.ModuleHeight        = 0;   % overridden for single-axis below
p.ModuleWidth         = 0;   % overridden for single-axis below
p.ModuleWeight        = 11.5;
p.CellToModule        = 0.98;
p.ModuleLabor         = 0.109;
p.ModuleElectricity   = 5.79;
p.ModuleDepreciation  = 0.118;
p.ModuleProfitMSP     = 0.15;
p.ModuleMarketPrice   = 336;

p.ESSWeight           = 9.4;
p.ESSLabor            = 0.22;
p.ESSElectricity      = 30.5;
p.ESSDepreciation     = 0.11;
p.ESSProfitMSP        = 0.25;
p.ESSMarketPrice      = 228;
p.Battery45X          = 10;
p.BatteryDuration     = 4;
p.ESS_ILR             = 1;

p.EBOSweight          = 2;
p.LaborBurdenRate     = 0.54;
p.ESSInstallLabor     = 2.15;

p.SalesTaxRate        = 0.058;
p.OverheadRate        = 0.01;
p.PartsLossRate       = 0.002;
p.ModuleLossRate      = 0.001;
p.InverterLossRate    = 0.09;
p.ESSLossRate         = 0.025;
p.LandArea            = 3000;   % m2 modules per ha
p.PropertyTaxRate     = 0.002;
p.InsuranceRate       = 0.0025;
p.Tariff301           = 0.25;
p.Tariff301high       = 0.50;
p.Tariff201           = 0.1425;

% Revenue / financial (not in source files — set per your model)
p.CapacityFactor      = 0.18;
p.ElectricityRate     = 0.10;
p.DiscountRate        = 0.07;
p.LifetimeYears       = 20;

if strcmp(mode, 'fixed')
    % --- Fixed-tilt parameters (FixedAxis .txt) ---
    p.ILR                    = 1.37;
    p.ModuleHeight           = 0;
    p.ModuleWidth            = 0;
    p.InverterWeight         = 0.789;
    p.InverterLabor          = 0.045;
    p.InverterElectricity    = 8.42;
    p.InverterDepreciation   = 0.125;
    p.InverterProfitMSP      = 0.15;
    p.InverterMarketPrice    = 75;     % string inverter
    p.FastenerWeight         = 0.2;
    p.TorqueTubeWeight       = 5.9;
    p.TorqueTube45X          = 0.87;   % 45X tax credit $/kg
    p.Fastener45X            = 2.28;
    p.SBOSshippingWeight     = 10;
    p.PiersPerRow            = 11;
    p.RowLength              = 86.3;
    p.PVElectricalLabor      = 0.55;   % hr/m2 — lower efficiency for commercial
    p.PVConstructionLabor    = 0.55;
    p.ElectricalWage         = 35;     % $/hr
    p.ConstructionWage       = 22.3;
    p.SitePrep_per_m2        = 8;
    p.EquipmentRental_per_m2 = 11;
    p.Officework_per_kwdc    = 3122.04;
    p.ContingencyRate        = 0.05;
    p.DeveloperProfitMSP     = 0.08;
    p.LandLease_per_ha       = 1500.6;
    p.OM_management_fixed    = 8853.30;
    p.Cleaning_per_m2        = 0.6;
    p.Inspection_per_m2      = 0.5389;
    p.StorageDuration        = 2;

else
    % --- Single-axis tracking parameters (SingleAxis .txt) ---
    p.ILR                    = 1.32;
    p.ModuleHeight           = 0;  
    p.ModuleWidth            = 0;
    p.InverterWeight         = 0.935;
    p.InverterLabor          = 0.045;
    p.InverterElectricity    = 8.42;
    p.InverterDepreciation   = 0.118;
    p.InverterProfitMSP      = 0.25;
    p.InverterMarketPrice    = 44;     % central inverter (utility scale)
    p.FastenerWeight         = 0.2;
    p.TorqueTubeWeight       = 4.8;
    p.TorqueTube45X          = 0.87;   % 45X tax credit $/kg
    p.Fastener45X            = 2.28;
    p.SBOSshippingWeight     = 10;
    p.PiersPerRow            = 11;
    p.RowLength              = 86.3;
    p.PVElectricalLabor      = 0.25;   % hr/m2 — higher efficiency at utility scale
    p.PVConstructionLabor    = 0.40;
    p.ElectricalWage         = 31;
    p.ConstructionWage       = 22;
    p.SitePrep_per_m2        = 8;
    p.EquipmentRental_per_m2 = 10;
    p.Officework_per_kwdc    = (32 + 66) / 2;  % midpoint of benchmark range
    p.ContingencyRate        = 0.025;  % lower contingency at utility scale
    p.DeveloperProfitMSP     = 0.05;   % lower margin at utility scale
    p.LandLease_per_ha       = 1800;
    p.OM_management_fixed    = 180000; % larger fixed management cost
    p.Cleaning_per_m2        = 0.54;
    p.Inspection_per_m2      = 0.64;
    p.StorageDuration        = 2.4;
end
end