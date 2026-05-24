function Cost = runEconomicModelNEW(params, height, width, modulelength, system_size)
% runEconomicModelNEW
% Small-plot agrivoltaic PV economic model.
%
% Inputs:
%   params       : agriParams struct
%   height       : module mounting height above ground, m
%   width        : module width, m
%   modulelength : module length, m
%   system_size  : rated DC system size, kWdc
%
% Output:
%   Cost         : present-value system cost over investigation period, $
%
% Notes:
%   - system_size is assumed to be kWdc.
%   - This model is intended for small agrivoltaic plots, not 100 MW
%     utility-scale systems.
%   - Cost = initial CapEx + discounted annual OpEx.
%   - No revenue, NPV, LCOE, or developer profit is included here.

    %% ------------------------------------------------------------
    % 0. Initialize and read key parameters
    % -------------------------------------------------------------
    Cost = 0.0;

    tracking = params.tracking_mode;
    efficiency = params.PV_n_p;
    system_size_kwdc = system_size;

    if system_size_kwdc <= 0 || efficiency <= 0
        Cost = 0.0;
        return;
    end

    discount_rate = get_param_default(params, 'discount_rate', 0.07);
    investigation_period = get_param_default(params, 'investigation_period', 20);
    startup_period = get_param_default(params, 'PV_startup_period', 1);

    projectduration = investigation_period;

    %% ------------------------------------------------------------
    % 1. Module cost parameters
    % Mostly from your friend's bottom-up model
    % -------------------------------------------------------------

    moduleyieldloss = 0.01;
    moduleweight = 12.56;       % kg/m2
    celltomodule = 0.98;
    moduledepreciation = 0.118;
    moduleprofitmsp = 0.25;

    % Intrinsic module cost elements
    cell = 87.72;               % $/m2
    frame = 1.62;               % $/m
    frontglass = 4.78;          % $/m2
    backsheet = 2.56;           % $/m2
    encapsulant = 2.43;         % $/m2
    junctionboxes = 3.74;       % $
    othermaterial = 4.41;       % $
    directlabor = 5.54;
    facilitiesoverhead = 0.88;
    businessoverhead = 2.95 + (28630000 / 1500000);
    depreciation = 56.87;
    shipping = 0.32;

    %% ------------------------------------------------------------
    % 2. Inverter cost parameters
    % -------------------------------------------------------------

    pcbassemblies = 186;
    enclosure = 48.45;
    directlaborinvert = 33.93;
    facilitiesoverheadinvert = 1.7;
    businessoverheadinvert = 2.7 + (58200000 / 2000000);
    depreciationinvert = 47.98;
    shippinginvert = 1.09;

    inverterweight = 0.935;     % kg/kWac
    inverterdepreciation = 0.118;
    inverteryieldloss = 0.01;
    inverterprofitmsp = 0.25;

    %% ------------------------------------------------------------
    % 3. Shared financial and BOS parameters
    % -------------------------------------------------------------

    salestaxrate = 0.058;
    contingencyrate = 0.0215;
    overheadrate = 0.01;

    % Loss / replacement / insurance rates
    ModuleLossRate   = 0.002;
    InverterLossRate = 0.09;
    SBOSLossRate     = 0.002;
    EBOSLossRate     = 0.002;
    InsuranceRate    = 0.0025;

    % EBOS / pad parameters
    EBOSPadCapacity = 50;       % kWac/m2 pad

    % Tariff
    Tariff301 = 0.25;

    %% ------------------------------------------------------------
    % 4. Tracking-specific assumptions
    % -------------------------------------------------------------

    if tracking == 0
        mode = 'fixed';

        % Fixed-axis systems often use slightly higher ILR
        ILR = 1.37;

        % Height scaling
        freightweightelasticity = 0.23;

        % O&M fixed adder for small plot
        % Keep this modest; crop/farm O&M is handled elsewhere.
        small_fixed_om_per_year = 500;   % $/yr

    else
        mode = 'single_axis';

        % Single-axis tracker ILR
        ILR = 1.34;

        freightweightelasticity = 0.20;

        % Slightly higher fixed O&M for tracking hardware
        small_fixed_om_per_year = 750;   % $/yr
    end

    %% ------------------------------------------------------------
    % 5. Module cost calculation
    %
    % This produces a unitized module cost approximately in $/kWdc.
    % -------------------------------------------------------------

    ModuleCostWithoutProfit = ...
        cell * ((1 / celltomodule) * (1 / (1 - moduleyieldloss))) + ...
        frame * (2 * (width + modulelength) / ...
        (width * modulelength * efficiency) * (1 / (1 - moduleyieldloss))) + ...
        (frontglass + backsheet + encapsulant) * ...
        ((1 / efficiency) * (1 / (1 - moduleyieldloss))) + ...
        (junctionboxes + othermaterial) * ...
        (1 / (efficiency * modulelength * width) * (1 / (1 - moduleyieldloss))) + ...
        (directlabor * facilitiesoverhead) * (1 / efficiency) + ...
        businessoverhead + ...
        depreciation * moduledepreciation + ...
        shipping * (moduleweight / efficiency);

    ModuleCost = ...
        moduleprofitmsp * ...
        (ModuleCostWithoutProfit - depreciation * moduledepreciation + depreciation) / 4 + ...
        ModuleCostWithoutProfit;

    %% ------------------------------------------------------------
    % 6. Inverter cost calculation
    %
    % Unitized inverter cost, approximately $/kWac.
    % -------------------------------------------------------------

    InverterCostWithoutProfit = ...
        (enclosure + pcbassemblies) * (1 / (1 - inverteryieldloss)) + ...
        directlaborinvert + ...
        facilitiesoverheadinvert + ...
        businessoverheadinvert;

    InverterCost = ...
        inverterprofitmsp * (InverterCostWithoutProfit + depreciation) / 4 + ...
        InverterCostWithoutProfit + ...
        depreciationinvert * inverterdepreciation + ...
        shippinginvert * inverterweight;

    %% ------------------------------------------------------------
    % 7. Agrivoltaic height scaling
    % -------------------------------------------------------------

    % Your design lower bound is 2.5 m. Treat extra height above that as
    % the agrivoltaic structural penalty.
    baseline_height = 2.5;
    heightdelta = max(0, height - baseline_height);

    rackingheightelasticity = 0.35;

    Piersweight = 5.68 * (1 + rackingheightelasticity * heightdelta^2);

    %% ------------------------------------------------------------
    % 8. SBOS, EBOS, fieldwork, officework
    %
    % These are unitized rates that ultimately feed into $/kWdc.
    % -------------------------------------------------------------

    if tracking == 0
        %% Fixed-axis / non-tracking structure

        Shippingcosting = 7 * (1 + freightweightelasticity * heightdelta);

        % SBOS without tracker-specific drive hardware
        SBOS = ...
            Piersweight + ...
            2.32 * 1.5 + ...
            1.07 + ...
            0.36 + ...
            346.4 * efficiency / (ILR * EBOSPadCapacity) + ...
            Shippingcosting;

        % EBOS
        EBOS = ...
            (30.01 + 5.46 + 6.71) * ILR + ...
            51.08 + ...
            4.69 + ...
            12.03 + ...
            56.59 + ...
            51.26 + ...
            2.26 + ...
            (15.14 + 1883000 / 74627) + ...
            131000 / 74627;

        laborcomplexityfactor = 1 + 0.25 * heightdelta;
        rentaldurationfactor  = 1 + 0.30 * heightdelta;

        Fieldwork = ...
            9.48 * laborcomplexityfactor + ...
            52.86 * 1.2 + ...
            (1 + 0.1 * heightdelta) * 10.54 + ...
            6.84 + ...
            9 * rentaldurationfactor + ...
            1.39;

        Officework = ...
            (ILR / efficiency) * (1.04 + 0.1 + 1.82) + ...
            (0.21 + 58.39 + 3.31) + ...
            (131200 + 478600 + 411400 + 32800) / 74627;

    else
        %% Single-axis tracking structure

        Shippingcosting = 10 * (1 + freightweightelasticity * heightdelta);

        % SBOS with tracker-specific hardware
        SBOS = ...
            7.44 + ...
            Piersweight + ...
            2.32 + ...
            1.07 + ...
            0.36 + ...
            3.63 + ...
            3.28 + ...
            1.83 + ...
            1.42 + ...
            346.4 * efficiency / (ILR * EBOSPadCapacity) + ...
            Shippingcosting;

        EBOS = ...
            (30.01 + 5.46 + 6.71) * ILR + ...
            51.08 + ...
            4.69 + ...
            12.03 + ...
            56.59 + ...
            51.26 + ...
            2.26 + ...
            (15.14 + 1883000 / 74627) + ...
            131000 / 74627;

        laborcomplexityfactor = 1 + 0.25 * heightdelta;
        rentaldurationfactor  = 1 + 0.30 * heightdelta;

        Fieldwork = ...
            9.48 * laborcomplexityfactor + ...
            52.86 * (efficiency / ILR) + ...
            (1 + 0.1 * heightdelta) * (6.92 + 1654000 / 456621) + ...
            6.84 + ...
            11.48 * rentaldurationfactor + ...
            1.39;

        Officework = ...
            (ILR / efficiency) * (1.04 + 0.1 + 1.82) + ...
            (0.21 + 58.39 + 3.31) + ...
            (131200 + 478600 + 411400 + 32800) / 74627;
    end

    %% ------------------------------------------------------------
    % 9. Other upfront project costs
    % -------------------------------------------------------------

    Salestax = ...
        (ModuleCost + InverterCost + SBOS + EBOS) * salestaxrate;

    Contingency = ...
        (ModuleCost + InverterCost + SBOS + EBOS + Fieldwork + Officework) * contingencyrate;

    Management = ...
        (ModuleCost + InverterCost + SBOS + EBOS + Fieldwork + Officework) * overheadrate + ...
        1363000 / 200000;

    Other = Salestax + Contingency + Management;

    %% ------------------------------------------------------------
    % 10. Unitized CapEx rate
    %
    % Approximate units: $/kWdc
    % -------------------------------------------------------------

    capex_rate_per_kwdc = ...
        ModuleCost + ...
        (InverterCost + EBOS + Officework) * (1 / ILR) + ...
        (SBOS + Fieldwork) * (1 / efficiency) + ...
        Other;

    % Optional user calibration multiplier
    % Use params.PV_cost_multiplier if you want to tune to quotes/literature.
    cost_multiplier = get_param_default(params, 'PV_cost_multiplier', 1.0);

    total_capex = system_size_kwdc * capex_rate_per_kwdc * cost_multiplier;

    %% ------------------------------------------------------------
    % 11. Annual OpEx
    %
    % This is separated from CapEx and discounted over the project period.
    % No large utility-scale fixed O&M charge is used.
    % -------------------------------------------------------------

    annual_om_rate_per_kwdc = ...
        (0.55 + 0.65) * (1 / efficiency) + ...
        ModuleCost * ModuleLossRate + ...
        InverterCost * InverterLossRate / ILR + ...
        SBOS * SBOSLossRate / efficiency + ...
        EBOS * EBOSLossRate / ILR + ...
        InsuranceRate * (SBOS + EBOS + ModuleCost + InverterCost);

    % Optional extra O&M per kWdc-year.
    % Default is zero because the bottom-up annual O&M already includes
    % replacement/inspection/insurance-like terms.
    extra_om_per_kwdc_yr = get_param_default(params, 'PV_extra_om_per_kwdc_yr', 0.0);

    annual_opex = ...
        small_fixed_om_per_year + ...
        system_size_kwdc * (annual_om_rate_per_kwdc + extra_om_per_kwdc_yr);

    %% ------------------------------------------------------------
    % 12. Discounted full-period cost
    % -------------------------------------------------------------

    exponent_value = startup_period:projectduration-1;

    discounted_opex = sum(annual_opex ./ ...
        (1 + discount_rate).^exponent_value);

    Cost = total_capex + discounted_opex;

end

%% ========================================================================
% Helper: parameter fallback
% ========================================================================
function value = get_param_default(params, fieldname, default_value)

    if isfield(params, fieldname)
        value = params.(fieldname);
    else
        value = default_value;
    end

end
% function Cost = runEconomicModelNEW(params, height, width, length, system_size)
% 
%     % 1. Determine Tracking Mode
%     if params.tracking_mode == 0
%         mode = 'fixed';
%     else
%         mode = 'single_axis';
%     end
% 
%     % 2. Load baseline NREL parameters
%     p = load_params(mode);
% 
%     % 3. Overwrite NREL defaults with YOUR simulation parameters
%     p.ModuleEfficiency = params.PV_n_p;
%     p.ModuleHeight     = length; % Length of panel in y_p direction
%     p.ModuleWidth      = width; % Width of panel in x_p direction
%     p.DiscountRate     = params.discount_rate;
%     p.LifetimeYears    = params.investigation_period;
%     p.ElectricityRate  = params.crop_elec_price;
%     p.OM_management_per_kwdc_yr = 20; %$/kW extra to help scale for larger size projects
% 
%     system_size_kwdc = system_size; 
%     module_area_m2 = system_size_kwdc / p.ModuleEfficiency;
% 
%     % If 0 panels (GA failed layout), return 0 costs to avoid NaN errors
%     if system_size_kwdc == 0
%         Cost = 0;
%         return;
%     end
% 
%     %% --- Agrivoltaic Cost Multipliers (Using agriVar.PV_z_p/height) ---
%     baseline_h = 1.5; % Standard utility height
%     steel_multiplier = max(1, height / baseline_h);
%     labor_multiplier = max(1, 1 + (height - baseline_h) * 0.15); 
% 
%     %% --- Derived sizing ---
%     kwac         = system_size_kwdc / p.ILR;
%     land_area_ha = module_area_m2   / p.LandArea;
% 
%     %% --- Component-level MSP ---
%     msp_module   = p.ModuleMarketPrice;  
%     msp_inverter = p.InverterMarketPrice;
% 
%     % SBOS — dynamically scaled by the steel_multiplier
%     if strcmp(mode, 'fixed')
%         sbos_torque_tubes  = p.TorqueTubeWeight * 2.5 * steel_multiplier;
%         sbos_driven_piers  = (p.PiersPerRow / (p.RowLength * p.ModuleHeight)) * 200 * steel_multiplier;
%         sbos_module_rails  = (1 / p.ModuleWidth) * 2.1 * (1 + p.Tariff301);
%         sbos_fasteners     = p.FastenerWeight    * 2.2 * (1 + p.Tariff301);
%         msp_sbos_per_m2    = sbos_torque_tubes + sbos_driven_piers + sbos_module_rails + sbos_fasteners;
%     else
%         sbos_torque_tubes  = (p.TorqueTubeWeight * 2.7 - 0.5 * p.TorqueTube45X * p.TorqueTubeWeight) * steel_multiplier;
%         sbos_driven_piers  = (p.PiersPerRow / (p.RowLength * p.ModuleHeight)) * 100 * steel_multiplier;
%         sbos_module_rails  = (1 / p.ModuleWidth) * 1.8 * (1 + p.Tariff301);
%         sbos_fasteners     = p.FastenerWeight    * 3.3 - 0.5 * p.Fastener45X * p.FastenerWeight;
%         sbos_slew_drive    = (1 / (p.ModuleHeight * p.RowLength)) * 300 * (1 + p.Tariff301);
%         sbos_dampers       = (p.PiersPerRow / (p.RowLength * p.ModuleHeight)) * 6.4 * (1 + p.Tariff301);
%         sbos_motor         = (1 / (p.ModuleHeight * p.RowLength)) * 440 * (1 + p.Tariff301);
%         sbos_controls      = (1 / (p.ModuleHeight * p.RowLength)) * 110 * (1 + p.Tariff301);
%         msp_sbos_per_m2    = sbos_torque_tubes + sbos_driven_piers + sbos_module_rails + ...
%                              sbos_fasteners    + sbos_slew_drive   + sbos_dampers      + ...
%                              sbos_motor        + sbos_controls;
%     end
% 
%     % EBOS 
%     if strcmp(mode, 'fixed')
%         msp_ebos_per_kwac = 90 + 30*(1+p.Tariff301) + 60 + 25 + 30*(1+p.Tariff301);
%     else
%         msp_ebos_per_kwac = 40 + 6.5*(1+p.Tariff301) + 15*(1+p.Tariff301) + 21 + ...
%                             9.6*(1+p.Tariff301) + 9.5 + 33*(1+p.Tariff301) + 17 + 66;
%     end
% 
%     % Fieldwork — scaled by labor_multiplier
%     ebos_labor_rate  = p.PVElectricalLabor  * (1 + p.LaborBurdenRate);
%     sbos_labor_rate  = p.PVConstructionLabor * (1 + p.LaborBurdenRate) * labor_multiplier;
% 
%     msp_fieldwork_per_m2 = ebos_labor_rate * p.ElectricalWage + ...
%                            sbos_labor_rate * p.ConstructionWage + ...
%                            p.SitePrep_per_m2 + p.EquipmentRental_per_m2;
% 
%     % Officework
%     msp_officework_per_kwdc = p.Officework_per_kwdc;
% 
%     %% --- System-level CapEx ---
%     capex_hardware       = (system_size_kwdc * msp_module) + (kwac * msp_inverter) + (module_area_m2 * msp_sbos_per_m2) + (kwac * msp_ebos_per_kwac);
%     capex_all_excl_other = capex_hardware + (module_area_m2 * msp_fieldwork_per_m2) + (system_size_kwdc * msp_officework_per_kwdc);
% 
%     capex_sales_tax   = capex_hardware       * p.SalesTaxRate;
%     capex_contingency = capex_all_excl_other * p.ContingencyRate;
%     capex_overhead    = capex_all_excl_other * p.OverheadRate;
%     capex_profit      = (capex_all_excl_other + capex_sales_tax + capex_contingency + capex_overhead) * p.DeveloperProfitMSP;
% 
%     total_capex = capex_all_excl_other + capex_sales_tax + capex_contingency + capex_overhead + capex_profit;
% 
%     %% --- Annual OpEx ---
%     om_cleaning      = module_area_m2   * (1/p.ModuleEfficiency) * p.Cleaning_per_m2;
%     om_inspection    = module_area_m2   * (1/p.ModuleEfficiency) * p.Inspection_per_m2;
%     om_new_bos       = module_area_m2   * msp_sbos_per_m2  * p.PartsLossRate;
%     om_new_modules   = system_size_kwdc * msp_module       * p.ModuleLossRate;
%     om_new_inverters = kwac             * msp_inverter     * p.InverterLossRate;
% 
%     annual_opex = om_cleaning + om_inspection + om_new_bos + om_new_modules + ...
%                   om_new_inverters + (land_area_ha * p.LandLease_per_ha) + ...
%                   (total_capex * p.PropertyTaxRate) + (total_capex * p.InsuranceRate) + p.OM_management_fixed + p.OM_management_per_kwdc_yr*system_size_kwdc;
% 
% exponent_value = params.PV_startup_period:p.LifetimeYears-1;
% 
% discounted_opex = sum(annual_opex ./ ...
%     (1 + p.DiscountRate).^exponent_value);
% 
% Cost = total_capex + discounted_opex;
% 
% end
% %%
% % 
% function p = load_params(mode)
% % load_params  Returns factor struct for the selected racking mode.
% % Values pulled directly from the PVSCM .txt parameter files.
% 
% % --- Shared factors (identical in both files) ---
% p.ModuleEfficiency    = 0.200;
% p.ModuleHeight        = 0;   % overridden for single-axis below
% p.ModuleWidth         = 0;   % overridden for single-axis below
% p.ModuleWeight        = 11.5;
% p.CellToModule        = 0.98;
% p.ModuleLabor         = 0.109;
% p.ModuleElectricity   = 5.79;
% p.ModuleDepreciation  = 0.118;
% p.ModuleProfitMSP     = 0.15;
% p.ModuleMarketPrice   = 336;
% 
% p.ESSWeight           = 9.4;
% p.ESSLabor            = 0.22;
% p.ESSElectricity      = 30.5;
% p.ESSDepreciation     = 0.11;
% p.ESSProfitMSP        = 0.25;
% p.ESSMarketPrice      = 228;
% p.Battery45X          = 10;
% p.BatteryDuration     = 4;
% p.ESS_ILR             = 1;
% 
% p.EBOSweight          = 2;
% p.LaborBurdenRate     = 0.54;
% p.ESSInstallLabor     = 2.15;
% 
% p.SalesTaxRate        = 0.058;
% p.OverheadRate        = 0.01;
% p.PartsLossRate       = 0.002;
% p.ModuleLossRate      = 0.001;
% p.InverterLossRate    = 0.09;
% p.ESSLossRate         = 0.025;
% p.LandArea            = 3000;   % m2 modules per ha
% p.PropertyTaxRate     = 0.002;
% p.InsuranceRate       = 0.0025;
% p.Tariff301           = 0.25;
% p.Tariff301high       = 0.50;
% p.Tariff201           = 0.1425;
% 
% % Revenue / financial (not in source files — set per your model)
% p.CapacityFactor      = 0.18;
% p.ElectricityRate     = 0.10;
% p.DiscountRate        = 0.07;
% p.LifetimeYears       = 20;
% 
% if strcmp(mode, 'fixed')
%     % --- Fixed-tilt parameters (FixedAxis .txt) ---
%     p.ILR                    = 1.37;
%     p.ModuleHeight           = 0;
%     p.ModuleWidth            = 0;
%     p.InverterWeight         = 0.789;
%     p.InverterLabor          = 0.045;
%     p.InverterElectricity    = 8.42;
%     p.InverterDepreciation   = 0.125;
%     p.InverterProfitMSP      = 0.15;
%     p.InverterMarketPrice    = 75;     % string inverter
%     p.FastenerWeight         = 0.2;
%     p.TorqueTubeWeight       = 5.9;
%     p.TorqueTube45X          = 0.87;   % 45X tax credit $/kg
%     p.Fastener45X            = 2.28;
%     p.SBOSshippingWeight     = 10;
%     p.PiersPerRow            = 11;
%     p.RowLength              = 86.3;
%     p.PVElectricalLabor      = 0.55;   % hr/m2 — lower efficiency for commercial
%     p.PVConstructionLabor    = 0.55;
%     p.ElectricalWage         = 35;     % $/hr
%     p.ConstructionWage       = 22.3;
%     p.SitePrep_per_m2        = 8;
%     p.EquipmentRental_per_m2 = 11;
%     p.Officework_per_kwdc    = 3122.04;
%     p.ContingencyRate        = 0.05;
%     p.DeveloperProfitMSP     = 0.08;
%     p.LandLease_per_ha       = 1500.6;
%     p.OM_management_fixed    = 8853.30;
%     p.Cleaning_per_m2        = 0.6;
%     p.Inspection_per_m2      = 0.5389;
%     p.StorageDuration        = 2;
% 
% else
%     % --- Single-axis tracking parameters (SingleAxis .txt) ---
%     p.ILR                    = 1.32;
%     p.ModuleHeight           = 0;  
%     p.ModuleWidth            = 0;
%     p.InverterWeight         = 0.935;
%     p.InverterLabor          = 0.045;
%     p.InverterElectricity    = 8.42;
%     p.InverterDepreciation   = 0.118;
%     p.InverterProfitMSP      = 0.25;
%     p.InverterMarketPrice    = 44;     % central inverter (utility scale)
%     p.FastenerWeight         = 0.2;
%     p.TorqueTubeWeight       = 4.8;
%     p.TorqueTube45X          = 0.87;   % 45X tax credit $/kg
%     p.Fastener45X            = 2.28;
%     p.SBOSshippingWeight     = 10;
%     p.PiersPerRow            = 11;
%     p.RowLength              = 86.3;
%     p.PVElectricalLabor      = 0.25;   % hr/m2 — higher efficiency at utility scale
%     p.PVConstructionLabor    = 0.40;
%     p.ElectricalWage         = 31;
%     p.ConstructionWage       = 22;
%     p.SitePrep_per_m2        = 8;
%     p.EquipmentRental_per_m2 = 10;
%     p.Officework_per_kwdc    = (32 + 66) / 2;  % midpoint of benchmark range
%     p.ContingencyRate        = 0.025;  % lower contingency at utility scale
%     p.DeveloperProfitMSP     = 0.05;   % lower margin at utility scale
%     p.LandLease_per_ha       = 1800;
%     p.OM_management_fixed    = 15000; % larger fixed management cost
%     p.Cleaning_per_m2        = 0.54;
%     p.Inspection_per_m2      = 0.64;
%     p.StorageDuration        = 2.4;
% end
% end