function print_value_breakdown(x_opt, agriParams)
    results = agrivoltaic_wrapper(x_opt, agriParams);

    E            = results(1);
    P            = results(2);
    r_el         = results(4);
    r_cr         = results(5);
    yearly_biomass = results(6);
    total_panels = results(7);

    land_m2    = agriParams.land_x * agriParams.land_y;
    land_ha    = land_m2 / 10000;
    land_acres = land_m2 / 4046.86;
    inv_yrs    = agriParams.investigation_period;
    disc       = agriParams.discount_rate * 100;

    E_tons       = E / 1000;
    scc          = 190 * E_tons;
    social_value = P + scc;


    fprintf('\n=== AGRIVOLTAICS VALUE BREAKDOWN ===\n');
    fprintf('NPV over %d years | Discount rate: %.1f%%\n', inv_yrs, disc);
    fprintf('Land: %.0f m2 (%.3f acres)\n\n', land_m2, land_acres);

    fprintf('MODEL OUTPUTS\n');
    fprintf('  Total profit (P)         : $%12.2f\n', P);
    fprintf('  PV electricity revenue   : $%12.2f\n', r_el);
    fprintf('  Raspberry crop revenue   : $%12.2f\n', r_cr);
    fprintf('  Emissions value          : $%12.2f  ($190/ton x %.1f tons)\n', scc, E_tons);
    fprintf('  Total social value       : $%12.2f\n', social_value);
    fprintf('  Total CO2 displaced      : %12.1f tons\n', E_tons);
    fprintf('  Yearly biomass           : %12.1f kg/yr (%.0f kg/ha/yr)\n', ...
        yearly_biomass, yearly_biomass / land_ha);
    fprintf('  Total panels             : %12.0f\n\n', total_panels);

    fprintf('NET\n');
    fprintf('  Total profit (P)         : $%12.2f\n', P);
    fprintf('  Total social value       : $%12.2f\n', social_value);
    fprintf('  Emissions value          : $%12.2f\n\n', scc);

    % Plot hourly tracking schedule if single-axis mode is active and tracking vars exist.
    if isfield(agriParams, 'tracking_mode') && agriParams.tracking_mode == 1
        try
            agriVar_plot = agriVarArray2Struct(x_opt, agriParams);
            if ~isfield(agriVar_plot, 'tracking_angles')
                error('tracking_angles not found in agriVar.');
            end

            tracking_deg = rad2deg(agriVar_plot.tracking_angles); % [4 x 24]
            max_tilt_deg = rad2deg(agriParams.PV_max_tilt);

            fprintf('TRACKING SUMMARY\n');
            fprintf('  Tracking mode          : single-axis (hourly)\n');
            fprintf('  Tilt range used        : %.2f to %.2f deg (mean %.2f deg)\n', ...
                min(tracking_deg(:)), max(tracking_deg(:)), mean(tracking_deg(:)));
            fprintf('  0 deg meaning          : panel horizontal/flat\n');
            fprintf('  Max tilt meaning       : tracker mechanical limit at %.1f deg\n\n', max_tilt_deg);

            if isfield(agriParams, 'weather') && ...
               isfield(agriParams.weather, 'spring') && ...
               isfield(agriParams.weather, 'summer') && ...
               isfield(agriParams.weather, 'fall') && ...
               isfield(agriParams.weather, 'winter')
                irr_mat = [
                    local_total_irradiance(agriParams.weather.spring).';
                    local_total_irradiance(agriParams.weather.summer).';
                    local_total_irradiance(agriParams.weather.fall).';
                    local_total_irradiance(agriParams.weather.winter).'
                ]; % [4 x 24]
            else
                irr_mat = nan(4,24);
            end

            % Flatten by season blocks: spring(1:24), summer(25:48), ...
            x_plot = 1:96;
            tracking_deg_vec = reshape(tracking_deg.', [1,96]);
            irr_vec = reshape(irr_mat.', [1,96]);

            figure('Name','Single-Axis Tracking Schedule','NumberTitle','off', ...
                   'Position',[100 120 1500 520]);

            yyaxis left;
            plot(x_plot, tracking_deg_vec, '-o', 'LineWidth', 1.3, 'MarkerSize', 3, ...
                 'DisplayName', 'Tracking tilt');
            ylabel('Tilt Magnitude (deg)');
            grid on;
            hold on;

            yline(0, 'k--', 'LineWidth', 1.1, ...
                'DisplayName', '0 deg (panel horizontal)');

            yline(max_tilt_deg, 'r--', 'LineWidth', 1.1, ...
                'DisplayName', sprintf('Max tilt %.1f deg (mechanical limit)', max_tilt_deg));

            for b = [24.5, 48.5, 72.5]
                xline(b, 'Color', [0.25, 0.25, 0.25], 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');
            end

            season_names = {'Spring','Summer','Fall','Winter'};
            season_mid = [12.5, 36.5, 60.5, 84.5];
            y_lim_left = ylim;
            y_text = y_lim_left(2) - 0.04 * (y_lim_left(2) - y_lim_left(1));
            for s = 1:4
                text(season_mid(s), y_text, season_names{s}, ...
                    'HorizontalAlignment', 'center', ...
                    'FontWeight', 'bold', ...
                    'Color', [0.1, 0.1, 0.1]);
            end

            yyaxis right;
            if any(~isnan(irr_vec))
                plot(x_plot, irr_vec, '-', 'Color', [0.90, 0.45, 0.10], ...
                    'LineWidth', 1.3, ...
                    'DisplayName', 'Hourly horizontal irradiance (DNI*sin(beta) + DHI)');
            end
            ylabel('Irradiance (W/m^2)');

            % Full 96-hour labeling: 12 AM ... 11 PM repeated for each season day
            hour_labels_24 = local_hour_labels();
            xticks(1:96);
            xticklabels(repmat(hour_labels_24, 1, 4));
            xtickangle(90);
            xlim([1, 96]);
            xlabel('Hour of day for each representative season day - Spring/Summer/Fall/Winter');

            title('Optimized Hourly Single-Axis Tracking Schedule by Season');
            legend('Location', 'northoutside', 'Orientation', 'horizontal');
            hold off;
        catch ME
            fprintf('Tracking plot skipped: %s\n\n', ME.message);
        end
    end
end

function labels = local_hour_labels()
    labels = cell(1,24);
    for h = 0:23
        if h == 0
            labels{h+1} = '12 AM';
        elseif h < 12
            labels{h+1} = sprintf('%d AM', h);
        elseif h == 12
            labels{h+1} = '12 PM';
        else
            labels{h+1} = sprintf('%d PM', h - 12);
        end
    end
end

function irr = local_total_irradiance(season_weather)
    dni = season_weather.DNI(:);
    dhi = season_weather.DHI(:);
    beta_s = season_weather.beta_s(:);

    if ~(numel(dni) == 24 && numel(dhi) == 24 && numel(beta_s) == 24)
        error('Expected 24 hourly values for DNI/DHI/beta_s.');
    end

    % Approximate global horizontal irradiance (GHI)
    irr = dni .* max(0, sin(beta_s)) + dhi;
end
