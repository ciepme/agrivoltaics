%% inspect_GA_results
% Edit result_file to inspect a saved agrivoltaics_GA_main run.

clear;
clc;
bdclose('all');
addpath(genpath(pwd));

% Recreate Simulink bus objects required by agrivoltaic_wrapper. The saved
% result file provides the actual optimized agriParams used below.
agrivoltaics_variable_definition;

% Leave empty to inspect the most recently modified GA result file.
result_file = '';
snapshot_season = 'summer';
snapshot_hour = 10;
snapshot_height_m = 0;
% make_shading_gif = default_make_shading_gif(); % Change to true to write the 96-hour GIF.
make_shading_gif = true;
gif_delay_seconds = 0.25;
use_current_row_centered_shading_neighborhood = true;

if isempty(result_file)
    result_files = dir('agrivoltaics_GA_main_data_w_selector_*.mat');
    if isempty(result_files)
        error('No agrivoltaics_GA_main_data_w_selector_*.mat files found in %s.', pwd);
    end
    [~, newest_idx] = max([result_files.datenum]);
    result_file = result_files(newest_idx).name;
end

loaded = load(result_file);
required_vars = {'ga_solve', 'agriParams'};
for i = 1:numel(required_vars)
    if ~isfield(loaded, required_vars{i})
        error('Result file %s is missing %s.', result_file, required_vars{i});
    end
end

ga_solve = loaded.ga_solve;
agriParams = loaded.agriParams;

if use_current_row_centered_shading_neighborhood && ...
        isfield(agriParams, 'geometry_mode') && agriParams.geometry_mode == 1
    agriParams.shading_optimization_row_offsets = -2:2;
    agriParams.shading_optimization_slice_offsets = -10:10;
    fprintf('Using current row-centered shading neighborhood for inspection: 5 x 21 source cells.\n');
    fprintf('Set use_current_row_centered_shading_neighborhood = false to reproduce saved result-file parameters exactly.\n');
end

ga_var = agriVarArray2Struct(ga_solve, agriParams);

fprintf('Inspecting GA result file: %s\n', result_file);
if isfield(loaded, 'exitflag')
    fprintf('Exit flag: %d\n', loaded.exitflag);
end
if isfield(loaded, 'fval')
    fprintf('Objective value: %.2f\n', -loaded.fval);
end

print_design_summary(ga_var, agriParams);
print_value_breakdown(ga_solve, agriParams);

[SF_spring, SF_summer, SF_fall, SF_winter] = calculate_shading(ga_var, agriParams);
plot_shading_profiles(SF_spring, SF_summer, SF_fall, SF_winter);

if isfield(agriParams, 'tracking_mode') && agriParams.tracking_mode == 1
    plot_tracking_angles(ga_var, agriParams);
end

if isfield(agriParams, 'geometry_mode') && agriParams.geometry_mode == 1
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    outdir = fullfile(pwd, 'results', ['ga_inspection_' timestamp]);
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    save_row_centered_snapshot(ga_var, agriParams, snapshot_season, snapshot_hour, snapshot_height_m, outdir);
    fprintf('Saved row-centered shadow snapshot to %s\n', outdir);

    if make_shading_gif
        gif_file = save_row_centered_shading_gif(ga_var, agriParams, snapshot_height_m, gif_delay_seconds, outdir);
        fprintf('Saved 96-hour row-centered shading GIF to %s\n', gif_file);
    end
end

function plot_shading_profiles(SF_spring, SF_summer, SF_fall, SF_winter)
    hours = 0:23;
    shading_mat = [SF_spring(:), SF_summer(:), SF_fall(:), SF_winter(:)];

    figure('Name', 'Seasonal Shading Fractions', 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [100 100 980 500]);
    plot(hours, shading_mat, '-o', 'LineWidth', 1.4, 'MarkerSize', 4);
    grid on;
    xlabel('Hour of representative day');
    ylabel('Shading fraction');
    ylim([0 1]);
    legend({'Spring', 'Summer', 'Fall', 'Winter'}, 'Location', 'best');
    title('Hourly shading fraction by representative season');
end

function plot_tracking_angles(var, params)
    tracking_deg = rad2deg(var.tracking_angles);
    max_tilt_deg = rad2deg(params.PV_max_tilt);

    figure('Name', 'Optimized Tracking Angles', 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [120 120 980 500]);
    imagesc(1:24, 1:4, tracking_deg);
    colorbar;
    clim([-max_tilt_deg, max_tilt_deg]);
    xlabel('Hour of representative day');
    ylabel('Season');
    yticks(1:4);
    yticklabels({'Spring', 'Summer', 'Fall', 'Winter'});
    title('Optimized single-axis tracking angles (deg)');
end

function save_row_centered_snapshot(var, params, season_name, hour_idx, plane_height, outdir)
    geom = build_row_centered_geometry(var, params);
    weather = params.weather.(season_name);
    sun_alt_deg = rad2deg(weather.beta_s(hour_idx));
    sun_az_deg = wrap_to_180_local(rad2deg(weather.phi_s(hour_idx)));

    if params.tracking_mode == 1
        tilt_deg = rad2deg(var.tracking_angles(season_index(season_name), hour_idx));
    else
        tilt_deg = rad2deg(var.PV_sigma);
    end

    snapshot = build_snapshot(geom, tilt_deg, sun_az_deg, sun_alt_deg, plane_height);
    fig = plot_snapshot(snapshot, geom, season_name, hour_idx, tilt_deg, sun_alt_deg, sun_az_deg);
    filename = sprintf('optimized_shadow_%s_hour%02d_z%.1fm.png', season_name, hour_idx, plane_height);
    saveas(fig, fullfile(outdir, filename));
    close(fig);
end

function gif_file = save_row_centered_shading_gif(var, params, plane_height, delay_seconds, outdir)
    seasons = {'spring', 'summer', 'fall', 'winter'};
    geom = build_row_centered_geometry(var, params);
    gif_file = fullfile(outdir, sprintf('optimized_row_centered_shading_96hr_z%.1fm.gif', plane_height));
    frame_idx = 0;

    for season_idx = 1:numel(seasons)
        season_name = seasons{season_idx};
        weather = params.weather.(season_name);

        for hour_idx = 1:24
            sun_alt_deg = rad2deg(weather.beta_s(hour_idx));
            sun_az_deg = wrap_to_180_local(rad2deg(weather.phi_s(hour_idx)));

            if params.tracking_mode == 1
                tilt_deg = rad2deg(var.tracking_angles(season_idx, hour_idx));
            else
                tilt_deg = rad2deg(var.PV_sigma);
            end

            snapshot = build_snapshot(geom, tilt_deg, sun_az_deg, sun_alt_deg, plane_height);
            fig = plot_snapshot(snapshot, geom, season_name, hour_idx, tilt_deg, sun_alt_deg, sun_az_deg);
            frame_file = fullfile(outdir, sprintf('gif_frame_%03d.png', frame_idx + 1));
            exportgraphics(fig, frame_file, 'Resolution', 120);
            [rgb_data, ~, alpha_data] = imread(frame_file);
            if ~isempty(alpha_data)
                rgb_data = apply_white_background(rgb_data, alpha_data);
            end
            [image_data, color_map] = rgb2ind(rgb_data, 256);

            frame_idx = frame_idx + 1;
            if frame_idx == 1
                imwrite(image_data, color_map, gif_file, 'gif', ...
                    'LoopCount', Inf, 'DelayTime', delay_seconds);
            else
                imwrite(image_data, color_map, gif_file, 'gif', ...
                    'WriteMode', 'append', 'DelayTime', delay_seconds);
            end

            close(fig);
            delete(frame_file);
        end
    end
end

function make_gif = default_make_shading_gif()
    make_gif = false;
end

function rgb_data = apply_white_background(rgb_data, alpha_data)
    alpha = double(alpha_data) ./ 255;
    if ismatrix(alpha)
        alpha = repmat(alpha, 1, 1, 3);
    end
    rgb_data = uint8(double(rgb_data) .* alpha + 255 .* (1 - alpha));
end

function snapshot = build_snapshot(geom, tilt_deg, sun_az_deg, sun_alt_deg, plane_height)
    shadow_unit = repmat(polyshape(), 1, size(geom.panel_centers, 1));
    shadow_hedge = repmat(polyshape(), 1, size(geom.panel_centers, 1));

    for idx = 1:size(geom.panel_centers, 1)
        raw = project_shadow(geom.panel_centers(idx, :)', geom.panel_span, geom.panel_length, ...
            tilt_deg, sun_az_deg, sun_alt_deg, geom.crop_azimuth_deg, plane_height);
        shadow_unit(idx) = intersect(raw, geom.unit_cell);
        shadow_hedge(idx) = intersect(shadow_unit(idx), geom.hedge_shape);
    end

    snapshot.unit_shadow = union(shadow_unit);
    snapshot.hedge_shadow = union(shadow_hedge);
    snapshot.SF_hedge = area(snapshot.hedge_shadow) / geom.hedge_area;
end

function fig = plot_snapshot(snapshot, geom, season_name, hour_idx, tilt_deg, sun_alt_deg, sun_az_deg)
    fig = figure('Name', 'Optimized Row-Centered Shadow Snapshot', 'Visible', 'off', ...
        'Color', 'w', 'Position', [100 100 900 560]);
    hold on;
    h_unit = plot(geom.unit_cell, 'FaceColor', [0.94 0.94 0.94], 'FaceAlpha', 0.18, ...
        'EdgeColor', [0.45 0.45 0.45], 'LineWidth', 1.2);
    h_hedge = plot(geom.hedge_shape, 'FaceColor', [0.55 0.78 0.45], 'FaceAlpha', 0.45, ...
        'EdgeColor', [0.10 0.40 0.10], 'LineWidth', 1.5);
    h_shadow = plot(snapshot.unit_shadow, 'FaceColor', [0.98 0.62 0.20], 'FaceAlpha', 0.35, ...
        'EdgeColor', [0.75 0.35 0.00], 'LineWidth', 1.2);
    h_overlap = plot(snapshot.hedge_shadow, 'FaceColor', [0.05 0.25 0.90], 'FaceAlpha', 0.75, ...
        'EdgeColor', [0.02 0.05 0.40], 'LineWidth', 1.5);
    in_slice = geom.panel_centers(:, 2) >= 0 & geom.panel_centers(:, 2) <= geom.slice_length;
    h_panel = scatter(geom.panel_centers(in_slice, 1), geom.panel_centers(in_slice, 2), 42, ...
        'Marker', 'v', 'MarkerFaceColor', [0.95 0.55 0.12], 'MarkerEdgeColor', [0.35 0.20 0.02]);
    axis equal;
    grid on;
    xlabel('Cross-row position (m)');
    ylabel('Along-row position (m)');
    title(sprintf('Optimized %s hour %02d: tilt %.1f deg, sun alt %.1f deg, sun az %.1f deg, SF %.3f', ...
        season_name, hour_idx, tilt_deg, sun_alt_deg, sun_az_deg, snapshot.SF_hedge));
    legend([h_unit, h_hedge, h_shadow, h_overlap, h_panel], ...
        {'Representative row slice', 'Raspberry hedge', 'Projected shadow in slice', 'Counted hedge shade', 'Panel source centers'}, ...
        'Location', 'eastoutside');
    hold off;
end

function shadow_poly = project_shadow(center_point, panel_span, panel_length, tilt_deg, sun_az_deg, sun_alt_deg, crop_azimuth_deg, plane_height)
    if sind(sun_alt_deg) <= 0
        shadow_poly = polyshape();
        return;
    end

    panel_2D = [
        panel_span / 2,  panel_length / 2, 0
       -panel_span / 2,  panel_length / 2, 0
       -panel_span / 2, -panel_length / 2, 0
        panel_span / 2, -panel_length / 2, 0
    ];

    Rz_crop = [-sind(crop_azimuth_deg), -cosd(crop_azimuth_deg), 0; cosd(crop_azimuth_deg), -sind(crop_azimuth_deg), 0; 0, 0, 1];
    Ry_pv = [cosd(tilt_deg), 0, sind(tilt_deg); 0, 1, 0; -sind(tilt_deg), 0, cosd(tilt_deg)];
    panel_row = repmat(center_point, 1, 4) + Ry_pv * panel_2D';
    panel_global = Rz_crop * panel_row;
    sun_vector = [cosd(sun_az_deg) * cosd(sun_alt_deg), sind(sun_az_deg) * cosd(sun_alt_deg), sind(sun_alt_deg)];
    projected = zeros(4, 3);

    for vertex = 1:4
        t = (plane_height - panel_global(3, vertex)) / sun_vector(3);
        projected(vertex, :) = [
            panel_global(1, vertex) + sun_vector(1) * t, ...
            panel_global(2, vertex) + sun_vector(2) * t, ...
            plane_height
        ];
    end

    projected_ENZ = transpose(Rz_crop) * projected';
    shadow_poly = polyshape(projected_ENZ(1, :), projected_ENZ(2, :));
end

function idx = season_index(season_name)
    season_names = {'spring', 'summer', 'fall', 'winter'};
    idx = find(strcmpi(season_name, season_names), 1);
    if isempty(idx)
        error('Unknown season name: %s', season_name);
    end
end

function wrapped = wrap_to_180_local(angle_deg)
    wrapped = mod(angle_deg + 180, 360) - 180;
end
