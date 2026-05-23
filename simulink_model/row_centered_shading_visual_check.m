%% row_centered_shading_visual_check
% Generate row-centered shadow snapshots under results/. Generated PNGs are
% diagnostics only and should not be committed by default.

base_dir_visual = fileparts(mfilename('fullpath'));
if isempty(base_dir_visual)
    base_dir_visual = pwd;
end
visual_outdir = fullfile(base_dir_visual, 'results', ['row_centered_visual_check_' datestr(now, 'yyyymmdd_HHMMSS')]);
if ~exist(visual_outdir, 'dir')
    mkdir(visual_outdir);
end
setenv('ROW_CENTERED_VISUAL_OUTDIR', visual_outdir);

bdclose('all');
agrivoltaics_variable_definition;
agriParams = configure_visual_case(agriParams, 0);
agriVar = configure_visual_var(agriVar, agriParams);
save_case_snapshot(agriVar, agriParams, 'fixed_axis');

agriParams = configure_visual_case(agriParams, 1);
agriVar = configure_visual_var(agriVar, agriParams);
save_case_snapshot(agriVar, agriParams, 'single_axis');

fprintf('Row-centered visual checks written to %s\n', getenv('ROW_CENTERED_VISUAL_OUTDIR'));

function agriParams = configure_visual_case(agriParams, tracking_mode)
    agriParams.geometry_mode = 1;
    agriParams.tracking_mode = tracking_mode;
end

function agriVar = configure_visual_var(agriVar, agriParams)
    agriVar.PV_l_p = agriParams.fixed_panel_length;
    agriVar.PV_x_p = 0;
    agriVar.PV_y_p = agriParams.row_pitch - agriVar.PV_w_p;

    if agriParams.tracking_mode == 1
        agriVar.PV_phi = agriParams.land_angle;
        agriVar.PV_sigma = 0;
        agriVar.tracking_angles = generate_physics_tracking(agriParams, agriVar);
    else
        agriVar.PV_phi = agriParams.row_centered_fixed_phi;
        if ~isfield(agriVar, 'PV_sigma')
            agriVar.PV_sigma = 0;
        end
        agriVar.tracking_angles = zeros(4, 24);
    end
end

function save_case_snapshot(agriVar, agriParams, label)
    season_name = 'summer';
    hour_idx = 10;
    plane_height = 0;
    geom = build_row_centered_geometry(agriVar, agriParams);
    weather = agriParams.weather.(season_name);
    sun_alt_deg = rad2deg(weather.beta_s(hour_idx));
    sun_az_deg = wrap_to_180_local(rad2deg(weather.phi_s(hour_idx)));

    if agriParams.tracking_mode == 1
        tilt_deg = rad2deg(agriVar.tracking_angles(2, hour_idx));
    else
        tilt_deg = rad2deg(agriVar.PV_sigma);
    end

    snapshot = build_snapshot(geom, tilt_deg, sun_az_deg, sun_alt_deg, plane_height);
    fig = plot_snapshot(snapshot, geom, label, season_name, hour_idx, tilt_deg, sun_alt_deg, sun_az_deg);
    outdir = getenv('ROW_CENTERED_VISUAL_OUTDIR');
    saveas(fig, fullfile(outdir, sprintf('%s_%s_hour%02d_z%.1fm.png', label, season_name, hour_idx, plane_height)));
    close(fig);
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

function fig = plot_snapshot(snapshot, geom, label, season_name, hour_idx, tilt_deg, sun_alt_deg, sun_az_deg)
    fig = figure('Name', 'Row-centered shadow snapshot', 'Visible', 'off', 'Color', 'w', 'Position', [100 100 900 560]);
    hold on;
    h_unit = plot(geom.unit_cell, 'FaceColor', [0.94 0.94 0.94], 'FaceAlpha', 0.18, 'EdgeColor', [0.45 0.45 0.45], 'LineWidth', 1.2);
    h_hedge = plot(geom.hedge_shape, 'FaceColor', [0.55 0.78 0.45], 'FaceAlpha', 0.45, 'EdgeColor', [0.10 0.40 0.10], 'LineWidth', 1.5);
    h_shadow = plot(snapshot.unit_shadow, 'FaceColor', [0.98 0.62 0.20], 'FaceAlpha', 0.35, 'EdgeColor', [0.75 0.35 0.00], 'LineWidth', 1.2);
    h_overlap = plot(snapshot.hedge_shadow, 'FaceColor', [0.05 0.25 0.90], 'FaceAlpha', 0.75, 'EdgeColor', [0.02 0.05 0.40], 'LineWidth', 1.5);
    in_slice = geom.panel_centers(:, 2) >= 0 & geom.panel_centers(:, 2) <= geom.slice_length;
    h_panel = scatter(geom.panel_centers(in_slice, 1), geom.panel_centers(in_slice, 2), 42, ...
        'Marker', 'v', 'MarkerFaceColor', [0.95 0.55 0.12], 'MarkerEdgeColor', [0.35 0.20 0.02]);
    axis equal;
    grid on;
    xlabel('Cross-row position (m)');
    ylabel('Along-row position (m)');
    title(sprintf('%s %s hour %02d: tilt %.1f deg, sun alt %.1f deg, sun az %.1f deg, SF %.3f', ...
        strrep(label, '_', ' '), season_name, hour_idx, tilt_deg, sun_alt_deg, sun_az_deg, snapshot.SF_hedge));
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

function wrapped = wrap_to_180_local(angle_deg)
    wrapped = mod(angle_deg + 180, 360) - 180;
end
