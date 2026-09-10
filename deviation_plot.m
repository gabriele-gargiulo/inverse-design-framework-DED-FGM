%% DEVIATION PLOT
% Compares two 3D composition files and plots the reference composition,
% evaluated composition, and pointwise deviation or correction magnitude.

%% INPUT FILES

reference_file = 'input_material_distribution_3d_blisk.csv'; % reference composition
evaluated_file = 'output_material_distribution_3d_free.csv'; % composition to evaluate

data_ref = readtable(reference_file);
data_eval = readtable(evaluated_file);

%% PROCESS PARAMETERS

h = 0.62; % height

%% USER OPTIONS

plot_mode = "deviation";       % "deviation" or "correction"
deviation_error_bar_max = 67.99; % positive number or "max"
correction_error_bar_max = "max"; % positive number or "max"

section = false;          % show Y-Z section plane
section_y_range = [-30, 30];
section_z_range = [0, 120];

three_quarter_blisk = false; % hide one quarter of the blisk

color_IN  = [123, 132, 181]/255; % IN plot color
color_MK  = [240, 121, 121]/255; % MK plot color
color_INC = [253, 243, 153]/255; % INC plot color

coordinate_mode = "cartesian";
cylinder_axis = 'Y';
shaft_center = [0, 0];      % center in the plane normal to cylinder_axis
theta_zero = 0;             % angular seam rotation [rad]
reference_radius_mode = "inner";
full_circumference = true;  % wrap theta across 0 and 2*pi

%% SELECT PLOT LABELS

switch plot_mode
    case "deviation"
        title_ref = 'reference composition';
        title_eval = 'output composition';
        title_error = 'deviation (%)';
        colorbar_label = 'Deviation (%)';
        error_bar_max = deviation_error_bar_max;

    case "correction"
        title_ref = 'reference composition';
        title_eval = 'corrected input composition';
        title_error = 'correction (%)';
        colorbar_label = 'Correction (%)';
        error_bar_max = correction_error_bar_max;

    otherwise
        error("Unknown plot_mode. Use 'deviation' or 'correction'.");
end

%% EXTRACT DATA

X_eval = data_eval.X;
Y_eval = data_eval.Y;
Z_eval = data_eval.Z;

ref_points = [data_ref.X, data_ref.Y, data_ref.Z];
eval_points = [X_eval, Y_eval, Z_eval];

C_ref_all = [data_ref.IN_frac, data_ref.MK_frac, data_ref.INC_frac];
C_eval = [data_eval.IN_frac, data_eval.MK_frac, data_eval.INC_frac];

if strcmpi(coordinate_mode, "cylindrical")
    [A_ref, Theta_ref, R_ref_points] = errorPlotCartesianToCylindrical( ...
        data_ref.X, data_ref.Y, data_ref.Z, cylinder_axis, shaft_center, theta_zero);
    [A_eval, Theta_eval, R_eval_points] = errorPlotCartesianToCylindrical( ...
        X_eval, Y_eval, Z_eval, cylinder_axis, shaft_center, theta_zero);

    if strcmpi(reference_radius_mode, "inner")
        R_period_ref = min(R_eval_points);
    else
        R_period_ref = mean(R_eval_points);
    end

    ref_points_model = [A_ref, R_period_ref*Theta_ref, R_ref_points];
    eval_points_model = [A_eval, R_period_ref*Theta_eval, R_eval_points];
else
    R_period_ref = NaN;
    ref_points_model = ref_points;
    eval_points_model = eval_points;
end

idx_ref = nearest_reference_point_by_layer( ...
    ref_points_model, eval_points_model, h, ...
    coordinate_mode, R_period_ref, full_circumference);
C_ref = C_ref_all(idx_ref,:);

%% NORMALIZE COMPOSITION

row_sums_ref = sum(C_ref,2);
row_sums_eval = sum(C_eval,2);

row_sums_ref(row_sums_ref == 0) = 1;
row_sums_eval(row_sums_eval == 0) = 1;

C_ref = C_ref ./ row_sums_ref;
C_eval = C_eval ./ row_sums_eval;

%% BUILD PLOT COLORS

RGB_ref = C_ref(:,1).*color_IN + ...
          C_ref(:,2).*color_MK + ...
          C_ref(:,3).*color_INC;

RGB_eval = C_eval(:,1).*color_IN + ...
           C_eval(:,2).*color_MK + ...
           C_eval(:,3).*color_INC;

points = eval_points;
axis_points = include_section_axis_points( ...
    points, section, section_y_range, section_z_range);
plot_mask = make_three_quarter_blisk_mask(X_eval, Z_eval, three_quarter_blisk);

%% PLOT REFERENCE COMPOSITION

figure; hold on;
scatter3(X_eval(plot_mask), Y_eval(plot_mask), Z_eval(plot_mask), ...
    25, RGB_ref(plot_mask,:), 'filled');
add_yz_section_patch(gca, section, section_y_range, section_z_range);
axis equal; view(3);
set_point_axis_limits(gca, axis_points);
ax = gca;
ax.YDir = 'reverse';
add_material_legend(color_IN, color_MK, color_INC);
zackificator_3d(gca, title_ref);

%% PLOT EVALUATED COMPOSITION

figure; hold on;
scatter3(X_eval(plot_mask), Y_eval(plot_mask), Z_eval(plot_mask), ...
    25, RGB_eval(plot_mask,:), 'filled');
add_yz_section_patch(gca, section, section_y_range, section_z_range);
axis equal; view(3);
set_point_axis_limits(gca, axis_points);
ax = gca;
ax.YDir = 'reverse';
add_material_legend(color_IN, color_MK, color_INC);
zackificator_3d(gca, title_eval);

%% COMPUTE DEVIATION

error_vec = C_eval - C_ref;
error_mag = sqrt(sum(error_vec.^2, 2));

% Convert to percentage of maximum possible error.
error_percent = 100 * error_mag / sqrt(2);

sq_error = sum(error_vec.^2, 2);
rmse = sqrt(mean(sq_error));
error_percent_rmse = 100 * rmse / sqrt(2);

%% PLOT DEVIATION MAP

figure; hold on;

cmap = error_plot_colormap(plot_mode);

if isstring(error_bar_max) || ischar(error_bar_max)
    if strcmpi(string(error_bar_max), "max")
        cmax = max(error_percent);
    else
        error('error_bar_max must be a positive number or "max".');
    end
elseif isnumeric(error_bar_max) && isscalar(error_bar_max) && error_bar_max > 0
    cmax = error_bar_max;
else
    error('error_bar_max must be a positive number or "max".');
end
cmin = 0;

scatter3(X_eval(plot_mask), Y_eval(plot_mask), Z_eval(plot_mask), ...
    25, error_percent(plot_mask), 'filled');
add_yz_section_patch(gca, section, section_y_range, section_z_range);
axis equal; view(3);
set_point_axis_limits(gca, axis_points);
colormap(cmap);
clim([cmin cmax]);
ax = gca;
ax.YDir = 'reverse';
cb = colorbar;
cb.Label.String = colorbar_label;
cb.Label.Interpreter = 'none';
zackificator_3d(gca, title_error);

%% PRINT METRICS

fprintf('Mean error: %.2f %%\n', mean(error_percent));
fprintf('Max error : %.2f %%\n', max(error_percent));
fprintf('Min error : %.2f %%\n', min(error_percent));
fprintf('RMSE: %.2f %%\n', error_percent_rmse);

%% LOCAL FUNCTIONS

function plot_mask = make_three_quarter_blisk_mask(X, Z, three_quarter_blisk)

plot_mask = true(size(X));

if ~three_quarter_blisk
    return
end

center_x = 0.5*(min(X) + max(X));
center_z = 0.5*(min(Z) + max(Z));

theta = atan2(Z - center_z, X - center_x);

% Remove the fourth quarter in the standard plotted view.
removed_quarter = theta >= pi/2 & theta <= pi;
plot_mask = ~removed_quarter;

end

function cmap = error_plot_colormap(plot_mode)

if strcmpi(string(plot_mode), "correction")
    cmap = hot(256);
else
    cmap = jet(256);
end

end

function add_yz_section_patch(ax, section, y_range, z_range)

if ~section
    return
end

x_section = 0;
y0 = y_range(1);
y1 = y_range(2);
z0 = z_range(1);
z1 = z_range(2);

patch(ax, ...
    x_section*[1 1 1 1], ...
    [y0 y1 y1 y0], ...
    [z0 z0 z1 z1], ...
    [0.15 0.15 0.15], ...
    'FaceAlpha', 0.16, ...
    'EdgeColor', [0.05 0.05 0.05], ...
    'LineWidth', 0.8, ...
    'HandleVisibility', 'off');

end

function axis_points = include_section_axis_points(points, section, y_range, z_range)

axis_points = points;

if ~section
    return
end

section_points = [ ...
    0, y_range(1), z_range(1); ...
    0, y_range(2), z_range(1); ...
    0, y_range(2), z_range(2); ...
    0, y_range(1), z_range(2)];

axis_points = [axis_points; section_points];

end

function idx_ref = nearest_reference_point_by_layer( ...
    ref_points_model, eval_points_model, h, coordinate_mode, reference_radius, full_circumference)

build_coord_eval = eval_points_model(:,3);
layer_of_eval = zeros(size(build_coord_eval));

if isempty(build_coord_eval)
    idx_ref = zeros(0,1);
    return
end

current_layer = 1;
layer_of_eval(1) = current_layer;

for q = 2:numel(build_coord_eval)
    if build_coord_eval(q) - build_coord_eval(q-1) > h/2
        current_layer = current_layer + 1;
    end

    layer_of_eval(q) = current_layer;
end

n_layers = current_layer;
ref_build_coord = ref_points_model(:,3);
idx_ref = zeros(size(build_coord_eval));
coord_tol = 1e-9;

for k = 1:n_layers
    idx_eval_layer = find(layer_of_eval == k);
    build_coord_layer = build_coord_eval(idx_eval_layer);
    layer_min = min(build_coord_layer) - coord_tol;
    layer_max = max(build_coord_layer) + coord_tol;
    idx_ref_layer = find( ...
        ref_build_coord >= layer_min & ...
        ref_build_coord <= layer_max);

    if isempty(idx_ref_layer)
        layer_center = mean(build_coord_layer);
        idx_ref_layer = find(abs(ref_build_coord - layer_center) <= h/2 + coord_tol);
    end

    if isempty(idx_ref_layer)
        error('No reference material points found on evaluated layer %d.', k);
    end

    if strcmpi(coordinate_mode, "cylindrical")
        [idx_local, ~] = nearest_cylindrical_points( ...
            ref_points_model(idx_ref_layer,:), ...
            eval_points_model(idx_eval_layer,:), ...
            reference_radius, full_circumference);
    else
        [idx_local, ~] = knnsearch( ...
            ref_points_model(idx_ref_layer,:), ...
            eval_points_model(idx_eval_layer,:));
    end

    idx_ref(idx_eval_layer) = idx_ref_layer(idx_local);
end

end

function [A, Theta, R] = errorPlotCartesianToCylindrical(X, Y, Z, cylinder_axis, shaft_center, theta_zero)

switch upper(cylinder_axis)
    case 'X'
        A = X;
        U = Y - shaft_center(1);
        V = Z - shaft_center(2);
    case 'Y'
        A = Y;
        U = X - shaft_center(1);
        V = Z - shaft_center(2);
    case 'Z'
        A = Z;
        U = X - shaft_center(1);
        V = Y - shaft_center(2);
    otherwise
        error('cylinder_axis must be ''X'', ''Y'', or ''Z''.');
end

Theta = mod(atan2(V, U) - theta_zero, 2*pi);
R = hypot(U, V);

end

function [idx_nearest, dist_nearest] = nearest_cylindrical_points( ...
    coords_ref_model, coords_query_model, reference_radius, full_circumference)

if full_circumference
    period = 2*pi*reference_radius;
    n_ref = size(coords_ref_model,1);
    coords_aug = [coords_ref_model; ...
                  coords_ref_model + [zeros(n_ref,1), period*ones(n_ref,1), zeros(n_ref,1)]; ...
                  coords_ref_model - [zeros(n_ref,1), period*ones(n_ref,1), zeros(n_ref,1)]];
    idx_aug_to_ref = [(1:n_ref)'; (1:n_ref)'; (1:n_ref)'];

    [idx_aug, dist_nearest] = knnsearch(coords_aug, coords_query_model);
    idx_nearest = idx_aug_to_ref(idx_aug);
else
    [idx_nearest, dist_nearest] = knnsearch(coords_ref_model, coords_query_model);
end

end

function set_point_axis_limits(ax, points)

margin = 0.02 * max([ ...
    max(points(:,1)) - min(points(:,1)), ...
    max(points(:,2)) - min(points(:,2)), ...
    max(points(:,3)) - min(points(:,3)), ...
    eps]);

xlim(ax, [min(points(:,1))-margin, max(points(:,1))+margin]);
ylim(ax, [min(points(:,2))-margin, max(points(:,2))+margin]);
zlim(ax, [min(points(:,3))-margin, max(points(:,3))+margin]);

end

function add_material_legend(color_IN, color_MK, color_INC)

h1 = scatter3(NaN, NaN, NaN, 55, color_IN,  'filled');
h2 = scatter3(NaN, NaN, NaN, 55, color_MK,  'filled');
h3 = scatter3(NaN, NaN, NaN, 55, color_INC, 'filled');

legend([h1 h2 h3], {'IN','MK','INC'});

end

function zackificator_3d(ax, ~)

set(gcf,'units','centimeters','position',[10,10,16,12]);
set(gcf,'color','w');

xlabel(ax,'');
ylabel(ax,'');
zlabel(ax,'');
title(ax,'');
axis(ax,'off');

ax.FontSize = 9;
ax.FontName = 'arial';
ax.LineWidth = 1;
ax.TickDir = 'in';
ax.TickLength = [0.015 0.025];
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'off';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.Title.FontName = 'arial';
ax.Title.FontSize = 10;
ax.Title.FontWeight = 'bold';
ax.XLabel.FontName = 'arial';
ax.XLabel.FontSize = 10;
ax.XLabel.FontWeight = 'normal';
ax.YLabel.FontName = 'arial';
ax.YLabel.FontSize = 10;
ax.YLabel.FontWeight = 'normal';
ax.ZLabel.FontName = 'arial';
ax.ZLabel.FontSize = 10;
ax.ZLabel.FontWeight = 'normal';

h_legend = findobj(gcf,'Type','Legend');
for lgd = reshape(h_legend,1,[])
    lgd.Location = 'northeastoutside';
    lgd.FontSize = 9;
    lgd.FontName = 'arial';
    lgd.Box = 'off';
end

h_colorbar = findall(gcf,'Type','ColorBar');
for cb = reshape(h_colorbar,1,[])
    cb.FontSize = 11;
    cb.FontName = 'arial';
    cb.LineWidth = 1;
    cb.TickDirection = 'in';
    cb.Label.FontSize = 13;
    cb.Label.FontName = 'arial';
    cb.Label.FontWeight = 'normal';
end

end
