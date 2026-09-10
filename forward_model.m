%clear; clc;

%% FORWARD 3D MIXING MODEL
% Computes the deposited composition along a 3D toolpath. The
% feed composition is sampled from the nearest material point, then updated
% with geometric remelting and in-track dragging.

%% INPUT FILES

input_material_file = 'input_material_distribution_3d_blisk.csv'; % feed composition points
toolpath_file = 'toolpath_3d_blisk.csv';                          % deposition toolpath

mat_data = readtable(input_material_file);
toolpath = readtable(toolpath_file);

%% PROCESS PARAMETERS

h = 0.62;          % height
w = 1.97;          % width
b = 0.21;          % penetration depth

% Point and hatch spacing are read directly from the toolpath coordinates.

d_spot = 1.6;      % laser spot diameter
L_melt = 1.5*d_spot; % melt pool length

%% USER OPTIONS

% Numerical resolution of the latest-bead-wins volume partition.
overlap_partition_samples = 16384;

% "cartesian" treats Z as the build direction.
% "cylindrical" treats radius from the shaft as the build direction.
coordinate_mode = "cartesian";
cylinder_axis = 'Y';
shaft_center = [0, 0];      % center in the plane normal to cylinder_axis
theta_zero = 0;             % angular seam rotation [rad]
reference_radius_mode = "inner";
full_circumference = true;  % wrap theta across 0 and 2*pi
track_break_distance = 2*w;

%% EXTRACT POINTS AND FEED COMPOSITION

X_mat = mat_data.X;
Y_mat = mat_data.Y;
Z_mat = mat_data.Z;

IN_mat = mat_data.IN_frac;
MK_mat = mat_data.MK_frac;

if ismember('INC_frac', mat_data.Properties.VariableNames)
    INC_mat = mat_data.INC_frac;
else
    INC_mat = zeros(height(mat_data),1);
end

X_path = toolpath.X;
Y_path = toolpath.Y;
Z_path = toolpath.Z;
path_points_real = [X_path, Y_path, Z_path];
material_points_real = [X_mat, Y_mat, Z_mat];

%% POINT CLOUD

free_points = path_points_real;
N = size(free_points,1);

if strcmpi(coordinate_mode, "cylindrical")

    [A_mat, Theta_mat, R_mat] = freeCartesianToCylindrical( ...
        X_mat, Y_mat, Z_mat, cylinder_axis, shaft_center, theta_zero);
    [A_path, Theta_path, R_path] = freeCartesianToCylindrical( ...
        X_path, Y_path, Z_path, cylinder_axis, shaft_center, theta_zero);

    if strcmpi(reference_radius_mode, "inner")
        R_ref = min(R_path);
    else
        R_ref = mean(R_path);
    end

    if R_ref <= 0
        error('reference radius must be positive. Check shaft_center and cylindrical coordinates.');
    end

    material_points_model = [A_mat, R_ref*Theta_mat, R_mat];
    coords_path = [A_path, R_ref*Theta_path, R_path];
    free_points_model = coords_path;

else

    R_ref = NaN;
    material_points_model = material_points_real;
    coords_path = [X_path, Y_path, Z_path];
    free_points_model = coords_path;
end

X_path = coords_path(:,1);
Y_path = coords_path(:,2);
Z_path = coords_path(:,3);

idx_toolpath = (1:N)';
dist_tool = zeros(N,1);

%% TOOLPATH SCAN STRATEGY

coord_tol = 1e-9;

build_coord_path = Z_path;
layer_of_path_point = zeros(length(Z_path),1);

if isempty(build_coord_path)
    n_layers_path = 0;
else
    current_layer = 1;
    layer_of_path_point(1) = current_layer;

    for q = 2:length(build_coord_path)
        if build_coord_path(q) - build_coord_path(q-1) > h/2
            current_layer = current_layer + 1;
        end

        layer_of_path_point(q) = current_layer;
    end

    n_layers_path = current_layer;
end

material_build_coord = material_points_model(:,3);
idx_feed_material = zeros(N,1);

for k = 1:n_layers_path
    idx_path_layer = find(layer_of_path_point == k);
    build_coord_layer = build_coord_path(idx_path_layer);
    layer_min = min(build_coord_layer) - coord_tol;
    layer_max = max(build_coord_layer) + coord_tol;
    idx_material_layer = find( ...
        material_build_coord >= layer_min & ...
        material_build_coord <= layer_max);

    if isempty(idx_material_layer)
        layer_center = mean(build_coord_layer);
        idx_material_layer = find(abs(material_build_coord - layer_center) <= h/2 + coord_tol);
    end

    if isempty(idx_material_layer)
        error('No feed material points found on toolpath layer %d.', k);
    end

    if strcmpi(coordinate_mode, "cylindrical")
        [idx_local, ~] = freeNearestToolpathCylindrical( ...
            material_points_model(idx_material_layer,:), ...
            free_points_model(idx_path_layer,:), ...
            R_ref, full_circumference);
    else
        [idx_local, ~] = knnsearch( ...
            material_points_model(idx_material_layer,:), ...
            free_points_model(idx_path_layer,:));
    end

    idx_feed_material(idx_path_layer) = idx_material_layer(idx_local);
end

C_material_feed = [IN_mat, MK_mat, INC_mat];
C_feed = C_material_feed(idx_feed_material,:);

scan_axis_path = ones(n_layers_path,1);   % 1 = tracks along X, 2 = tracks along Y

for k = 1:n_layers_path

    idx_layer = find(layer_of_path_point == k);
    xy_layer = [X_path(idx_layer), Y_path(idx_layer)];

    if size(xy_layer,1) > 1
        dx_total = sum(abs(diff(xy_layer(:,1))));
        dy_total = sum(abs(diff(xy_layer(:,2))));
    else
        dx_total = 0;
        dy_total = 0;
    end

    if dy_total > dx_total
        scan_axis_path(k) = 2;
    else
        scan_axis_path(k) = 1;
    end

end

% Segment the toolpath into physical tracks.
path_track_id = zeros(length(Z_path),1);
track_layer = zeros(length(Z_path),1);
track_scan_axis = zeros(length(Z_path),1);
track_scan_dir = zeros(length(Z_path),1);
track_hatch_dir = zeros(length(Z_path),1);
track_hatch_coord = zeros(length(Z_path),1);
next_track_id = 0;

for k = 1:n_layers_path

    idx_layer = find(layer_of_path_point == k);

    if isempty(idx_layer)
        continue
    end

    scan_axis = scan_axis_path(k);
    hatch_axis = 3 - scan_axis;
    xy_layer = [X_path(idx_layer), Y_path(idx_layer)];
    hatch_values_layer = xy_layer(:,hatch_axis);

    track_start_local = 1;

    for q = 2:length(idx_layer)
        toolpath_jump = norm( ...
            path_points_real(idx_layer(q),:) - path_points_real(idx_layer(q-1),:));

        if abs(hatch_values_layer(q) - hatch_values_layer(q-1)) > coord_tol || ...
           toolpath_jump > track_break_distance

            idx_track_path = idx_layer(track_start_local:q-1);
            next_track_id = next_track_id + 1;
            path_track_id(idx_track_path) = next_track_id;

            scan_values_track = [X_path(idx_track_path), Y_path(idx_track_path)];
            scan_values_track = scan_values_track(:,scan_axis);
            ds_track = diff(scan_values_track);
            ds_track = ds_track(abs(ds_track) > coord_tol);

            if ~isempty(ds_track)
                scan_dir = sign(ds_track(1));
            elseif numel(scan_values_track) > 1 && ...
                   abs(scan_values_track(end)-scan_values_track(1)) > coord_tol
                scan_dir = sign(scan_values_track(end)-scan_values_track(1));
            else
                scan_dir = 1;
            end

            track_layer(next_track_id,1) = k;
            track_scan_axis(next_track_id,1) = scan_axis;
            track_scan_dir(next_track_id,1) = scan_dir;
            track_hatch_coord(next_track_id,1) = hatch_values_layer(track_start_local);

            track_start_local = q;
        end
    end

    idx_track_path = idx_layer(track_start_local:end);
    next_track_id = next_track_id + 1;
    path_track_id(idx_track_path) = next_track_id;

    scan_values_track = [X_path(idx_track_path), Y_path(idx_track_path)];
    scan_values_track = scan_values_track(:,scan_axis);
    ds_track = diff(scan_values_track);
    ds_track = ds_track(abs(ds_track) > coord_tol);

    if ~isempty(ds_track)
        scan_dir = sign(ds_track(1));
    elseif numel(scan_values_track) > 1 && ...
           abs(scan_values_track(end)-scan_values_track(1)) > coord_tol
        scan_dir = sign(scan_values_track(end)-scan_values_track(1));
    else
        scan_dir = 1;
    end

    track_layer(next_track_id,1) = k;
    track_scan_axis(next_track_id,1) = scan_axis;
    track_scan_dir(next_track_id,1) = scan_dir;
    track_hatch_coord(next_track_id,1) = hatch_values_layer(track_start_local);
end

for k = 1:n_layers_path

    tracks_in_layer = find(track_layer == k);

    for q = 1:length(tracks_in_layer)

        t = tracks_in_layer(q);

        if q < length(tracks_in_layer)
            hatch_step = track_hatch_coord(tracks_in_layer(q+1)) - track_hatch_coord(t);
        elseif q > 1
            hatch_step = track_hatch_coord(t) - track_hatch_coord(tracks_in_layer(q-1));
        else
            hatch_step = 1;
        end

        if abs(hatch_step) > coord_tol
            track_hatch_dir(t,1) = sign(hatch_step);
        else
            track_hatch_dir(t,1) = 1;
        end
    end
end

%% ASSIGN POINTS TO TOOLPATH LAYERS AND TRACKS

layer_of_free_point = layer_of_path_point(idx_toolpath);

track_id = path_track_id(idx_toolpath);

unassigned_track = track_id == 0;
if any(unassigned_track)
    track_id(unassigned_track) = 1;
end

%% BUILD TRUE PRINT ORDER FROM TOOLPATH

sort_matrix = [idx_toolpath(:), dist_tool(:)];
[~, sort_idx] = sortrows(sort_matrix, [1 2]);

free_points_print = free_points_model(sort_idx,:);
free_points_real_print = free_points(sort_idx,:);
C_feed_print = C_feed(sort_idx,:);
track_id_print = track_id(sort_idx);
layer_print = layer_of_free_point(sort_idx);
scan_axis_print = track_scan_axis(track_id_print);
scan_dir_print = track_scan_dir(track_id_print);
hatch_dir_print = track_hatch_dir(track_id_print);

inverse_sort_idx = zeros(N,1);
inverse_sort_idx(sort_idx) = 1:N; 

track_ids_unique = unique(track_id_print, 'stable');
track_indices = cell(numel(track_ids_unique),1);
track_group_print = zeros(N,1);
position_in_track = zeros(N,1);
track_distance_print = zeros(N,1);

for q = 1:numel(track_ids_unique)
    idx_track = find(track_id_print == track_ids_unique(q));
    track_indices{q} = idx_track;
    track_group_print(idx_track) = q;
    position_in_track(idx_track) = 1:numel(idx_track);

    P_track = free_points_real_print(idx_track,:);
    if numel(idx_track) > 1
        segment_lengths = sqrt(sum(diff(P_track,1,1).^2,2));
        track_distance_print(idx_track) = [0; cumsum(segment_lengths)];
    end
end

%% TOOLPATH-DRIVEN POINT MIXING

C_real_print = zeros(N,3);
C_standard_print = zeros(N,3);

idx_lat_print = zeros(N,1);
idx_vert_print = zeros(N,1);
alpha_lat_print = zeros(N,1);
alpha_vert_print = zeros(N,1);
n_lat_neighbors = zeros(N,1);
n_vert_neighbors = zeros(N,1);

ellipsoid_width = w;
ellipsoid_height = h + b;
neighbor_search_radius = sqrt(ellipsoid_width^2 + ellipsoid_height^2);

if exist('rangesearch','file') == 2
    neighbor_candidates_all = freeRangeSearchCandidates( ...
        free_points_print, neighbor_search_radius, ...
        coordinate_mode, R_ref, full_circumference);
else
    neighbor_candidates_all = {};
end

for i = 1:N

    point_curr = free_points_print(i,:);
    C_f = C_feed_print(i,:);

    % Only already-printed points can remelt the current point.
    if isempty(neighbor_candidates_all)
        prev_idx = 1:i-1;
    else
        prev_idx = neighbor_candidates_all{i};
        prev_idx = prev_idx(prev_idx < i);
    end

    % Dense points on the current continuous track represent the same
    % moving meltpool and are handled by dragging, not by remelting.
    if ~isempty(prev_idx)
        same_track = track_id_print(prev_idx) == track_id_print(i);
        distance_behind = track_distance_print(i) - track_distance_print(prev_idx);
        prev_idx(same_track & distance_behind < w) = [];
    end

    C_standard_print(i,:) = C_f;

    if ~isempty(prev_idx)

        if strcmpi(coordinate_mode, "cylindrical")
            [dxy, dz_signed, center_xy] = freeCylindricalLocalDistances( ...
                free_points_real_print(prev_idx,:), ...
                free_points_real_print(i,:), ...
                cylinder_axis, shaft_center);
        else
            delta = free_points_print(prev_idx,:) - point_curr;
            center_xy = delta(:,1:2);
            dxy = hypot(delta(:,1), delta(:,2));
            dz_signed = delta(:,3);
        end

        meltpool_intersects = ...
            dxy <= ellipsoid_width + coord_tol & ...
            freeBeadVerticalRangesOverlap(dz_signed, h, b, coord_tol);

        same_layer = layer_print(prev_idx) == layer_print(i);
        lower_layer = layer_print(prev_idx) < layer_print(i);
        eligible = meltpool_intersects & (same_layer | lower_layer);

        if any(eligible)
            candidates = prev_idx(eligible);
            candidate_centers = center_xy(eligible,:);
            candidate_dz = dz_signed(eligible);

            % Chronological order is required by the latest-bead-wins rule.
            [candidates, order] = sort(candidates);
            candidate_centers = candidate_centers(order,:);
            candidate_dz = candidate_dz(order);

            weights = freeExclusiveOverlapFractions( ...
                candidate_centers, candidate_dz, w, h, b, ...
                overlap_partition_samples);
            positive = weights > 0;

            candidates = candidates(positive);
            weights = weights(positive);

            if ~isempty(candidates)
                alpha_total = sum(weights);
                fresh_fraction = max(0, 1-alpha_total);
                C_standard_print(i,:) = ...
                    fresh_fraction*C_f + ...
                    sum(C_real_print(candidates,:) .* weights, 1);

                candidate_layers = layer_print(candidates);
                lat_mask = candidate_layers == layer_print(i);
                vert_mask = candidate_layers < layer_print(i);

                alpha_lat_print(i) = sum(weights(lat_mask));
                alpha_vert_print(i) = sum(weights(vert_mask));
                n_lat_neighbors(i) = sum(lat_mask);
                n_vert_neighbors(i) = sum(vert_mask);

                if any(lat_mask)
                    idx_lat_print(i) = candidates(find(lat_mask,1,'last'));
                end
                if any(vert_mask)
                    idx_vert_print(i) = candidates(find(vert_mask,1,'last'));
                end
            end
        end
    end

    track_idx = track_indices{track_group_print(i)};
    kk = position_in_track(i);
    idx_track_so_far = track_idx(1:kk);
    C_real_print(i,:) = freeApplyTrackDraggingAtPoint( ...
        C_standard_print(idx_track_so_far,:), ...
        free_points_real_print(idx_track_so_far,:), ...
        L_melt);
end

%% RESTORE ORIGINAL INPUT ORDER

C_real = zeros(size(C_real_print));
C_real(sort_idx,:) = C_real_print;

%% SAVE OUTPUT

X_free = free_points(:,1);
Y_free = free_points(:,2);
Z_free = free_points(:,3);

output_table = table( ...
    X_free, Y_free, Z_free, ...
    round(C_real(:,1),3), round(C_real(:,2),3), round(C_real(:,3),3), ...
    'VariableNames', ...
    {'X','Y','Z','IN_frac','MK_frac','INC_frac'});

writetable(output_table, 'output_material_distribution_3d_free.csv');

disp('Free-point 3D material distribution saved.');

%% LOCAL FUNCTIONS

function fractions = freeExclusiveOverlapFractions( ...
    center_xy, dz_signed, w, h, b, n_samples)

center_xy = reshape(center_xy, [], 2);
dz_signed = dz_signed(:);
n_candidates = size(center_xy,1);
fractions = zeros(n_candidates,1);

if n_candidates == 0
    return
end

unit_points = freeUnitBallSamples(n_samples);
a = w/2;
sample_xyz = zeros(size(unit_points));
sample_xyz(:,1:2) = a*unit_points(:,1:2);

top_half = unit_points(:,3) >= 0;
sample_xyz(top_half,3) = h*unit_points(top_half,3);
sample_xyz(~top_half,3) = b*unit_points(~top_half,3);

% The asymmetric z scaling has a different Jacobian in each half.
sample_weights = b*ones(size(unit_points,1),1);
sample_weights(top_half) = h;
total_weight = sum(sample_weights);
unclaimed = true(size(sample_weights));

% Assign every sample to at most one candidate, newest first.
for j = n_candidates:-1:1
    dx = (sample_xyz(:,1)-center_xy(j,1))/a;
    dy = (sample_xyz(:,2)-center_xy(j,2))/a;
    z_relative = sample_xyz(:,3)-dz_signed(j);
    z_scale = b*ones(size(z_relative));
    z_scale(z_relative >= 0) = h;
    inside = ...
        dx.^2 + dy.^2 + (z_relative./z_scale).^2 <= 1 + 10*eps;
    claimed = inside & unclaimed;

    fractions(j) = sum(sample_weights(claimed)) / total_weight;
    unclaimed(claimed) = false;
end

fraction_sum = sum(fractions);
if fraction_sum > 1
    fractions = fractions / fraction_sum;
end

end

function points = freeUnitBallSamples(n_samples)

persistent cached_count cached_points

n_samples = max(1024, round(n_samples));
if ~isempty(cached_count) && cached_count == n_samples
    points = cached_points;
    return
end

sample_index = (1:n_samples)';
u_radius = freeRadicalInverse(sample_index, 2);
u_polar = freeRadicalInverse(sample_index, 3);
u_azimuth = freeRadicalInverse(sample_index, 5);

radius = u_radius.^(1/3);
cos_polar = 2*u_polar-1;
sin_polar = sqrt(max(0, 1-cos_polar.^2));
azimuth = 2*pi*u_azimuth;

points = [ ...
    radius.*sin_polar.*cos(azimuth), ...
    radius.*sin_polar.*sin(azimuth), ...
    radius.*cos_polar];

cached_count = n_samples;
cached_points = points;

end

function values = freeRadicalInverse(sample_index, base)

values = zeros(size(sample_index));
factor = 1/base;

while any(sample_index > 0)
    values = values + factor*mod(sample_index,base);
    sample_index = floor(sample_index/base);
    factor = factor/base;
end

end

function has_overlap = freeBeadVerticalRangesOverlap(dz_signed, h, b, coord_tol)

previous_bottom = dz_signed - b;
previous_top = dz_signed + h;
current_bottom = -b;
current_top = h;

z0 = max(current_bottom, previous_bottom);
z1 = min(current_top, previous_top);
has_overlap = z1 >= z0 - coord_tol;

end

function [A, Theta, R] = freeCartesianToCylindrical(X, Y, Z, cylinder_axis, shaft_center, theta_zero)

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

function [idx_nearest, dist_nearest] = freeNearestToolpathCylindrical( ...
    coords_path_model, coords_free_model, reference_radius, full_circumference)

if full_circumference
    period = 2*pi*reference_radius;
    n_path = size(coords_path_model,1);
    coords_aug = [coords_path_model; ...
                  coords_path_model + [zeros(n_path,1), period*ones(n_path,1), zeros(n_path,1)]; ...
                  coords_path_model - [zeros(n_path,1), period*ones(n_path,1), zeros(n_path,1)]];
    idx_aug_to_path = [(1:n_path)'; (1:n_path)'; (1:n_path)'];

    [idx_aug, dist_nearest] = knnsearch(coords_aug, coords_free_model);
    idx_nearest = idx_aug_to_path(idx_aug);
else
    [idx_nearest, dist_nearest] = knnsearch(coords_path_model, coords_free_model);
end

end

function neighbor_candidates_all = freeRangeSearchCandidates( ...
    coords_model, search_radius, coordinate_mode, reference_radius, full_circumference)

if strcmpi(coordinate_mode, "cylindrical") && full_circumference
    period = 2*pi*reference_radius;
    n_points = size(coords_model,1);
    theta_shift = [zeros(n_points,1), period*ones(n_points,1), zeros(n_points,1)];

    coords_aug = [coords_model; ...
                  coords_model + theta_shift; ...
                  coords_model - theta_shift];
    idx_aug_to_point = [(1:n_points)'; (1:n_points)'; (1:n_points)'];

    [idx_aug_all, ~] = rangesearch(coords_aug, coords_model, search_radius);
    neighbor_candidates_all = cell(size(idx_aug_all));

    for q = 1:numel(idx_aug_all)
        idx_mapped = idx_aug_to_point(idx_aug_all{q});
        neighbor_candidates_all{q} = unique(idx_mapped(:), 'stable');
    end
else
    [neighbor_candidates_all, ~] = rangesearch( ...
        coords_model, coords_model, search_radius);

    for q = 1:numel(neighbor_candidates_all)
        neighbor_candidates_all{q} = neighbor_candidates_all{q}(:);
    end
end

end

function [dxy, dz_signed, center_xy] = freeCylindricalLocalDistances( ...
    P_previous, P_current, cylinder_axis, shaft_center)

rel = bsxfun(@minus, P_previous, P_current);

switch upper(cylinder_axis)
    case 'X'
        U = P_current(2) - shaft_center(1);
        V = P_current(3) - shaft_center(2);
        radius = hypot(U, V);
        if radius <= eps
            error('cylindrical local frame is undefined on the shaft axis.');
        end
        e_radial = [0, U/radius, V/radius];
        e_tangent = [0, -V/radius, U/radius];
        d_axis = rel(:,1);

    case 'Y'
        U = P_current(1) - shaft_center(1);
        V = P_current(3) - shaft_center(2);
        radius = hypot(U, V);
        if radius <= eps
            error('cylindrical local frame is undefined on the shaft axis.');
        end
        e_radial = [U/radius, 0, V/radius];
        e_tangent = [-V/radius, 0, U/radius];
        d_axis = rel(:,2);

    case 'Z'
        U = P_current(1) - shaft_center(1);
        V = P_current(2) - shaft_center(2);
        radius = hypot(U, V);
        if radius <= eps
            error('cylindrical local frame is undefined on the shaft axis.');
        end
        e_radial = [U/radius, V/radius, 0];
        e_tangent = [-V/radius, U/radius, 0];
        d_axis = rel(:,3);

    otherwise
        error('cylinder_axis must be ''X'', ''Y'', or ''Z''.');
end

dz_signed = rel * e_radial(:);
d_tangent = rel * e_tangent(:);
center_xy = [d_axis, d_tangent];
dxy = hypot(d_axis, d_tangent);

end

function C_filtered_current = freeApplyTrackDraggingAtPoint(C_standard_track, P_track, L_melt)

k = size(C_standard_track,1);

if k == 1
    C_filtered_current = C_standard_track(1,:);
    return
end

segment_lengths = zeros(k,1);

for q = 2:k
    segment_lengths(q) = norm(P_track(q,:) - P_track(q-1,:));
end

cumulative_distance = cumsum(segment_lengths);
j_oldest = k;

while j_oldest > 1
    tail_distance = cumulative_distance(k) - cumulative_distance(j_oldest-1);
    tail_weight = exp(-tail_distance / L_melt);
    segment_weight = 1 - exp(-segment_lengths(j_oldest) / L_melt);

    if tail_weight / max(segment_weight, eps) <= 0.1
        break
    end

    j_oldest = j_oldest - 1;
end

C_filtered_current = zeros(1,size(C_standard_track,2));
distance_oldest = cumulative_distance(k) - cumulative_distance(j_oldest);
C_filtered_current = C_filtered_current + ...
    exp(-distance_oldest / L_melt) * C_standard_track(j_oldest,:);

for j = (j_oldest+1):k
    distance_back = cumulative_distance(k) - cumulative_distance(j);
    coef = 1 - exp(-segment_lengths(j) / L_melt);
    weight = coef * exp(-distance_back / L_melt);
    C_filtered_current = C_filtered_current + weight * C_standard_track(j,:);
end

end
