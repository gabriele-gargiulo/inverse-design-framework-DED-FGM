%clear; clc;

%% INVERSE 3D MIXING MODEL
% Finds the required input composition for a desired 3D output.
% The forward operator is replayed without assembling the matrix, using the
% same remelting and in-track dragging logic as the forward model.

tic;

%% INPUT FILES

material_file = 'input_material_distribution_3d_blisk.csv'; % material composition points
toolpath_file = 'toolpath_3d_blisk.csv';                    % deposition toolpath

material_data = readtable(material_file);
toolpath = readtable(toolpath_file);

%% PROCESS PARAMETERS

h = 0.62;          % height
w = 1.97;          % width
b = 0.21;          % penetration depth

% Point and hatch spacing are read directly from the toolpath coordinates.

d_spot = 1.6;      % laser spot diameter
L_melt = 1.5*d_spot; % melt pool length

%% USER OPTIONS

material_names = {'IN','MK','INC'}; % material labels
enable_in_track_grading = true;     % enable correction using in-track grading

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

solver_display = 'iter';        % fmincon display mode
solver_algorithm = 'interior-point'; % fmincon algorithm
solver_optimality_tolerance = 1e-10; % optimality tolerance
solver_constraint_tolerance = 1e-10; % constraint tolerance
solver_step_tolerance = 1e-12;       % step tolerance
solver_max_iterations = 1000;        % maximum fmincon iterations

%% EXTRACT DATA

X_mat = material_data.X;
Y_mat = material_data.Y;
Z_mat = material_data.Z;

IN_des = material_data.IN_frac;
MK_des = material_data.MK_frac;

if ismember('INC_frac', material_data.Properties.VariableNames)
    INC_des = material_data.INC_frac;
else
    INC_des = zeros(height(material_data),1);
end

X_path = toolpath.X;
Y_path = toolpath.Y;
Z_path = toolpath.Z;
path_points_real = [X_path, Y_path, Z_path];
material_points_real = [X_mat, Y_mat, Z_mat];

%% DEPOSITION POINT CLOUD

% The inverse unknowns live on the fixed toolpath/deposition points.
% The desired material distribution is sampled onto these points below.
free_points = path_points_real;
N = size(free_points,1);

if strcmpi(coordinate_mode, "cylindrical")

    [A_mat, Theta_mat, R_mat] = inverse3dFreeCartesianToCylindrical( ...
        X_mat, Y_mat, Z_mat, cylinder_axis, shaft_center, theta_zero);
    [A_path, Theta_path, R_path] = inverse3dFreeCartesianToCylindrical( ...
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
idx_target_material = zeros(N,1);

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
        error('No desired material points found on toolpath layer %d.', k);
    end

    if strcmpi(coordinate_mode, "cylindrical")
        [idx_local, ~] = inverse3dFreeNearestToolpathCylindrical( ...
            material_points_model(idx_material_layer,:), ...
            free_points_model(idx_path_layer,:), ...
            R_ref, full_circumference);
    else
        [idx_local, ~] = knnsearch( ...
            material_points_model(idx_material_layer,:), ...
            free_points_model(idx_path_layer,:));
    end

    idx_target_material(idx_path_layer) = idx_material_layer(idx_local);
end

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
track_id_print = track_id(sort_idx);
layer_print = layer_of_free_point(sort_idx);

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

%% PRECOMPUTE POINT MIXING AND DRAGGING OPERATOR

fprintf('Precomputing free-point mixing operator...\n');

ellipsoid_width = w;
ellipsoid_height = h + b;
neighbor_search_radius = sqrt(ellipsoid_width^2 + ellipsoid_height^2);

if exist('rangesearch','file') == 2
    neighbor_candidates_all = inverse3dFreeRangeSearchCandidates( ...
        free_points_print, neighbor_search_radius, ...
        coordinate_mode, R_ref, full_circumference);
else
    neighbor_candidates_all = {};
end

mix_prev_idx = cell(N,1);
mix_prev_coef = cell(N,1);
mix_u_coef = ones(N,1);

for i = 1:N

    point_curr = free_points_print(i,:);

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

    if ~isempty(prev_idx)

        if strcmpi(coordinate_mode, "cylindrical")
            [dxy, dz_signed, center_xy] = inverse3dFreeCylindricalLocalDistances( ...
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
            inverse3dFreeBeadVerticalRangesOverlap(dz_signed, h, b, coord_tol);

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

            weights = inverse3dFreeExclusiveOverlapFractions( ...
                candidate_centers, candidate_dz, w, h, b, ...
                overlap_partition_samples);
            positive = weights > 0;

            candidates = candidates(positive);
            weights = weights(positive);

            if ~isempty(candidates)
                mix_prev_idx{i} = candidates(:);
                mix_prev_coef{i} = weights(:);
                mix_u_coef(i) = max(0, 1-sum(weights));
            end
        end
    end
end

drag_idx = cell(N,1);
drag_coef = cell(N,1);

if enable_in_track_grading

    for q = 1:numel(track_indices)
        idx_track = track_indices{q};

        for kk = 1:numel(idx_track)
            i = idx_track(kk);
            idx_track_so_far = idx_track(1:kk);
            [drag_idx{i}, drag_coef{i}] = inverse3dFreeDragWeights( ...
                idx_track_so_far, free_points_real_print(idx_track_so_far,:), L_melt);
        end
    end

else

    for i = 1:N
        drag_idx{i} = i;
        drag_coef{i} = 1;
    end
end

operator_data = struct( ...
    'N', N, ...
    'mix_prev_idx', {mix_prev_idx}, ...
    'mix_prev_coef', {mix_prev_coef}, ...
    'mix_u_coef', mix_u_coef, ...
    'drag_idx', {drag_idx}, ...
    'drag_coef', {drag_coef});

fprintf('Free-point inverse operator ready.\n');

%% BUILD TARGET VECTOR y

% Sample the desired material field onto the toolpath/deposition points.
C_material_desired = [IN_des, MK_des, INC_des];
C_desired = C_material_desired(idx_target_material,:);
C_desired_print = C_desired(sort_idx,:);

y_IN = C_desired_print(:,1);
y_MK = C_desired_print(:,2);
y_INC = C_desired_print(:,3);
Y_desired_cols = [y_IN, y_MK, y_INC];

zero_material = [ ...
    all(abs(y_IN) < 1e-12), ...
    all(abs(y_MK) < 1e-12), ...
    all(abs(y_INC) < 1e-12)];

zero_material_idx = find(zero_material);
active_material_idx = find(~zero_material);
two_material_input = numel(active_material_idx) == 2;

if two_material_input
    fprintf('Two-material input detected on sampled toolpath target: %s fixed to zero.\n', ...
        material_names{zero_material_idx});
elseif numel(zero_material_idx) > 1
    fprintf('Multiple unused materials detected on sampled toolpath target: fixing zero target materials to zero.\n');
end

%% SOLVE INVERSE PROBLEM

fprintf('Solving inverse problem...\n');

if isscalar(active_material_idx)

    C_required_print = zeros(N,3);
    C_required_print(:,active_material_idx) = 1;
    exitflag = 1;
    output.message = 'Single-material target detected; no optimization needed.';
    lambda = struct();

elseif two_material_input

    mat_a = active_material_idx(1);
    mat_b = active_material_idx(2);

    y_a = Y_desired_cols(:,mat_a);
    y_b = Y_desired_cols(:,mat_b);

    lb_qp = zeros(N,1);
    ub_qp = ones(N,1);
    z0 = min(max(y_a,0),1);

    opts = optimoptions('fmincon', ...
        'Display', solver_display, ...
        'Algorithm', solver_algorithm, ...
        'SpecifyObjectiveGradient',true, ...
        'HessianApproximation','lbfgs', ...
        'OptimalityTolerance', solver_optimality_tolerance, ...
        'ConstraintTolerance', solver_constraint_tolerance, ...
        'StepTolerance', solver_step_tolerance, ...
        'MaxIterations', solver_max_iterations);

    objective = @(z) inverse3dFreeObjectiveTwoMaterial( ...
        z, y_a, y_b, operator_data);

    [z,~,exitflag,output,lambda] = fmincon( ...
        objective, z0, [], [], [], [], lb_qp, ub_qp, [], opts);

    C_required_print = zeros(N,3);
    C_required_print(:,mat_a) = z;
    C_required_print(:,mat_b) = 1 - z;

else

    Aineq = [speye(N), speye(N)];
    bineq = ones(N,1);
    lb_fmin = zeros(2*N,1);
    ub_fmin = ones(2*N,1);

    z0 = [min(max(y_IN,0),1); min(max(y_MK,0),1)];
    sum_z0 = z0(1:N) + z0(N+1:2*N);
    idx_over_simplex = find(sum_z0 > 1);
    z0(idx_over_simplex) = z0(idx_over_simplex) ./ sum_z0(idx_over_simplex);
    z0(N + idx_over_simplex) = z0(N + idx_over_simplex) ./ sum_z0(idx_over_simplex);

    opts = optimoptions('fmincon', ...
        'Display', solver_display, ...
        'Algorithm', solver_algorithm, ...
        'SpecifyObjectiveGradient',true, ...
        'HessianApproximation','lbfgs', ...
        'OptimalityTolerance', solver_optimality_tolerance, ...
        'ConstraintTolerance', solver_constraint_tolerance, ...
        'StepTolerance', solver_step_tolerance, ...
        'MaxIterations', solver_max_iterations);

    objective = @(z) inverse3dFreeObjectiveThreeMaterial( ...
        z, y_IN, y_MK, y_INC, operator_data);

    [z,~,exitflag,output,lambda] = fmincon( ...
        objective, z0, Aineq, bineq, [], [], lb_fmin, ub_fmin, [], opts);

    IN_req = z(1:N);
    MK_req = z(N+1:2*N);
    INC_req = 1 - IN_req - MK_req;
    C_required_print = [IN_req MK_req INC_req];
end

residual = [ ...
    inverse3dFreeApplyA(C_required_print(:,1), operator_data) - y_IN; ...
    inverse3dFreeApplyA(C_required_print(:,2), operator_data) - y_MK; ...
    inverse3dFreeApplyA(C_required_print(:,3), operator_data) - y_INC];

resnorm = sum(residual.^2);

fprintf('Inverse problem solved.\n');
fprintf('exitflag = %d\n', exitflag);
fprintf('resnorm  = %.12e\n\n', resnorm);

%% RECOVER REQUIRED INPUT DISTRIBUTION

C_required = zeros(size(C_required_print));
C_required(sort_idx,:) = C_required_print;

%% SAVE REQUIRED INPUT

X_free = free_points(:,1);
Y_free = free_points(:,2);
Z_free = free_points(:,3);

required_input = table( ...
    X_free, ...
    Y_free, ...
    Z_free, ...
    C_required(:,1), ...
    C_required(:,2), ...
    C_required(:,3), ...
    'VariableNames', ...
    {'X','Y','Z','IN_frac','MK_frac','INC_frac'});

writetable(required_input, ...
    'required_input_material_distribution_3d_free.csv');

fprintf('Required free-point input distribution saved.\n');

elapsed_time = toc;
fprintf('\nTotal runtime = %.3f seconds\n', elapsed_time);

%% VALIDATION

fprintf('Running validation...\n');

C_out_print = [ ...
    inverse3dFreeApplyA(C_required_print(:,1), operator_data), ...
    inverse3dFreeApplyA(C_required_print(:,2), operator_data), ...
    inverse3dFreeApplyA(C_required_print(:,3), operator_data)];

C_out = zeros(size(C_out_print));
C_out(sort_idx,:) = C_out_print;

C_target_original = zeros(size(C_desired_print));
C_target_original(sort_idx,:) = C_desired_print;

err = C_out - C_target_original;
rmse = sqrt(mean(err(:).^2));

fprintf('\n');
fprintf('====================================\n');
fprintf('RMSE = %.8f\n',rmse);
fprintf('====================================\n');

%% SAVE VALIDATED OUTPUT

validated_output = table( ...
    X_free, ...
    Y_free, ...
    Z_free, ...
    C_out(:,1), ...
    C_out(:,2), ...
    C_out(:,3), ...
    'VariableNames', ...
    {'X','Y','Z','IN_frac','MK_frac','INC_frac'});

writetable(validated_output, ...
    'validated_output_material_distribution_3d_free.csv');

fprintf('Validation output saved.\n');

%% LOCAL FUNCTIONS

function [f,g] = inverse3dFreeObjectiveTwoMaterial(z,y_a,y_b,op)

r_a = inverse3dFreeApplyA(z,op) - y_a;
r_b = inverse3dFreeApplyA(1 - z,op) - y_b;

f = sum(r_a.^2) + sum(r_b.^2);

if nargout > 1
    g = 2*inverse3dFreeApplyAT(r_a,op) - 2*inverse3dFreeApplyAT(r_b,op);
end

end

function [f,g] = inverse3dFreeObjectiveThreeMaterial(z,y_IN,y_MK,y_INC,op)

N = op.N;

IN = z(1:N);
MK = z(N+1:2*N);
INC = 1 - IN - MK;

r_IN = inverse3dFreeApplyA(IN,op) - y_IN;
r_MK = inverse3dFreeApplyA(MK,op) - y_MK;
r_INC = inverse3dFreeApplyA(INC,op) - y_INC;

f = sum(r_IN.^2) + sum(r_MK.^2) + sum(r_INC.^2);

if nargout > 1
    g_IN = 2*inverse3dFreeApplyAT(r_IN,op) - 2*inverse3dFreeApplyAT(r_INC,op);
    g_MK = 2*inverse3dFreeApplyAT(r_MK,op) - 2*inverse3dFreeApplyAT(r_INC,op);
    g = [g_IN; g_MK];
end

end

function y = inverse3dFreeApplyA(u,op)

N = op.N;
y = zeros(N,1);
y_standard = zeros(N,1);

for i = 1:N

    y_standard(i) = op.mix_u_coef(i) * u(i);

    idx_prev = op.mix_prev_idx{i};

    if ~isempty(idx_prev)
        idx_prev = idx_prev(:);
        y_standard(i) = y_standard(i) + ...
            sum(op.mix_prev_coef{i}(:) .* y(idx_prev));
    end

    idx_drag = op.drag_idx{i};
    idx_drag = idx_drag(:);
    y(i) = sum(op.drag_coef{i}(:) .* y_standard(idx_drag));
end

end

function g = inverse3dFreeApplyAT(r,op)

N = op.N;
g = zeros(N,1);
adj_y = r;
adj_standard = zeros(N,1);

for i = N:-1:1

    adj = adj_y(i);

    if adj ~= 0
        idx_drag = op.drag_idx{i};
        idx_drag = idx_drag(:);
        adj_standard(idx_drag) = adj_standard(idx_drag) + ...
            op.drag_coef{i}(:) * adj;
        adj_y(i) = 0;
    end

    adj = adj_standard(i);

    if adj ~= 0
        g(i) = g(i) + op.mix_u_coef(i) * adj;

        idx_prev = op.mix_prev_idx{i};

        if ~isempty(idx_prev)
            idx_prev = idx_prev(:);
            adj_y(idx_prev) = adj_y(idx_prev) + op.mix_prev_coef{i}(:) * adj;
        end

        adj_standard(i) = 0;
    end
end

end

function [A, Theta, R] = inverse3dFreeCartesianToCylindrical(X, Y, Z, cylinder_axis, shaft_center, theta_zero)

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

function [idx_nearest, dist_nearest] = inverse3dFreeNearestToolpathCylindrical( ...
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

function neighbor_candidates_all = inverse3dFreeRangeSearchCandidates( ...
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

function [dxy, dz_signed, center_xy] = inverse3dFreeCylindricalLocalDistances( ...
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

function [idx_global, weights] = inverse3dFreeDragWeights(idx_track_so_far, P_track, L_melt)

k = numel(idx_track_so_far);

if k == 1
    idx_global = idx_track_so_far(1);
    weights = 1;
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

idx_global = idx_track_so_far(j_oldest:k);
idx_global = idx_global(:);
weights = zeros(numel(idx_global),1);

distance_oldest = cumulative_distance(k) - cumulative_distance(j_oldest);
weights(1) = exp(-distance_oldest / L_melt);

for j = (j_oldest+1):k
    local_pos = j - j_oldest + 1;
    distance_back = cumulative_distance(k) - cumulative_distance(j);
    coef = 1 - exp(-segment_lengths(j) / L_melt);
    weights(local_pos) = coef * exp(-distance_back / L_melt);
end

end

function fractions = inverse3dFreeExclusiveOverlapFractions( ...
    center_xy, dz_signed, w, h, b, n_samples)

center_xy = reshape(center_xy, [], 2);
dz_signed = dz_signed(:);
n_candidates = size(center_xy,1);
fractions = zeros(n_candidates,1);

if n_candidates == 0
    return
end

unit_points = inverse3dFreeUnitBallSamples(n_samples);
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


function points = inverse3dFreeUnitBallSamples(n_samples)

persistent cached_count cached_points

n_samples = max(1024, round(n_samples));
if ~isempty(cached_count) && cached_count == n_samples
    points = cached_points;
    return
end

sample_index = (1:n_samples)';
u_radius = inverse3dFreeRadicalInverse(sample_index, 2);
u_polar = inverse3dFreeRadicalInverse(sample_index, 3);
u_azimuth = inverse3dFreeRadicalInverse(sample_index, 5);

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

function values = inverse3dFreeRadicalInverse(sample_index, base)

values = zeros(size(sample_index));
factor = 1/base;

while any(sample_index > 0)
    values = values + factor*mod(sample_index,base);
    sample_index = floor(sample_index/base);
    factor = factor/base;
end

end

function has_overlap = inverse3dFreeBeadVerticalRangesOverlap(dz_signed, h, b, coord_tol)

previous_bottom = dz_signed - b;
previous_top = dz_signed + h;
current_bottom = -b;
current_top = h;

z0 = max(current_bottom, previous_bottom);
z1 = min(current_top, previous_top);
has_overlap = z1 >= z0 - coord_tol;

end
