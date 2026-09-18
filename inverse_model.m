%clear; clc;

%% INVERSE 3D MIXING MODEL
% Finds the required input composition for a desired 3D output.
% The forward operator is replayed without assembling the matrix, using the
% same remelting and in-track dragging logic as the forward model.
%
% Model sequence:
%   1. Sample the target composition onto the ordered deposition points.
%   2. Recover layers/tracks and precompute the causal forward operator.
%   3. Minimize the summed squared composition residual subject to physical
%      material-fraction bounds and the unit-sum constraint.
%   4. Reapply the forward operator to the optimized input for validation.
%   5. Save both the required input and the corresponding validated output.
%
% The forward map is linear in each material fraction, but it is applied as a
% sequence of remelting and dragging operations rather than stored as a large
% dense matrix. Its transpose is evaluated by an explicit reverse/adjoint pass.
% Coordinates and geometric process parameters are expressed in millimetres.

% Start a wall-clock timer for precomputation, optimization, and validation.
tic;

%% INPUT FILES

material_file = 'input_material_distribution_3d_blisk.csv'; % desired output-composition field
toolpath_file = 'toolpath_3d_blisk.csv';                    % ordered deposition coordinates

% Import both CSV inputs as tables so columns can be addressed by name.
material_data = readtable(material_file);
toolpath = readtable(toolpath_file);

%% PROCESS PARAMETERS

h = 0.62;          % deposited bead height above its centre plane [mm]
w = 1.97;          % full lateral bead/melt-pool width [mm]
b = 0.21;          % penetration depth below the bead centre plane [mm]

% Point spacing and inter-track travel distances are read directly from the
% ordered toolpath coordinates.

d_spot = 1.6;        % laser spot diameter [mm]
L_melt = 1.5*d_spot; % characteristic melt-pool dragging length [mm]

%% USER OPTIONS

material_names = {'IN','MK','INC'}; % labels corresponding to composition columns
enable_in_track_grading = true;     % include in-track dragging in inverse map

% Number of deterministic quasi-Monte-Carlo samples used to integrate each
% asymmetric ellipsoidal overlap. Higher values increase accuracy and cost.
overlap_partition_samples = 16384;

% "cartesian" treats Z as the build direction.
% "cylindrical" unwraps the part and treats radius as the build direction.
coordinate_mode = "cartesian";
cylinder_axis = 'Y';        % shaft axis used only in cylindrical mode
shaft_center = [0, 0];      % shaft centre in the plane normal to the axis [mm]
theta_zero = 0;             % angular location of the unwrapped seam [rad]
reference_radius_mode = "inner"; % "inner" radius or mean toolpath radius
full_circumference = true;  % make angular neighbour searches periodic
track_break_distance = 2*w; % consecutive-point jump that starts a new track

% fmincon settings governing console output, convergence, and work limit.
solver_display = 'iter';        % fmincon display mode
solver_algorithm = 'interior-point'; % fmincon algorithm
solver_optimality_tolerance = 1e-10; % optimality tolerance
solver_constraint_tolerance = 1e-10; % constraint tolerance
solver_step_tolerance = 1e-12;       % step tolerance
solver_max_iterations = 1000;        % maximum fmincon iterations

%% EXTRACT DATA

% Coordinates at which the desired output composition is prescribed.
X_mat = material_data.X;
Y_mat = material_data.Y;
Z_mat = material_data.Z;

% The first two material fractions are required by the input-file format.
IN_des = material_data.IN_frac;
MK_des = material_data.MK_frac;

% Permit legacy two-material targets by creating a zero INC fraction.
if ismember('INC_frac', material_data.Properties.VariableNames)
    INC_des = material_data.INC_frac;
% Missing third-material data represent a two-material target field.
else
    INC_des = zeros(height(material_data),1);
end

% Toolpath rows are already ordered chronologically by deposition time.
X_path = toolpath.X;
Y_path = toolpath.Y;
Z_path = toolpath.Z;
path_points_real = [X_path, Y_path, Z_path];
material_points_real = [X_mat, Y_mat, Z_mat];

%% DEPOSITION POINT CLOUD

% The inverse unknowns live on the fixed toolpath/deposition points.
% The desired material distribution is sampled onto these points below.
free_points = path_points_real;
N = size(free_points,1); % number of unknown feed vectors/deposition points

% Build an unwrapped metric coordinate system when modelling a cylindrical
% part. The model coordinates are axial distance, arc length, and radius.
if strcmpi(coordinate_mode, "cylindrical")

    [A_mat, Theta_mat, R_mat] = inverse3dFreeCartesianToCylindrical( ...
        X_mat, Y_mat, Z_mat, cylinder_axis, shaft_center, theta_zero);
    [A_path, Theta_path, R_path] = inverse3dFreeCartesianToCylindrical( ...
        X_path, Y_path, Z_path, cylinder_axis, shaft_center, theta_zero);

    % Select the radius that converts angular separation into arc length.
    if strcmpi(reference_radius_mode, "inner")
        R_ref = min(R_path);
    % The alternative uses one representative radius over the whole path.
    else
        R_ref = mean(R_path);
    end

    % A zero/negative radius would make the angular metric undefined.
    if R_ref <= 0
        error('reference radius must be positive. Check shaft_center and cylindrical coordinates.');
    end

    % Replace theta by R_ref*theta so all model coordinates have mm units.
    material_points_model = [A_mat, R_ref*Theta_mat, R_mat];
    coords_path = [A_path, R_ref*Theta_path, R_path];
    free_points_model = coords_path;

% Retain XYZ geometry when no cylindrical unwrapping was requested.
else

    % Cartesian data require no coordinate transformation.
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

%% TOOLPATH LAYERS AND TRACKS

% Infer layer and track membership directly from the ordered coordinates.
% Track detection is independent of the global X and Y axes, allowing
% straight, diagonal, and curved paths. The forward model uses the same rule.

coord_tol = 1e-9;

build_coord_path = Z_path; % Z in Cartesian mode, radius in cylindrical mode
layer_of_path_point = zeros(length(Z_path),1);

% Treat an empty toolpath explicitly; otherwise start the first layer at 1.
if isempty(build_coord_path)
    n_layers_path = 0;
% Nonempty paths always contain at least the initial layer.
else
    current_layer = 1;
    layer_of_path_point(1) = current_layer;

    % A positive build-coordinate jump larger than half a bead height marks
    % the first point of a newly deposited layer.
    for q = 2:length(build_coord_path)
        % Increment the layer counter only at an upward build jump.
        if build_coord_path(q) - build_coord_path(q-1) > h/2
            current_layer = current_layer + 1;
        end

        layer_of_path_point(q) = current_layer;
    end

    n_layers_path = current_layer;
end

% Sample the target field by nearest neighbour within each detected layer;
% layer restriction prevents geometrically close adjacent layers being mixed.
material_build_coord = material_points_model(:,3);
idx_target_material = zeros(N,1);

% Process each layer independently during target-field sampling.
for k = 1:n_layers_path
    idx_path_layer = find(layer_of_path_point == k);
    build_coord_layer = build_coord_path(idx_path_layer);
    layer_min = min(build_coord_layer) - coord_tol;
    layer_max = max(build_coord_layer) + coord_tol;
    idx_material_layer = find( ...
        material_build_coord >= layer_min & ...
        material_build_coord <= layer_max);

    % If coordinate noise prevents an exact range match, accept material
    % points within half a nominal bead height of the layer centre.
    if isempty(idx_material_layer)
        layer_center = mean(build_coord_layer);
        idx_material_layer = find(abs(material_build_coord - layer_center) <= h/2 + coord_tol);
    end

    % Stop with a specific layer number rather than silently using bad data.
    if isempty(idx_material_layer)
        error('No desired material points found on toolpath layer %d.', k);
    end

    % Use a seam-periodic nearest-neighbour lookup for a full cylinder.
    if strcmpi(coordinate_mode, "cylindrical")
        [idx_local, ~] = inverse3dFreeNearestToolpathCylindrical( ...
            material_points_model(idx_material_layer,:), ...
            free_points_model(idx_path_layer,:), ...
            R_ref, full_circumference);
    % Cartesian matching uses ordinary Euclidean nearest neighbours.
    else
        [idx_local, ~] = knnsearch( ...
            material_points_model(idx_material_layer,:), ...
            free_points_model(idx_path_layer,:));
    end

    idx_target_material(idx_path_layer) = idx_material_layer(idx_local);
end

% Segment the ordered toolpath into physical tracks. Every layer begins a new
% track. Within a layer, consecutive points remain on the same track unless
% their real-space separation exceeds the prescribed travel-jump threshold.
% No assumption is made about track orientation or curvature.
path_track_id = zeros(length(Z_path),1);
next_track_id = 0;

% Walk through every detected layer in chronological order.
for k = 1:n_layers_path

    idx_layer = find(layer_of_path_point == k);

    % Defensive guard for a layer label with no associated path rows.
    if isempty(idx_layer)
        continue
    end

    % The first point in every layer necessarily starts a new track.
    next_track_id = next_track_id + 1;
    path_track_id(idx_layer(1)) = next_track_id;

    % Inspect every transition between consecutive points in this layer.
    for q = 2:numel(idx_layer)
        idx_prev = idx_layer(q-1);
        idx_curr = idx_layer(q);

        toolpath_jump = norm( ...
            path_points_real(idx_curr,:) - path_points_real(idx_prev,:));

        % A sufficiently large travel jump starts a new physical track.
        if toolpath_jump > track_break_distance
            next_track_id = next_track_id + 1;
        end

        path_track_id(idx_curr) = next_track_id;
    end
end

%% ASSIGN POINTS TO TOOLPATH LAYERS AND TRACKS

% Transfer path-derived layer and track IDs to the inverse unknowns.
layer_of_free_point = layer_of_path_point(idx_toolpath);

track_id = path_track_id(idx_toolpath);

% Normally every point is assigned; this fallback avoids invalid zero indices
% if a malformed toolpath happens to escape the segmentation logic.
unassigned_track = track_id == 0;
% Apply the defensive fallback only when at least one zero ID is present.
if any(unassigned_track)
    track_id(unassigned_track) = 1;
end

%% BUILD TRUE PRINT ORDER FROM TOOLPATH

% Sort by toolpath row and then distance along that row. For the present
% one-point-per-row representation, this retains the original CSV order.
sort_matrix = [idx_toolpath(:), dist_tool(:)];
[~, sort_idx] = sortrows(sort_matrix, [1 2]);

free_points_print = free_points_model(sort_idx,:);
free_points_real_print = free_points(sort_idx,:);
track_id_print = track_id(sort_idx);
layer_print = layer_of_free_point(sort_idx);

% Store the inverse permutation for diagnostics and order restoration.
inverse_sort_idx = zeros(N,1);
inverse_sort_idx(sort_idx) = 1:N;

track_ids_unique = unique(track_id_print, 'stable');
track_indices = cell(numel(track_ids_unique),1);
track_group_print = zeros(N,1);
position_in_track = zeros(N,1);
track_distance_print = zeros(N,1);

% Cache membership, local index, and cumulative physical distance for every
% track; the dragging kernel depends on distance rather than sample count.
for q = 1:numel(track_ids_unique)
    idx_track = find(track_id_print == track_ids_unique(q));
    track_indices{q} = idx_track;
    track_group_print(idx_track) = q;
    position_in_track(idx_track) = 1:numel(idx_track);

    P_track = free_points_real_print(idx_track,:);
    % A single-point track has zero cumulative distance by construction.
    if numel(idx_track) > 1
        segment_lengths = sqrt(sum(diff(P_track,1,1).^2,2));
        track_distance_print(idx_track) = [0; cumsum(segment_lengths)];
    end
end

%% PRECOMPUTE POINT MIXING AND DRAGGING OPERATOR

fprintf('Precomputing free-point mixing operator...\n');

ellipsoid_width = w;      % total transverse footprint used for broad search
ellipsoid_height = h + b; % full top-to-bottom bead extent
neighbor_search_radius = sqrt(ellipsoid_width^2 + ellipsoid_height^2);

% Use a range search when available; an empty-cell fallback below evaluates
% all previous points and therefore changes performance, not model physics.
if exist('rangesearch','file') == 2
    neighbor_candidates_all = inverse3dFreeRangeSearchCandidates( ...
        free_points_print, neighbor_search_radius, ...
        coordinate_mode, R_ref, full_circumference);
% Without rangesearch, an empty marker triggers exhaustive prior-point scans.
else
    neighbor_candidates_all = {};
end

% For point i, mix_prev_idx/mix_prev_coef store already realized output
% compositions that are remelted. mix_u_coef stores the fresh-feed fraction.
mix_prev_idx = cell(N,1);
mix_prev_coef = cell(N,1);
mix_u_coef = ones(N,1);

% March forward through print time to encode the causal remelting dependency.
for i = 1:N

    point_curr = free_points_print(i,:);

    % Select either every prior point or only precomputed nearby candidates.
    if isempty(neighbor_candidates_all)
        prev_idx = 1:i-1;
    % Range-search candidates still must be restricted to earlier print times.
    else
        prev_idx = neighbor_candidates_all{i};
        prev_idx = prev_idx(prev_idx < i);
    end

    % Dense points on the current continuous track represent the same
    % moving meltpool and are handled by dragging, not by remelting.
    % Exclude nearby points in the same moving melt pool. They are handled by
    % the dragging operator and must not also contribute as remelted beads.
    if ~isempty(prev_idx)
        same_track = track_id_print(prev_idx) == track_id_print(i);
        distance_behind = track_distance_print(i) - track_distance_print(prev_idx);
        prev_idx(same_track & distance_behind < w) = [];
    end

    % Continue with exact overlap tests only when prior candidates remain.
    if ~isempty(prev_idx)

        % Express offsets in a local bead frame. Cylindrical mode resolves
        % axial/tangential/radial components at the current deposition point.
        if strcmpi(coordinate_mode, "cylindrical")
            [dxy, dz_signed, center_xy] = inverse3dFreeCylindricalLocalDistances( ...
                free_points_real_print(prev_idx,:), ...
                free_points_real_print(i,:), ...
                cylinder_axis, shaft_center);
        % In Cartesian mode, XYZ differences already form the local bead frame.
        else
            delta = free_points_print(prev_idx,:) - point_curr;
            center_xy = delta(:,1:2);
            dxy = hypot(delta(:,1), delta(:,2));
            dz_signed = delta(:,3);
        end

        % Broad geometric eligibility: transverse proximity plus intersection
        % of the asymmetric build-direction bead intervals [-b,h].
        meltpool_intersects = ...
            dxy <= ellipsoid_width + coord_tol & ...
            inverse3dFreeBeadVerticalRangesOverlap(dz_signed, h, b, coord_tol);

        same_layer = layer_print(prev_idx) == layer_print(i);
        lower_layer = layer_print(prev_idx) < layer_print(i);
        eligible = meltpool_intersects & (same_layer | lower_layer);

        % Partition the current bead only when at least one candidate passes.
        if any(eligible)
            candidates = prev_idx(eligible);
            candidate_centers = center_xy(eligible,:);
            candidate_dz = dz_signed(eligible);

            % Chronological order is required by the latest-bead-wins rule.
            % Chronological sorting allows reverse volume assignment inside
            % the partition routine, implementing latest-bead-wins ownership.
            [candidates, order] = sort(candidates);
            candidate_centers = candidate_centers(order,:);
            candidate_dz = candidate_dz(order);

            % Each weight is one candidate's exclusive fraction of current
            % bead volume after resolving overlaps among prior candidates.
            weights = inverse3dFreeExclusiveOverlapFractions( ...
                candidate_centers, candidate_dz, w, h, b, ...
                overlap_partition_samples);
            positive = weights > 0;

            candidates = candidates(positive);
            weights = weights(positive);

            % Cache sparse dependencies; the remaining fraction is new feed.
            if ~isempty(candidates)
                mix_prev_idx{i} = candidates(:);
                mix_prev_coef{i} = weights(:);
                mix_u_coef(i) = max(0, 1-sum(weights));
            end
        end
    end
end

% drag_idx/drag_coef define the linear exponential filter applied to each
% point's sequence of post-remelting compositions within its current track.
drag_idx = cell(N,1);
drag_coef = cell(N,1);

% Either reproduce full in-track dragging or reduce the filter to identity.
if enable_in_track_grading

    % Precompute a causal dragging stencil for every prefix of every track.
    for q = 1:numel(track_indices)
        idx_track = track_indices{q};

        % Build a stencil for every successive current point on this track.
        for kk = 1:numel(idx_track)
            i = idx_track(kk);
            idx_track_so_far = idx_track(1:kk);
            [drag_idx{i}, drag_coef{i}] = inverse3dFreeDragWeights( ...
                idx_track_so_far, free_points_real_print(idx_track_so_far,:), L_melt);
        end
    end

% Disabling in-track grading replaces every dragging stencil with identity.
else

    % Identity weights allow remelting-only inverse studies with one switch.
    for i = 1:N
        drag_idx{i} = i;
        drag_coef{i} = 1;
    end
end

% Bundle all sparse recurrence data into one object passed to the optimizer.
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
% Gather the target composition that was mapped to each deposition point.
C_material_desired = [IN_des, MK_des, INC_des];
C_desired = C_material_desired(idx_target_material,:);
C_desired_print = C_desired(sort_idx,:);

y_IN = C_desired_print(:,1);
y_MK = C_desired_print(:,2);
y_INC = C_desired_print(:,3);
Y_desired_cols = [y_IN, y_MK, y_INC];

% Detect materials absent everywhere in the sampled target. Eliminating those
% variables reduces problem size and exactly preserves their zero fraction.
zero_material = [ ...
    all(abs(y_IN) < 1e-12), ...
    all(abs(y_MK) < 1e-12), ...
    all(abs(y_INC) < 1e-12)];

zero_material_idx = find(zero_material);
active_material_idx = find(~zero_material);
two_material_input = numel(active_material_idx) == 2;

% Report which lower-dimensional composition simplex will be optimized.
if two_material_input
    fprintf('Two-material input detected on sampled toolpath target: %s fixed to zero.\n', ...
        material_names{zero_material_idx});
% More than one absent material means a single-material target.
elseif numel(zero_material_idx) > 1
    fprintf('Multiple unused materials detected on sampled toolpath target: fixing zero target materials to zero.\n');
end

%% SOLVE INVERSE PROBLEM

fprintf('Solving inverse problem...\n');

% A single active material has the unique feasible input fraction of one, so
% no numerical optimization is necessary.
if isscalar(active_material_idx)

    C_required_print = zeros(N,3);
    C_required_print(:,active_material_idx) = 1;
    exitflag = 1;
    output.message = 'Single-material target detected; no optimization needed.';
    lambda = struct();

% Two active materials reduce the pointwise simplex to one scalar variable.
elseif two_material_input

    % In a two-material system, optimize only material A; material B is 1-z.
    mat_a = active_material_idx(1);
    mat_b = active_material_idx(2);

    y_a = Y_desired_cols(:,mat_a);
    y_b = Y_desired_cols(:,mat_b);

    % Box bounds enforce both fractions in [0,1]. The clipped target is a
    % feasible and physically meaningful starting composition.
    lb_qp = zeros(N,1);
    ub_qp = ones(N,1);
    z0 = min(max(y_a,0),1);

    % Supply the analytical gradient and use a limited-memory Hessian model.
    opts = optimoptions('fmincon', ...
        'Display', solver_display, ...
        'Algorithm', solver_algorithm, ...
        'SpecifyObjectiveGradient',true, ...
        'HessianApproximation','lbfgs', ...
        'OptimalityTolerance', solver_optimality_tolerance, ...
        'ConstraintTolerance', solver_constraint_tolerance, ...
        'StepTolerance', solver_step_tolerance, ...
        'MaxIterations', solver_max_iterations);

    % The objective includes residuals of both complementary materials.
    objective = @(z) inverse3dFreeObjectiveTwoMaterial( ...
        z, y_a, y_b, operator_data);

    [z,~,exitflag,output,lambda] = fmincon( ...
        objective, z0, [], [], [], [], lb_qp, ub_qp, [], opts);

    C_required_print = zeros(N,3);
    C_required_print(:,mat_a) = z;
    C_required_print(:,mat_b) = 1 - z;

% Otherwise solve the full three-material constrained problem.
else

    % For three materials, optimize IN and MK explicitly and recover
    % INC = 1-IN-MK. Aineq enforces the nonnegative INC constraint.
    Aineq = [speye(N), speye(N)];
    bineq = ones(N,1);
    lb_fmin = zeros(2*N,1);
    ub_fmin = ones(2*N,1);

    % Start from the clipped target fractions, then project any row whose IN
    % and MK sum exceeds one back onto the unit-simplex boundary.
    z0 = [min(max(y_IN,0),1); min(max(y_MK,0),1)];
    sum_z0 = z0(1:N) + z0(N+1:2*N);
    idx_over_simplex = find(sum_z0 > 1);
    z0(idx_over_simplex) = z0(idx_over_simplex) ./ sum_z0(idx_over_simplex);
    z0(N + idx_over_simplex) = z0(N + idx_over_simplex) ./ sum_z0(idx_over_simplex);

    % Supply the analytical gradient and use a limited-memory Hessian model.
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

% Reapply the operator independently to all material columns and concatenate
% residuals exactly as they are measured by the least-squares objective.
residual = [ ...
    inverse3dFreeApplyA(C_required_print(:,1), operator_data) - y_IN; ...
    inverse3dFreeApplyA(C_required_print(:,2), operator_data) - y_MK; ...
    inverse3dFreeApplyA(C_required_print(:,3), operator_data) - y_INC];

resnorm = sum(residual.^2);

fprintf('Inverse problem solved.\n');
fprintf('exitflag = %d\n', exitflag);
fprintf('resnorm  = %.12e\n\n', resnorm);

%% RECOVER REQUIRED INPUT DISTRIBUTION

% Undo the print-order permutation before constructing the output table.
C_required = zeros(size(C_required_print));
C_required(sort_idx,:) = C_required_print;

%% SAVE REQUIRED INPUT

% Save original physical coordinates, including for cylindrical calculations.
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

% The fixed filename is used by deviation_plot.m in correction mode.
writetable(required_input, ...
    'required_input_material_distribution_3d_free.csv');

fprintf('Required free-point input distribution saved.\n');

% Stop and report the timer after the required input has been saved.
elapsed_time = toc;
fprintf('\nTotal runtime = %.3f seconds\n', elapsed_time);

%% VALIDATION

fprintf('Running validation...\n');

% Forward-propagate the optimized fractions without rerunning geometry setup.
C_out_print = [ ...
    inverse3dFreeApplyA(C_required_print(:,1), operator_data), ...
    inverse3dFreeApplyA(C_required_print(:,2), operator_data), ...
    inverse3dFreeApplyA(C_required_print(:,3), operator_data)];

% Restore both prediction and target to their original toolpath row order.
C_out = zeros(size(C_out_print));
C_out(sort_idx,:) = C_out_print;

C_target_original = zeros(size(C_desired_print));
C_target_original(sort_idx,:) = C_desired_print;

% RMSE here is over every point-material scalar, not over vector magnitudes.
err = C_out - C_target_original;
rmse = sqrt(mean(err(:).^2));

fprintf('\n');
fprintf('====================================\n');
fprintf('RMSE = %.8f\n',rmse);
fprintf('====================================\n');

%% SAVE VALIDATED OUTPUT

% Store the forward prediction produced by the optimized/corrected input.
validated_output = table( ...
    X_free, ...
    Y_free, ...
    Z_free, ...
    C_out(:,1), ...
    C_out(:,2), ...
    C_out(:,3), ...
    'VariableNames', ...
    {'X','Y','Z','IN_frac','MK_frac','INC_frac'});

% This fixed filename is used by deviation_plot.m for corrected-output error.
writetable(validated_output, ...
    'validated_output_material_distribution_3d_free.csv');

fprintf('Validation output saved.\n');

%% LOCAL FUNCTIONS

function [f,g] = inverse3dFreeObjectiveTwoMaterial(z,y_a,y_b,op)

%INVERSE3DFREEOBJECTIVETWOMATERIAL Two-component least-squares objective.
% z is the input fraction of material A and 1-z is material B. The function
% returns the summed squared forward-model residual and, when requested, its
% exact gradient computed with the transpose/adjoint operator.

% Propagate both complementary input fractions through the forward operator.
r_a = inverse3dFreeApplyA(z,op) - y_a;
r_b = inverse3dFreeApplyA(1 - z,op) - y_b;

% Penalize target mismatch for both materials symmetrically.
f = sum(r_a.^2) + sum(r_b.^2);

% Avoid the adjoint work when fmincon requests only the objective value.
if nargout > 1
    g = 2*inverse3dFreeApplyAT(r_a,op) - 2*inverse3dFreeApplyAT(r_b,op);
end

end

function [f,g] = inverse3dFreeObjectiveThreeMaterial(z,y_IN,y_MK,y_INC,op)

%INVERSE3DFREEOBJECTIVETHREEMATERIAL Three-component least-squares objective.
% Only IN and MK are independent optimization variables; INC is recovered by
% the per-point unit-sum constraint. The inequality constraints supplied to
% fmincon ensure all three reconstructed fractions remain nonnegative.

N = op.N;

% Split the stacked decision vector and recover the dependent third fraction.
IN = z(1:N);
MK = z(N+1:2*N);
INC = 1 - IN - MK;

% Forward-propagate each material fraction and form target residuals.
r_IN = inverse3dFreeApplyA(IN,op) - y_IN;
r_MK = inverse3dFreeApplyA(MK,op) - y_MK;
r_INC = inverse3dFreeApplyA(INC,op) - y_INC;

% The scalar objective is the total SSE over all points and materials.
f = sum(r_IN.^2) + sum(r_MK.^2) + sum(r_INC.^2);

% Chain-rule signs are negative for INC because INC = 1-IN-MK.
if nargout > 1
    g_IN = 2*inverse3dFreeApplyAT(r_IN,op) - 2*inverse3dFreeApplyAT(r_INC,op);
    g_MK = 2*inverse3dFreeApplyAT(r_MK,op) - 2*inverse3dFreeApplyAT(r_INC,op);
    g = [g_IN; g_MK];
end

end

function y = inverse3dFreeApplyA(u,op)

%INVERSE3DFREEAPPLYA Apply the matrix-free forward operator to one material.
% The first recurrence forms the post-remelting value y_standard(i) from new
% feed u(i) and already realized outputs y(previous). The second recurrence
% applies the precomputed in-track dragging stencil to those standard values.

N = op.N;
y = zeros(N,1);
y_standard = zeros(N,1);

% Forward print-order traversal is required by the causal remelting recurrence.
for i = 1:N

    % Fresh material contributes through the unremelted bead-volume fraction.
    y_standard(i) = op.mix_u_coef(i) * u(i);

    idx_prev = op.mix_prev_idx{i};

    % Add already deposited compositions occupying remelted volume.
    if ~isempty(idx_prev)
        idx_prev = idx_prev(:);
        y_standard(i) = y_standard(i) + ...
            sum(op.mix_prev_coef{i}(:) .* y(idx_prev));
    end

    % Exponentially average the retained post-remelting track history.
    idx_drag = op.drag_idx{i};
    idx_drag = idx_drag(:);
    y(i) = sum(op.drag_coef{i}(:) .* y_standard(idx_drag));
end

end

function g = inverse3dFreeApplyAT(r,op)

%INVERSE3DFREEAPPLYAT Apply the exact transpose of the matrix-free operator.
% Reverse accumulation propagates output sensitivities first through dragging
% and then through remelting. This produces A'*r without assembling A and is
% used to give fmincon analytical least-squares gradients.

N = op.N;
g = zeros(N,1);
adj_y = r;
adj_standard = zeros(N,1);

% Reverse print order is the adjoint of the causal forward recurrence.
for i = N:-1:1

    adj = adj_y(i);

    % Scatter output sensitivity through the dragging weights.
    if adj ~= 0
        idx_drag = op.drag_idx{i};
        idx_drag = idx_drag(:);
        adj_standard(idx_drag) = adj_standard(idx_drag) + ...
            op.drag_coef{i}(:) * adj;
        adj_y(i) = 0;
    end

    adj = adj_standard(i);

    % Scatter post-remelting sensitivity to fresh feed and previous outputs.
    if adj ~= 0
        g(i) = g(i) + op.mix_u_coef(i) * adj;

        idx_prev = op.mix_prev_idx{i};

        % Earlier outputs receive sensitivity proportional to overlap volume.
        if ~isempty(idx_prev)
            idx_prev = idx_prev(:);
            adj_y(idx_prev) = adj_y(idx_prev) + op.mix_prev_coef{i}(:) * adj;
        end

        adj_standard(i) = 0;
    end
end

end

function [A, Theta, R] = inverse3dFreeCartesianToCylindrical(X, Y, Z, cylinder_axis, shaft_center, theta_zero)

%INVERSE3DFREECARTESIANTOCYLINDRICAL Convert XYZ to axial/angular/radial form.
% A is the coordinate along the selected shaft axis; U and V span the normal
% plane; Theta is wrapped to [0,2*pi); R is distance from the shaft axis.

% Select the axial coordinate and transverse plane for the requested axis.
switch upper(cylinder_axis)
    % Shaft along X leaves Y-Z as the transverse plane.
    case 'X'
        A = X;
        U = Y - shaft_center(1);
        V = Z - shaft_center(2);
    % Shaft along Y leaves X-Z as the transverse plane.
    case 'Y'
        A = Y;
        U = X - shaft_center(1);
        V = Z - shaft_center(2);
    % Shaft along Z leaves X-Y as the transverse plane.
    case 'Z'
        A = Z;
        U = X - shaft_center(1);
        V = Y - shaft_center(2);
    % Reject misspelled axes before returning meaningless geometry.
    otherwise
        error('cylinder_axis must be ''X'', ''Y'', or ''Z''.');
end

% Rotate the angular origin to theta_zero and wrap across the chosen seam.
Theta = mod(atan2(V, U) - theta_zero, 2*pi);
R = hypot(U, V);

end

function [idx_nearest, dist_nearest] = inverse3dFreeNearestToolpathCylindrical( ...
    coords_path_model, coords_free_model, reference_radius, full_circumference)

%INVERSE3DFREENEARESTTOOLPATHCYLINDRICAL Find periodic nearest neighbours.
% Full cylinders duplicate the reference field one circumference to either
% side so points adjacent across the unwrapped seam remain close.

% Add angular-period images only for a complete circumference.
if full_circumference
    period = 2*pi*reference_radius;
    n_path = size(coords_path_model,1);
    coords_aug = [coords_path_model; ...
                  coords_path_model + [zeros(n_path,1), period*ones(n_path,1), zeros(n_path,1)]; ...
                  coords_path_model - [zeros(n_path,1), period*ones(n_path,1), zeros(n_path,1)]];
    idx_aug_to_path = [(1:n_path)'; (1:n_path)'; (1:n_path)'];

    % Translate augmented-image indices back to original material rows.
    [idx_aug, dist_nearest] = knnsearch(coords_aug, coords_free_model);
    idx_nearest = idx_aug_to_path(idx_aug);
% Partial cylinders have no periodic seam and use direct nearest neighbours.
else
    [idx_nearest, dist_nearest] = knnsearch(coords_path_model, coords_free_model);
end

end

function neighbor_candidates_all = inverse3dFreeRangeSearchCandidates( ...
    coords_model, search_radius, coordinate_mode, reference_radius, full_circumference)

%INVERSE3DFREERANGESEARCHCANDIDATES Precompute nearby deposition-point IDs.
% This is only a performance filter. Exact bead intersection and layer tests
% are still evaluated when the sparse forward recurrence is constructed.

% Periodic images make the unwrapped cylindrical seam transparent to search.
if strcmpi(coordinate_mode, "cylindrical") && full_circumference
    period = 2*pi*reference_radius;
    n_points = size(coords_model,1);
    theta_shift = [zeros(n_points,1), period*ones(n_points,1), zeros(n_points,1)];

    coords_aug = [coords_model; ...
                  coords_model + theta_shift; ...
                  coords_model - theta_shift];
    idx_aug_to_point = [(1:n_points)'; (1:n_points)'; (1:n_points)'];

    % Search the augmented coordinates and map results to original point IDs.
    [idx_aug_all, ~] = rangesearch(coords_aug, coords_model, search_radius);
    neighbor_candidates_all = cell(size(idx_aug_all));

    % Stable uniqueness removes duplicate images without reordering point IDs.
    for q = 1:numel(idx_aug_all)
        idx_mapped = idx_aug_to_point(idx_aug_all{q});
        neighbor_candidates_all{q} = unique(idx_mapped(:), 'stable');
    end
% Use the original coordinate set when periodic images are unnecessary.
else
    % Cartesian and partial-cylinder cases require no periodic duplication.
    [neighbor_candidates_all, ~] = rangesearch( ...
        coords_model, coords_model, search_radius);

    % Store candidate lists as column vectors for consistent later indexing.
    for q = 1:numel(neighbor_candidates_all)
        neighbor_candidates_all{q} = neighbor_candidates_all{q}(:);
    end
end

end

function [dxy, dz_signed, center_xy] = inverse3dFreeCylindricalLocalDistances( ...
    P_previous, P_current, cylinder_axis, shaft_center)

%INVERSE3DFREECYLINDRICALLOCALDISTANCES Resolve offsets in a local frame.
% center_xy contains axial and tangential offsets; dxy is their resultant;
% dz_signed is the radial/build-direction offset. The local frame is evaluated
% at the current point because it rotates around the cylinder.

% Vector from the current deposition point to each previous candidate.
rel = bsxfun(@minus, P_previous, P_current);

% Construct radial and tangential unit vectors for the selected shaft axis.
switch upper(cylinder_axis)
    % Resolve the local basis around an X-oriented shaft.
    case 'X'
        U = P_current(2) - shaft_center(1);
        V = P_current(3) - shaft_center(2);
        radius = hypot(U, V);
        % A radial/tangential frame is undefined exactly on the shaft axis.
        if radius <= eps
            error('cylindrical local frame is undefined on the shaft axis.');
        end
        e_radial = [0, U/radius, V/radius];
        e_tangent = [0, -V/radius, U/radius];
        d_axis = rel(:,1);

    % Resolve the local basis around a Y-oriented shaft.
    case 'Y'
        U = P_current(1) - shaft_center(1);
        V = P_current(3) - shaft_center(2);
        radius = hypot(U, V);
        % A radial/tangential frame is undefined exactly on the shaft axis.
        if radius <= eps
            error('cylindrical local frame is undefined on the shaft axis.');
        end
        e_radial = [U/radius, 0, V/radius];
        e_tangent = [-V/radius, 0, U/radius];
        d_axis = rel(:,2);

    % Resolve the local basis around a Z-oriented shaft.
    case 'Z'
        U = P_current(1) - shaft_center(1);
        V = P_current(2) - shaft_center(2);
        radius = hypot(U, V);
        % A radial/tangential frame is undefined exactly on the shaft axis.
        if radius <= eps
            error('cylindrical local frame is undefined on the shaft axis.');
        end
        e_radial = [U/radius, V/radius, 0];
        e_tangent = [-V/radius, U/radius, 0];
        d_axis = rel(:,3);

    % Reject an invalid axis before projecting candidate offsets.
    otherwise
        error('cylinder_axis must be ''X'', ''Y'', or ''Z''.');
end

% Project global offsets onto the local orthonormal directions.
dz_signed = rel * e_radial(:);
d_tangent = rel * e_tangent(:);
center_xy = [d_axis, d_tangent];
dxy = hypot(d_axis, d_tangent);

end

function [idx_global, weights] = inverse3dFreeDragWeights(idx_track_so_far, P_track, L_melt)

%INVERSE3DFREEDRAGWEIGHTS Construct one causal exponential dragging stencil.
% Returned global indices reference post-remelting values on the current track;
% weights depend on physical distance, so they are independent of point density.

k = numel(idx_track_so_far);

% At the first track point the dragging operator is the identity.
if k == 1
    idx_global = idx_track_so_far(1);
    weights = 1;
    return
end

% Measure successive real-space distances along the deposited track.
segment_lengths = zeros(k,1);

% Compute the distance contributed by each newly added segment.
for q = 2:k
    segment_lengths(q) = norm(P_track(q,:) - P_track(q-1,:));
end

% Search backward for the oldest point with a non-negligible tail weight.
cumulative_distance = cumsum(segment_lengths);
j_oldest = k;

% Expand the retained history backward until its omitted tail is negligible.
while j_oldest > 1
    tail_distance = cumulative_distance(k) - cumulative_distance(j_oldest-1);
    tail_weight = exp(-tail_distance / L_melt);
    segment_weight = 1 - exp(-segment_lengths(j_oldest) / L_melt);

    % Truncate when the remaining upstream tail is at most 10% of the local
    % segment contribution, matching the direct forward implementation.
    if tail_weight / max(segment_weight, eps) <= 0.1
        break
    end

    j_oldest = j_oldest - 1;
end

% Retain only the active part of track history in this sparse stencil.
idx_global = idx_track_so_far(j_oldest:k);
idx_global = idx_global(:);
weights = zeros(numel(idx_global),1);

% The oldest retained value receives the residual exponential tail entering
% the finite window from all omitted upstream history.
distance_oldest = cumulative_distance(k) - cumulative_distance(j_oldest);
weights(1) = exp(-distance_oldest / L_melt);

% Newer values receive exact integrated exponential segment weights.
for j = (j_oldest+1):k
    local_pos = j - j_oldest + 1;
    distance_back = cumulative_distance(k) - cumulative_distance(j);
    coef = 1 - exp(-segment_lengths(j) / L_melt);
    weights(local_pos) = coef * exp(-distance_back / L_melt);
end

end

function fractions = inverse3dFreeExclusiveOverlapFractions( ...
    center_xy, dz_signed, w, h, b, n_samples)

%INVERSE3DFREEEXCLUSIVEOVERLAPFRACTIONS Partition remelted bead volume.
% Each returned fraction belongs exclusively to one previous candidate under
% the latest-bead-wins rule; the unclaimed remainder is the fresh-feed volume.

% Standardize array shapes for scalar and multi-candidate calls.
center_xy = reshape(center_xy, [], 2);
dz_signed = dz_signed(:);
n_candidates = size(center_xy,1);
fractions = zeros(n_candidates,1);

% No prior candidates imply no remelted volume contribution.
if n_candidates == 0
    return
end

% Deform deterministic unit-ball samples into the asymmetric bead ellipsoid.
unit_points = inverse3dFreeUnitBallSamples(n_samples);
a = w/2; % transverse semi-width [mm]
sample_xyz = zeros(size(unit_points));
sample_xyz(:,1:2) = a*unit_points(:,1:2);

% Apply separate upper-height and lower-penetration scales.
top_half = unit_points(:,3) >= 0;
sample_xyz(top_half,3) = h*unit_points(top_half,3);
sample_xyz(~top_half,3) = b*unit_points(~top_half,3);

% Weight by the appropriate Jacobian so the two differently scaled halves
% represent their correct physical volumes.
sample_weights = b*ones(size(unit_points,1),1);
sample_weights(top_half) = h;
total_weight = sum(sample_weights);
unclaimed = true(size(sample_weights));

% Assign every sample to at most one candidate, newest first; candidate order
% is chronological, so reverse traversal implements latest-bead ownership.
for j = n_candidates:-1:1
    dx = (sample_xyz(:,1)-center_xy(j,1))/a;
    dy = (sample_xyz(:,2)-center_xy(j,2))/a;
    z_relative = sample_xyz(:,3)-dz_signed(j);
    % Test whether each still-unclaimed sample lies inside candidate j.
    z_scale = b*ones(size(z_relative));
    z_scale(z_relative >= 0) = h;
    inside = ...
        dx.^2 + dy.^2 + (z_relative./z_scale).^2 <= 1 + 10*eps;
    claimed = inside & unclaimed;

    % Convert newly claimed weighted volume to a fraction of current bead.
    fractions(j) = sum(sample_weights(claimed)) / total_weight;
    unclaimed(claimed) = false;
end

% Renormalize only to remove a possible floating-point overshoot above one.
fraction_sum = sum(fractions);
% Correct only an overshoot; valid sub-unit sums retain fresh-feed volume.
if fraction_sum > 1
    fractions = fractions / fraction_sum;
end

end


function points = inverse3dFreeUnitBallSamples(n_samples)

%INVERSE3DFREEUNITBALLSAMPLES Reproducible low-discrepancy volume samples.
% Radical-inverse sequences replace random draws, ensuring the precomputed
% inverse operator and the direct forward model remain deterministic.

persistent cached_count cached_points

% Enforce a minimum integration resolution and reuse an identical cached set.
n_samples = max(1024, round(n_samples));
% Return the cached deterministic point set when its resolution matches.
if ~isempty(cached_count) && cached_count == n_samples
    points = cached_points;
    return
end

% Bases 2, 3, and 5 provide independent radius/polar/azimuth coordinates;
% the cubic root of the radius coordinate gives uniform ball-volume density.
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

% Cache because every geometric-overlap call uses the same point set.
cached_count = n_samples;
cached_points = points;

end

function values = inverse3dFreeRadicalInverse(sample_index, base)

%INVERSE3DFREERADICALINVERSE Compute a Van der Corput sequence in one base.
% Each pass consumes one integer digit and reflects it across the radix point.

values = zeros(size(sample_index));
factor = 1/base;

% Continue until no sample index has unprocessed base-b digits.
while any(sample_index > 0)
    values = values + factor*mod(sample_index,base);
    sample_index = floor(sample_index/base);
    factor = factor/base;
end

end

function has_overlap = inverse3dFreeBeadVerticalRangesOverlap(dz_signed, h, b, coord_tol)

%INVERSE3DFREEBEADVERTICALRANGESOVERLAP Test build-interval intersection.
% dz_signed locates each previous bead centre relative to the current centre;
% each bead occupies the interval [-b,+h] around its own centre.

previous_bottom = dz_signed - b;
previous_top = dz_signed + h;
current_bottom = -b;
current_top = h;

% Intersection exists when the common lower bound is no greater than the
% common upper bound, within the numerical coordinate tolerance.
z0 = max(current_bottom, previous_bottom);
z1 = min(current_top, previous_top);
has_overlap = z1 >= z0 - coord_tol;

end
