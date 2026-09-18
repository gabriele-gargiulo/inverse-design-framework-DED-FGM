%clear; clc;

%% FORWARD 3D MIXING MODEL
% Computes the deposited composition along a 3D toolpath. The
% feed composition is sampled from the nearest material point, then updated
% with geometric remelting and in-track dragging.
%
% Model sequence:
%   1. Map the prescribed feed-composition field to each toolpath point.
%   2. Recover layers and continuous tracks from the ordered toolpath.
%   3. At every deposition point, partition the remelted volume among all
%      eligible previously deposited beads using a latest-bead-wins rule.
%   4. Apply exponential in-track dragging over the moving melt-pool length.
%   5. Restore the original CSV row order and save the predicted composition.
%
% Each row of a composition array contains the fractions [IN, MK, INC].
% Coordinates and geometric process parameters are expressed in millimetres.

%% INPUT FILES

input_material_file = 'input_material_distribution_3d_blisk.csv'; % spatial feed-composition field
toolpath_file = 'toolpath_3d_blisk.csv';                          % ordered deposition coordinates

% Import both CSV files as tables so columns can be addressed by name.
mat_data = readtable(input_material_file);
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

%% EXTRACT POINTS AND FEED COMPOSITION

% Material-field coordinates at which the prescribed feed is defined.
X_mat = mat_data.X;
Y_mat = mat_data.Y;
Z_mat = mat_data.Z;

% The first two material columns are required by the input-file format.
IN_mat = mat_data.IN_frac;
MK_mat = mat_data.MK_frac;

% Permit legacy two-material files by creating a zero INC fraction.
if ismember('INC_frac', mat_data.Properties.VariableNames)
    INC_mat = mat_data.INC_frac;
% Missing third-material data represent a two-material composition field.
else
    INC_mat = zeros(height(mat_data),1);
end

% Toolpath rows are already ordered chronologically by deposition time.
X_path = toolpath.X;
Y_path = toolpath.Y;
Z_path = toolpath.Z;
path_points_real = [X_path, Y_path, Z_path];
material_points_real = [X_mat, Y_mat, Z_mat];

%% POINT CLOUD

% Every toolpath row is one model/deposition point in the final point cloud.
free_points = path_points_real;
N = size(free_points,1); % total number of deposition points

% Build an unwrapped metric coordinate system when modelling a cylindrical
% part. The three model coordinates are axial distance, arc length, radius.
if strcmpi(coordinate_mode, "cylindrical")

    [A_mat, Theta_mat, R_mat] = freeCartesianToCylindrical( ...
        X_mat, Y_mat, Z_mat, cylinder_axis, shaft_center, theta_zero);
    [A_path, Theta_path, R_path] = freeCartesianToCylindrical( ...
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
% straight, diagonal, and curved deposition paths.

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

% Assign a feed-field row to every path point, but search only within the
% corresponding layer so neighbouring layers cannot be confused spatially.
material_build_coord = material_points_model(:,3);
idx_feed_material = zeros(N,1);

% Process each layer independently during composition-field sampling.
for k = 1:n_layers_path
    idx_path_layer = find(layer_of_path_point == k);
    build_coord_layer = build_coord_path(idx_path_layer);
    layer_min = min(build_coord_layer) - coord_tol;
    layer_max = max(build_coord_layer) + coord_tol;
    idx_material_layer = find( ...
        material_build_coord >= layer_min & ...
        material_build_coord <= layer_max);

    % If coordinate noise prevents an exact layer-range match, accept all
    % material points within half a nominal bead height of the layer centre.
    if isempty(idx_material_layer)
        layer_center = mean(build_coord_layer);
        idx_material_layer = find(abs(material_build_coord - layer_center) <= h/2 + coord_tol);
    end

    % Stop with a specific layer number rather than silently using bad data.
    if isempty(idx_material_layer)
        error('No feed material points found on toolpath layer %d.', k);
    end

    % Use a seam-periodic nearest-neighbour lookup for a full cylinder.
    if strcmpi(coordinate_mode, "cylindrical")
        [idx_local, ~] = freeNearestToolpathCylindrical( ...
            material_points_model(idx_material_layer,:), ...
            free_points_model(idx_path_layer,:), ...
            R_ref, full_circumference);
    % Cartesian matching uses ordinary Euclidean nearest neighbours.
    else
        [idx_local, ~] = knnsearch( ...
            material_points_model(idx_material_layer,:), ...
            free_points_model(idx_path_layer,:));
    end

    idx_feed_material(idx_path_layer) = idx_material_layer(idx_local);
end

% Gather the three-component feed vector associated with every path point.
C_material_feed = [IN_mat, MK_mat, INC_mat];
C_feed = C_material_feed(idx_feed_material,:);

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

% Transfer the metadata from toolpath rows to the model point-cloud rows.
layer_of_free_point = layer_of_path_point(idx_toolpath);

track_id = path_track_id(idx_toolpath);

% This fallback prevents invalid zero indexing if a malformed path escaped
% the segmentation logic; normally every point already has a track ID.
unassigned_track = track_id == 0;
% Apply the defensive fallback only when at least one zero ID is present.
if any(unassigned_track)
    track_id(unassigned_track) = 1;
end

%% BUILD TRUE PRINT ORDER FROM TOOLPATH

% Sort primarily by toolpath row and secondarily by distance along that row.
% For the current one-point-per-row representation this preserves CSV order.
sort_matrix = [idx_toolpath(:), dist_tool(:)];
[~, sort_idx] = sortrows(sort_matrix, [1 2]);

free_points_print = free_points_model(sort_idx,:);
free_points_real_print = free_points(sort_idx,:);
C_feed_print = C_feed(sort_idx,:);
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

% Cache point membership, local position, and cumulative physical distance
% for each track; dragging later uses distance rather than point count.
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

%% TOOLPATH-DRIVEN POINT MIXING

% C_standard is the composition immediately after remelting; C_real is the
% final composition after subsequent in-track dragging at the same point.
C_real_print = zeros(N,3);
C_standard_print = zeros(N,3);

% Diagnostic arrays retain the most recent lateral/vertical contributor,
% their total fractions, and the number of contributors of each type.
idx_lat_print = zeros(N,1);
idx_vert_print = zeros(N,1);
alpha_lat_print = zeros(N,1);
alpha_vert_print = zeros(N,1);
n_lat_neighbors = zeros(N,1);
n_vert_neighbors = zeros(N,1);

ellipsoid_width = w;      % total transverse footprint used for broad search
ellipsoid_height = h + b; % full top-to-bottom bead extent
neighbor_search_radius = sqrt(ellipsoid_width^2 + ellipsoid_height^2);

% Use the Statistics Toolbox range search when available. The empty-cell
% fallback below performs an exhaustive search and gives identical physics.
if exist('rangesearch','file') == 2
    neighbor_candidates_all = freeRangeSearchCandidates( ...
        free_points_print, neighbor_search_radius, ...
        coordinate_mode, R_ref, full_circumference);
% Without rangesearch, an empty marker triggers exhaustive prior-point scans.
else
    neighbor_candidates_all = {};
end

% March forward through print time because remelting is causal: point i may
% depend only on already printed points 1,...,i-1.
for i = 1:N

    point_curr = free_points_print(i,:);
    C_f = C_feed_print(i,:);

    % Only already-printed points can remelt the current point.
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
    % Exclude nearby samples from the same continuous melt pool. Their effect
    % is represented once by the dragging filter, avoiding double counting.
    if ~isempty(prev_idx)
        same_track = track_id_print(prev_idx) == track_id_print(i);
        distance_behind = track_distance_print(i) - track_distance_print(prev_idx);
        prev_idx(same_track & distance_behind < w) = [];
    end

    % With no eligible overlap, the newly deposited composition equals feed.
    C_standard_print(i,:) = C_f;

    % Evaluate remelting only if at least one earlier point remains.
    if ~isempty(prev_idx)

        % Express prior-point offsets in the local bead frame. Cylindrical
        % mode uses axial/tangential/radial distances at the current point.
        if strcmpi(coordinate_mode, "cylindrical")
            [dxy, dz_signed, center_xy] = freeCylindricalLocalDistances( ...
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

        % Broad geometric eligibility: transverse proximity plus overlap of
        % the asymmetric vertical bead intervals [-b,h].
        meltpool_intersects = ...
            dxy <= ellipsoid_width + coord_tol & ...
            freeBeadVerticalRangesOverlap(dz_signed, h, b, coord_tol);

        same_layer = layer_print(prev_idx) == layer_print(i);
        lower_layer = layer_print(prev_idx) < layer_print(i);
        eligible = meltpool_intersects & (same_layer | lower_layer);

        % Partition the current bead only when at least one candidate passes.
        if any(eligible)
            candidates = prev_idx(eligible);
            candidate_centers = center_xy(eligible,:);
            candidate_dz = dz_signed(eligible);

            % Chronological order is required by the latest-bead-wins rule.
            % Oldest-to-newest ordering lets the partition routine traverse
            % in reverse so newly deposited material claims shared volume.
            [candidates, order] = sort(candidates);
            candidate_centers = candidate_centers(order,:);
            candidate_dz = candidate_dz(order);

            % Each weight is the exclusive fraction of the current bead
            % occupied by one prior bead after latest-bead-wins partitioning.
            weights = freeExclusiveOverlapFractions( ...
                candidate_centers, candidate_dz, w, h, b, ...
                overlap_partition_samples);
            positive = weights > 0;

            candidates = candidates(positive);
            weights = weights(positive);

            % Blend prior realized compositions with the unremelted fresh feed.
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

                % Save the most recently printed lateral contributor.
                if any(lat_mask)
                    idx_lat_print(i) = candidates(find(lat_mask,1,'last'));
                end
                % Save the most recently printed lower-layer contributor.
                if any(vert_mask)
                    idx_vert_print(i) = candidates(find(vert_mask,1,'last'));
                end
            end
        end
    end

    % Apply the causal exponential dragging filter to the current track from
    % its first retained upstream point through the present point.
    track_idx = track_indices{track_group_print(i)};
    kk = position_in_track(i);
    idx_track_so_far = track_idx(1:kk);
    C_real_print(i,:) = freeApplyTrackDraggingAtPoint( ...
        C_standard_print(idx_track_so_far,:), ...
        free_points_real_print(idx_track_so_far,:), ...
        L_melt);
end

%% RESTORE ORIGINAL INPUT ORDER

% Undo the print-order permutation so output rows match the input toolpath.
C_real = zeros(size(C_real_print));
C_real(sort_idx,:) = C_real_print;

%% SAVE OUTPUT

% Retain the original physical coordinates, including for cylindrical runs.
X_free = free_points(:,1);
Y_free = free_points(:,2);
Z_free = free_points(:,3);

% Round only the saved composition values; all internal calculations use the
% full floating-point precision accumulated by the model.
output_table = table( ...
    X_free, Y_free, Z_free, ...
    round(C_real(:,1),3), round(C_real(:,2),3), round(C_real(:,3),3), ...
    'VariableNames', ...
    {'X','Y','Z','IN_frac','MK_frac','INC_frac'});

% The fixed output name is consumed directly by the deviation plotting script.
writetable(output_table, 'output_material_distribution_3d_free.csv');

disp('Free-point 3D material distribution saved.');

%% LOCAL FUNCTIONS

function fractions = freeExclusiveOverlapFractions( ...
    center_xy, dz_signed, w, h, b, n_samples)

%FREEEXCLUSIVEOVERLAPFRACTIONS Partition the current bead's remelted volume.
% center_xy(j,:) and dz_signed(j) locate prior candidate j relative to the
% current bead. The returned fraction for each candidate is exclusive: a
% sample already claimed by a newer candidate cannot be claimed again by an
% older one. Consequently, the fractions sum to at most one.

% Standardize shapes so scalar and multi-candidate calls behave identically.
center_xy = reshape(center_xy, [], 2);
dz_signed = dz_signed(:);
n_candidates = size(center_xy,1);
fractions = zeros(n_candidates,1);

% No candidate means no remelted contribution.
if n_candidates == 0
    return
end

% Generate deterministic samples in a unit ball, then deform the ball into
% an asymmetric ellipsoid with upper height h and lower penetration b.
unit_points = freeUnitBallSamples(n_samples);
a = w/2; % semi-width of the ellipsoidal cross-section [mm]
sample_xyz = zeros(size(unit_points));
sample_xyz(:,1:2) = a*unit_points(:,1:2);

% Scale the positive and negative halves independently in the build direction.
top_half = unit_points(:,3) >= 0;
sample_xyz(top_half,3) = h*unit_points(top_half,3);
sample_xyz(~top_half,3) = b*unit_points(~top_half,3);

% The asymmetric z scaling has a different Jacobian in each half; weighting
% samples by h or b makes the volume estimate unbiased between the two halves.
sample_weights = b*ones(size(unit_points,1),1);
sample_weights(top_half) = h;
total_weight = sum(sample_weights);
unclaimed = true(size(sample_weights));

% Assign every sample to at most one candidate, newest first. Candidate order
% was made chronological by the caller, so reverse iteration implements the
% latest-bead-wins material ownership rule.
for j = n_candidates:-1:1
    dx = (sample_xyz(:,1)-center_xy(j,1))/a;
    dy = (sample_xyz(:,2)-center_xy(j,2))/a;
    z_relative = sample_xyz(:,3)-dz_signed(j);
    % Test membership in candidate j's own asymmetric ellipsoid.
    z_scale = b*ones(size(z_relative));
    z_scale(z_relative >= 0) = h;
    inside = ...
        dx.^2 + dy.^2 + (z_relative./z_scale).^2 <= 1 + 10*eps;
    claimed = inside & unclaimed;

    % Convert newly claimed weighted volume to a fraction of the current bead.
    fractions(j) = sum(sample_weights(claimed)) / total_weight;
    unclaimed(claimed) = false;
end

% Protect against tiny floating-point overshoots above unit total volume.
fraction_sum = sum(fractions);
% Correct only an overshoot; valid sub-unit sums leave fresh-feed volume.
if fraction_sum > 1
    fractions = fractions / fraction_sum;
end

end

function points = freeUnitBallSamples(n_samples)

%FREEUNITBALLSAMPLES Return reproducible low-discrepancy 3D volume samples.
% Radical-inverse sequences in bases 2, 3, and 5 replace random sampling, so
% repeated runs produce the same overlap fractions and output compositions.

persistent cached_count cached_points

% Enforce a useful minimum integration resolution.
n_samples = max(1024, round(n_samples));
% Reuse the previous point set when this function is called at the same size.
if ~isempty(cached_count) && cached_count == n_samples
    points = cached_points;
    return
end

% Independent low-discrepancy coordinates supply radius, polar angle, and
% azimuth. Taking the cubic root yields uniform density over ball volume.
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

% Cache the deterministic result because every overlap call uses it again.
cached_count = n_samples;
cached_points = points;

end

function values = freeRadicalInverse(sample_index, base)

%FREERADICALINVERSE Compute the Van der Corput radical-inverse sequence.
% Each loop consumes one base-b digit and reflects it across the radix point.

values = zeros(size(sample_index));
factor = 1/base;

% Continue until every integer sample index has no remaining base-b digits.
while any(sample_index > 0)
    values = values + factor*mod(sample_index,base);
    sample_index = floor(sample_index/base);
    factor = factor/base;
end

end

function has_overlap = freeBeadVerticalRangesOverlap(dz_signed, h, b, coord_tol)

%FREEBEADVERTICALRANGESOVERLAP Test build-direction interval intersection.
% dz_signed is the previous bead centre relative to the current bead centre;
% both beads span from -b below their centre to +h above their centre.

previous_bottom = dz_signed - b;
previous_top = dz_signed + h;
current_bottom = -b;
current_top = h;

% The two intervals overlap when their common lower bound does not exceed
% their common upper bound, allowing the stated numerical tolerance.
z0 = max(current_bottom, previous_bottom);
z1 = min(current_top, previous_top);
has_overlap = z1 >= z0 - coord_tol;

end

function [A, Theta, R] = freeCartesianToCylindrical(X, Y, Z, cylinder_axis, shaft_center, theta_zero)

%FREECARTESIANTOCYLINDRICAL Convert global XYZ to axial/angular/radial form.
% A is the coordinate along the selected shaft axis. U and V span its normal
% plane, Theta is wrapped to [0,2*pi), and R is distance from the shaft axis.

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
    % Reject misspelled axes before returning physically meaningless geometry.
    otherwise
        error('cylinder_axis must be ''X'', ''Y'', or ''Z''.');
end

% Rotate the angular origin to theta_zero and wrap across the chosen seam.
Theta = mod(atan2(V, U) - theta_zero, 2*pi);
R = hypot(U, V);

end

function [idx_nearest, dist_nearest] = freeNearestToolpathCylindrical( ...
    coords_path_model, coords_free_model, reference_radius, full_circumference)

%FREENEARESTTOOLPATHCYLINDRICAL Find nearest points in unwrapped coordinates.
% For a complete circumference, duplicate reference points one angular period
% to either side so points adjacent across the seam remain nearest neighbours.

% Add periodic images only when the part wraps through a full revolution.
if full_circumference
    period = 2*pi*reference_radius;
    n_path = size(coords_path_model,1);
    coords_aug = [coords_path_model; ...
                  coords_path_model + [zeros(n_path,1), period*ones(n_path,1), zeros(n_path,1)]; ...
                  coords_path_model - [zeros(n_path,1), period*ones(n_path,1), zeros(n_path,1)]];
    idx_aug_to_path = [(1:n_path)'; (1:n_path)'; (1:n_path)'];

    % Map augmented-image indices back to the original reference rows.
    [idx_aug, dist_nearest] = knnsearch(coords_aug, coords_free_model);
    idx_nearest = idx_aug_to_path(idx_aug);
% Partial cylinders have no periodic seam and use direct nearest neighbours.
else
    [idx_nearest, dist_nearest] = knnsearch(coords_path_model, coords_free_model);
end

end

function neighbor_candidates_all = freeRangeSearchCandidates( ...
    coords_model, search_radius, coordinate_mode, reference_radius, full_circumference)

%FREERANGESEARCHCANDIDATES Precompute geometrically nearby point indices.
% This is a performance filter only: exact asymmetric-bead intersection and
% layer eligibility are evaluated later in the main marching loop.

% Periodic image copies make the unwrapped cylindrical seam transparent.
if strcmpi(coordinate_mode, "cylindrical") && full_circumference
    period = 2*pi*reference_radius;
    n_points = size(coords_model,1);
    theta_shift = [zeros(n_points,1), period*ones(n_points,1), zeros(n_points,1)];

    coords_aug = [coords_model; ...
                  coords_model + theta_shift; ...
                  coords_model - theta_shift];
    idx_aug_to_point = [(1:n_points)'; (1:n_points)'; (1:n_points)'];

    % Search all images, then translate results back to original point IDs.
    [idx_aug_all, ~] = rangesearch(coords_aug, coords_model, search_radius);
    neighbor_candidates_all = cell(size(idx_aug_all));

    % Stable uniqueness removes duplicate images without reordering print IDs.
    for q = 1:numel(idx_aug_all)
        idx_mapped = idx_aug_to_point(idx_aug_all{q});
        neighbor_candidates_all{q} = unique(idx_mapped(:), 'stable');
    end
% Use the original coordinate set when periodic images are unnecessary.
else
    % Cartesian or partial-cylinder cases require no periodic duplication.
    [neighbor_candidates_all, ~] = rangesearch( ...
        coords_model, coords_model, search_radius);

    % Store every candidate list as a column vector for consistent indexing.
    for q = 1:numel(neighbor_candidates_all)
        neighbor_candidates_all{q} = neighbor_candidates_all{q}(:);
    end
end

end

function [dxy, dz_signed, center_xy] = freeCylindricalLocalDistances( ...
    P_previous, P_current, cylinder_axis, shaft_center)

%FREECYLINDRICALLOCALDISTANCES Resolve offsets in the current local frame.
% center_xy contains axial and tangential offsets, dxy is their resultant,
% and dz_signed is radial/build-direction offset. The local basis changes
% around the cylinder and must therefore be evaluated at P_current.

% Vector from the current deposition point to every previous candidate.
rel = bsxfun(@minus, P_previous, P_current);

% Construct radial and tangential unit vectors for the selected shaft axis.
switch upper(cylinder_axis)
    % Resolve the local basis around an X-oriented shaft.
    case 'X'
        U = P_current(2) - shaft_center(1);
        V = P_current(3) - shaft_center(2);
        radius = hypot(U, V);
        % A radial/tangential frame cannot be defined exactly on the axis.
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
        % A radial/tangential frame cannot be defined exactly on the axis.
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
        % A radial/tangential frame cannot be defined exactly on the axis.
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

% Project global offsets onto the orthonormal local directions.
dz_signed = rel * e_radial(:);
d_tangent = rel * e_tangent(:);
center_xy = [d_axis, d_tangent];
dxy = hypot(d_axis, d_tangent);

end

function C_filtered_current = freeApplyTrackDraggingAtPoint(C_standard_track, P_track, L_melt)

%FREEAPPLYTRACKDRAGGINGATPOINT Apply the causal exponential dragging model.
% The filter combines upstream post-remelting compositions according to their
% physical distance behind the current point. L_melt is the exponential decay
% length, so the result is independent of toolpath sampling density.

k = size(C_standard_track,1);

% At the first point there is no upstream melt-pool history to mix.
if k == 1
    C_filtered_current = C_standard_track(1,:);
    return
end

% Measure each toolpath segment in real Cartesian space.
segment_lengths = zeros(k,1);

% Compute the distance contributed by each newly added segment.
for q = 2:k
    segment_lengths(q) = norm(P_track(q,:) - P_track(q-1,:));
end

% Search backward for the oldest upstream point with a non-negligible weight.
cumulative_distance = cumsum(segment_lengths);
j_oldest = k;

% Expand the retained history backward until its omitted tail is negligible.
while j_oldest > 1
    tail_distance = cumulative_distance(k) - cumulative_distance(j_oldest-1);
    tail_weight = exp(-tail_distance / L_melt);
    segment_weight = 1 - exp(-segment_lengths(j_oldest) / L_melt);

    % Stop when the remaining upstream tail is at most 10% of the local
    % segment contribution; this truncates an otherwise infinite memory.
    if tail_weight / max(segment_weight, eps) <= 0.1
        break
    end

    j_oldest = j_oldest - 1;
end

% The oldest retained point represents all composition entering the retained
% window from farther upstream and therefore receives the residual tail weight.
C_filtered_current = zeros(1,size(C_standard_track,2));
distance_oldest = cumulative_distance(k) - cumulative_distance(j_oldest);
C_filtered_current = C_filtered_current + ...
    exp(-distance_oldest / L_melt) * C_standard_track(j_oldest,:);

% Add the exact integrated exponential contribution of every newer segment.
for j = (j_oldest+1):k
    distance_back = cumulative_distance(k) - cumulative_distance(j);
    coef = 1 - exp(-segment_lengths(j) / L_melt);
    weight = coef * exp(-distance_back / L_melt);
    C_filtered_current = C_filtered_current + weight * C_standard_track(j,:);
end

end
