%% Mesh Pre-processing (MeshLab)
% Simpify point cloud:
% Filters-Cleaning and Repairing-Merge Close Vertices
% Filters-Sampling-Point Cloud Simplification-Number of Samples: 100000 OR Filters-Sampling-Poisson-disk sampling-Toggle Base mesh subsampling, Number of Samples: 100000
% Surface reconstruction:
% Filters-Remeshing, Simplification and Reconstruction-Screened Poisson Surface Reconstruction
% Clean mesh:
% Filters-Selection-Select faces with edge longer than...-Delete selected faces and vertices
% Filters-Selection-Small component selection-Delete selected faces and vertices

%% Mesh Post-processing (Scale calibration)
% How to get the real scale of scene in world units? 
% Option 1: Measure the distance between two camera locations (either
% manually or connect two cameras with know distance), then this distance
% actually becomes a "virtual" calibration ruler
% Option 2: Attach a motion sensor with camera that can record the moving distance between each photo

%% Libraries related:
% Graph & Mesh: https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph 
% Queue/Stack (Cqueue.m:96 bug fix): https://www.mathworks.com/matlabcentral/fileexchange/28922-list-queue-stack
addpath(genpath('toolbox_graph'));
addpath(genpath('datastructure'));

%% Control panel
close all;
name = 'clean_mesh';
global PLOT PLOT_FIG plot_mesh_original plot_mesh_curvature;
PLOT = true;
PLOT_FIG = 1;
    plot_mesh_original = 1;
    plot_mesh_curvature = 1;

%% Read Mesh and Pre-compute Mesh Info
global vertex faces nvertex nface face_rings vertex_rings face_normals face_centers face_colors;
% Denote V = No. of vertices, F = No. of faces
% vertex:           3 x V matrix for coordinates (x, y, z) of all V vertices, vertex(:, i)
% faces:            3 x F matrix for vertex index (I1, I2, I3) of all F faces, faces(:,i)
% face_rings:       length = F list, face_rings{i}: adjacent face indices (edge adjacency) of ith face
% vertex_rings:     length = F list, vertex_rings{i}: adjacent face indices (vertex adjacency) of ith face
% face_normals:     3 x F matrix, face_normals(:, i): normalized normal vector of ith face
% face_centers:     3 x F matrix, face_centers(:, i): center coordinates of ith face

tic

[vertex,faces] = read_mesh(name);

% Pre-compute Mesh Info
nvertex = size(vertex, 2);
nface = size(faces, 2);
face_rings = compute_face_ring(faces); 
vertex_rings{nface} = [];              
face_normals = zeros(size(faces));     
face_centers = zeros(size(faces));     

% Compute vertex rings (we can't accelerate this step because it's cell structure)
vertex_face_rings = compute_vertex_face_ring(faces); % length = V list, vertex_face_rings{i}: adjacent face indices of ith vertex
for i = 1 : nface
    face = faces(:, i);
    all_neighbors = unique([vertex_face_rings{[face(1) face(2) face(3)]}]);
    vertex_rings{i} = all_neighbors(all_neighbors ~= i); % exclude the current face
end

% Compute face normals and face centers (accelerated by matrix operation)
v1 = vertex(:, faces(1, :));
v2 = vertex(:, faces(2, :));
v3 = vertex(:, faces(3, :));
% Normal vectors
face_normals = cross(v2 - v1, v3 - v1, 1); % cross product along dimension = 1, assume a counter-clockwise order convention of vertices, then the cross product should also be counter-clockwise, v2-v1-->v3-v1, right-hand principle
face_normals = face_normals ./ vecnorm(face_normals); % normalize vector along dimension = 1
% Center coordinates
face_centers = (v1 + v2 + v3) / 3;

fprintf('Pre-compute mesh: %f seconds\n', toc);

global_list = who('global');
clearvars('-except', global_list{:});

%% Display original mesh
if PLOT && plot_mesh_original
    face_colors = 0.3 * ones(nface, 1); % 0.3-red
    plot_face_color('Original Mesh', 0);
end

%% Compute face angle
tic
global face_angles;
K = 3; % order of neighbor search ring
face_angles = compute_face_angle(K);
face_angles = (face_angles - min(face_angles)) ./ (max(face_angles) - min(face_angles));
fprintf('Compute face angle: %f seconds\n', toc);

%% Display curvature mesh
if PLOT && plot_mesh_curvature
    face_colors = face_angles;
    plot_face_color('Curvature Mesh', 0);
end

%% Breadth-First Search Traversal
% Label all faces into one of the following categories in array 'state':
% Cat 0: Unvisited face (black-colored)
% Cat 1: Object face    (red-colored)
% Cat 2: Boundary face  (white-colored)
state = zeros(size(faces, 2), 1); % array for recording face status: 0-unvisited face(black), 1-object face(red), 2-boundary face(white)
q = CQueue();                     % queue for BFS, storing the face index, i.e. ith face. Any face that entered into this queue is labelled 'object face'.

% Good seeds
seed = 2400;
% bad seed
seed = 42885;
rng('shuffle');
%seed = randi(size(faces, 2));

% clean_mesh, local jiayi = 0.8, face angle jiayi 3rd neighbors
% good seed
seed = 10767; % 4
seed = 49500; % 5
seed = 42885; % 15
seed = 46018; % 7
seed = 63104; % 3
seed = 36193; % 17
seed = 66780; % dabian
% seed = 5620; % xiaobian
% seed = 62000; % xiaoxiaobian
% seed = 3674; % shuzhede

% Thresholds and Controls
threshold1 = 0.6;         % local criterion for convexity/curvature (adjustable, and since both vectors are normalized, this threshold value is actually a theta angle that you can specify)
threshold2 = 0.15;         % global criterion for centroid facing (dynamically adjustable, this criteria doesn't start functioning much until a significant amount of object faces are identified b/c the calculated centroid rely on the current faces are sufficient to form a particle body)
global_trigger = false;   % flag for start applying the global criterion
global_trigger_val = 0;  % trigger 'center_change_ratio(%)' value for the global criterion. 10 means the center change drops below 10%, which shows stabilize
% Deciding when to start applying global criterion is still vague:
% maybe, developing a metric of "percentage of shape completion" can help, 
% and will also help the later completion step; 
% maybe, a surface/volume ratio metric can help; 
% maybe, the sum/avg of perpendicular distance from centroid to all object faces can help
neighboring = 1; % 1-face neighboring (edge adjacency), 2-vertex neighboring (vertex adjacency)

% Initialize
placeholder = 0; % BFS placeholder (a face index can never be 0 b/c Matlab is 1-based)
% seed = randi(size(faces, 2)); % or push a random starting face index % TODO: what if the start face is a boundary face? during BFS we can add a check on all 3 adjacent faces (no matter visited or unvisited) and if the face contradicts with all adjacent faces, remove it from the object face set
q.push(seed); q.push(placeholder);
state(seed) = 1; % label starting face as 'object'
centroid = face_centers(:, seed); % initialize dynamic centroid as current face
update_times = 0; % record times of centroid updates
center_sum = zeros(3, 1); face_count = 0; % used for accumulating face coords between update loops
center_change_sum = zeros(3, 1); % used for accumulating the denominator of the center change ratio
% Dynamic centroid of all current object faces is updated asynchronously:
% When in a normal loop, 'center_sum' accumulates every visited center coords; 'face_count' accumulates the count of object faces.
% When in a placeholder loop, 'centroid' is updated by 'center_sum/face_count'

while q.isempty() ~= 1
    % Pop a face from the queue
    curr_id = q.pop();
    
    % If a placeholder is met, update the current centroid; otherwise, just keep accumulating
    if curr_id == placeholder
        if update_times > 0 % skip the first placeholder (centroid is just the starting seed face, no change)
            centroid_old = centroid;
            centroid = center_sum / face_count;
            center_change_sum = center_change_sum + abs(centroid - centroid_old);
            center_change_ratio = norm(centroid - centroid_old) / norm(center_change_sum) * 100;
            center_change_ratios(update_times) = center_change_ratio; % 'center_change_ratios' records the change ratio at each update
            if center_change_ratio < global_trigger_val
                global_trigger = true; % trigger global criterion
            end
        end
        update_times = update_times + 1;
        if q.isempty() ~= 1
            q.push(placeholder); % if there are unvisited faces in queue, push the placeholder
        else
            break; % otherwise, this placeholder is the last element, exit (otherwise you will infinite loop)
        end
        continue; % for placeholder step, just skip to the next loop after updating
    end
    
    % Accumulate the dynamic centroid
    curr_normal = - face_normals(:, curr_id); % flip sign, now normal points to the interior of a particle
    curr_center = face_centers(:, curr_id);
    center_sum = center_sum + curr_center;
    face_count = face_count + 1;
    
    % Label neighboring faces into object face (1) OR boundary face (2) 
    if neighboring == 1
        neighbors = face_rings{curr_id};
    else
        neighbors = vertex_rings{curr_id};
    end
    for f = 1 : length(neighbors)
        adj_id = neighbors(f);
        adj_normal = - face_normals(:, adj_id); % flip sign, now normal points to the interior of a particle
        adj_center = face_centers(:, adj_id); 
        if state(adj_id) == 0 % skip all visited faces
            % Criteria 1: Local curvature criteria (between a face and its adjacent faces)
            local_center_diff = adj_center - curr_center;
            local_center_diff = local_center_diff / norm(local_center_diff); % normalize
            local_measure = dot(curr_normal, local_center_diff + adj_normal);
            % Criteria 2: Global centroid criteria (between a face and the global centroid)
            global_center_diff = adj_center - centroid;
            global_measure = dot(adj_normal, global_center_diff); % Attempt 1: dot product without normalization
            global_center_diff = global_center_diff / norm(global_center_diff);
            global_measure = dot(-adj_normal, global_center_diff); % Attempt 2: projection (global_center_diff is normalized, so |b|*cos(theta) = projection, the threshold criteria is also essentially a theta angle)

            if face_angles(adj_id) >= 0.8 && (global_trigger == false || global_measure > threshold2)
            %if local_measure > threshold1 && (global_trigger == false || global_measure > threshold2) % if trigger is false, ignore the global criterion; otherwise check it
                % object face (push to queue)
                state(adj_id) = 1;
                q.push(adj_id); 
            else
                % boundary face (do not push to queue, BFS ends at this face)
                state(adj_id) = 2;
                % label its adjacency as well
                if neighboring == 1
                    boundary_neighbors = face_rings{adj_id};
                else
                    boundary_neighbors = vertex_rings{adj_id};
                end
                for b = 1 : length(boundary_neighbors)
                    boundary_id = boundary_neighbors(b);
                    if state(boundary_id) == 0
                        state(boundary_id) = 2;
                    end
                end
            end
        end
    end
end

% Display raw mesh
% Per-face coloring: Non-visited face-0 (black), object vertex-0.3 (red), boundary vertex-1.0 (white)
% close all;
% figure(2); 
% title('Raw Mesh');
% face_colors = zeros(size(faces, 2), 1);
% objects = state == 1; face_colors(objects) = 0.3;
% boundaries = state == 2; face_colors(boundaries) = 1.0;
% DISPLAY_HHH = false;
% if DISPLAY_HHH
%     % My manual display settings
%     h = patch('vertices',vertex','faces',faces','FaceVertexCData', face_colors, 'FaceColor', 'flat');
%     lighting phong;
%     camproj('perspective');
%     axis square; 
%     axis off;
%     %cameramenu;
%     cameratoolbar;
%     axis tight;
%     axis equal;
%     shading faceted;
%     colormap hot;
%     caxis([0 1]); % fix colormap range
%     camlight;
% else
%     % Use graph toolbox display settings
%     options.face_vertex_color = face_colors;
%     plot_mesh(vertex, faces, options);
%     shading faceted; % or shading interp; % for display smooth surface
%     colormap hot;
%     caxis([0 1]); % fix colormap range
% end
% 
% % Plot normal for a specific face index
% hold on;
% id = seed;
% quiver3(face_centers(1,id), face_centers(2,id), face_centers(3,id), face_normals(1,id), face_normals(2,id), face_normals(3,id), 'r'); % draw normal

% Mesh Cleaning by Detecting connected components in unvisited faces (BFS)
state2 = zeros(size(faces, 2), 1);
state2(state ~= 0) = 1; % 0-unvisited, 1-visited
connect = 0; % No. of connected components
while true
    unvisited = find(state2 == 0); % index of all unvisited faces
    if (length(unvisited) == 0) 
        break; % if no more unvisited faces, exit
    end
    starter = unvisited(randi(length(unvisited)));
    face_set = zeros(size(faces, 2), 1);
    q.empty();
    q.push(starter);
    state2(starter) = 1;
    face_set(starter) = 1; % don't forget!
    count = 0;
    while q.isempty() ~= 1
        curr_id = q.pop();
        count = count + 1;
        neighbors = face_rings{curr_id};
        for f = 1 : length(neighbors)
            adj_id = neighbors(f);
            if state2(adj_id) == 0
                state2(adj_id) = 1;
                face_set(adj_id) = 1;
                q.push(adj_id);
            end
        end
    end
    connect = connect + 1;
    component{connect} = face_set; % component{i} is the logical face index array of the ith component
    fprintf("Connected Component %d: %d faces, starter %d\n", connect, count, starter);
end

% Find the connected component that includes the most faces
max_set = 1;
max_face = 0;
for c = 1 : connect
    if length(component{c}) > max_face
        max_set = c;
        max_face = length(component{c});
    end
end

% Clean the mesh by taking the complement of the largest component
objects = component{max_set} == 0;

% Extract boundary faces from the cleaned mesh
boundaries = boundary_extract(faces, objects);

% Display cleaned mesh
figure(2), clf;
title('Cleaned Mesh');
face_colors = zeros(size(faces, 2), 1);
face_colors(objects) = 0.3;
face_colors(boundaries) = 1.0;
options.face_vertex_color = face_colors;
plot_mesh(vertex, faces, options);
shading faceted;  
colormap hot;

% Optimize boundary (shortcut)
E_local = [];
E_global = [];
for i = 1 : length(boundaries)
    neighbors = vertex_rings{boundaries(i)};
    neighbors = intersect(neighbors, boundaries);
    temp = [repmat(boundaries(i),length(neighbors),1) neighbors];
    E_local = [E_local; temp];
end

object_faces = find(objects == 1); % edge for object faces only
% object_faces = 1 : length(faces); % edge for all faces
for i = 1 : length(object_faces)
    neighbors = vertex_rings{object_faces(i)};
    temp = [repmat(object_faces(i),length(neighbors),1) neighbors'];
    E_global = [E_global; temp];   
end

num_change = 3;
total_faces = length(object_faces);
for i = 1 : num_change
    boundaries = optimize_with_direct_dist(E_local, E_global, face_centers', boundaries', 1); 
    % BFS
    visited = zeros(size(faces, 2), 1);
    fence = zeros(size(faces, 2), 1);
    q.empty();
    q.push(seed);
    visited(seed) = 1;
    fence(boundaries) = 1;
    while q.isempty() ~= 1
        curr_id = q.pop();
        neighbors = face_rings{curr_id};
        for f = 1 : length(neighbors)
            adj_id = neighbors(f);
            if visited(adj_id) == 0 % skip all visited faces
                if fence(adj_id) == 0
                    % object face (push to queue)
                    visited(adj_id) = 1;
                    q.push(adj_id); 
                end
            end
        end
    end
    objects_temp = visited == 1;
    if sum(objects_temp(:)) <= 0.6 * total_faces % should add seed
        break;
    end
    objects = objects_temp;
    boundaries = boundary_extract(faces, objects);
    
    % Update E_local
    E_local = [];
    for j = 1 : length(boundaries)
        neighbors = vertex_rings{boundaries(j)};
        neighbors = intersect(neighbors, boundaries);
        temp = [repmat(boundaries(j),length(neighbors),1) neighbors];
        E_local = [E_local; temp];
    end
end

% Display optimized mesh
% Display optimized boundary
figure(3);
title('Optimized Mesh_0');
face_colors = zeros(size(faces, 2), 1);
face_colors(objects) = 0.3;
face_colors(boundaries) = 1.0;
options.face_vertex_color = face_colors;
plot_mesh(vertex, faces, options);
shading faceted;  
colormap hot;
caxis([0 1]);