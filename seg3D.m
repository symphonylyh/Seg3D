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
% Q1. How to get the real scale of scene in world units? 
% Option 1: Measure the distance between two camera locations (either
% manually or connect two cameras with know distance), then this distance
% actually becomes a "virtual" calibration ruler
% Option 2: Attach a motion sensor with camera that can record the moving distance between each photo

% Q2. How to estimate the "percentage of shape completion" of a particle?
% Option 1: poll on normal vector directions and fill a sphere surface (details in ground removal algorithm)
% Option 2: surface/volume ratio

%% Control panel
% Libraries related:
% Graph & Mesh: https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph 
% Queue/Stack (Cqueue.m:96 bug fix): https://www.mathworks.com/matlabcentral/fileexchange/28922-list-queue-stack
addpath(genpath('toolbox_graph'));
addpath(genpath('datastructure'));

close all;
name = 'clean_mesh';
global PLOT PLOT_FIG plot_mesh_original plot_mesh_curvature plot_mesh_raw plot_mesh_clean plot_mesh_optimized plot_particle plot_volume;
PLOT = true;
PLOT_FIG = 1;
    plot_mesh_original = 0;
    plot_mesh_curvature = 0;
    plot_mesh_raw = 0;
    plot_mesh_clean = 0;
    plot_mesh_optimized = 0;
    plot_particle = 1;
    plot_volume = 0;

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

% Clean workspace
global_list = who('global');
clearvars('-except', global_list{:});

% Display original mesh
if PLOT && plot_mesh_original
    face_colors = 0.3 * ones(nface, 1); % 0.3-red
    plot_face_color('Original Mesh', 0);
end

%% Compute face angle
tic
global face_angles;
K = 3; % order of neighbor search ring
face_angles = compute_face_angle(K);
% larger = face_angles > t;
% sum(larger(:))/length(face_angles) 
% jiayi *1.1 is wanmei number! can be used to compute the 0.9 constant
%face_angles = (face_angles - min(face_angles)) ./ (max(face_angles) - min(face_angles));
fprintf('Compute face angle: %f seconds\n', toc);

% Display curvature mesh
if PLOT && plot_mesh_curvature
    face_colors = face_angles;
    plot_face_color('Curvature Mesh', 0);
end

%% Breadth-First Search Traversal (Curvature criterion)
% Dynamic thresholding
threshold = 0.7; % curvature for convexity/concavity (adjustable, and since both vectors are normalized, this threshold value is actually a critical angle that you can specify)
increment = 0.01;
INCREMENT = true;
if INCREMENT
    fprintf('Incremental mode...\n');
else
    fprintf('Non-Incremental mode...\n');
end
max_faces = 0;
max_objects = 0;
optimal_threshold = threshold; % cache optimal results
optimal_object_no = 0;
optimal_object_set = logical(zeros(nface, []));
while true
    face_remain = logical(ones(nface, 1)); % 1-remain faces
    face_segmented = logical(zeros(nface, 1)); % 1-segmented faces
    object_no = 0; % total No. of particles segmented
    object_set = logical(zeros(nface, [])); % object_set(:,i) is all faces that belong to ith particle
    
    % Thresholds and Controls
    neighboring = 2; % 1-face neighboring (edge adjacency), 2-vertex neighboring (vertex adjacency)
    placeholder = 0; % BFS placeholder (a face index can never be 0 since Matlab is 1-based)
    failure = 0;     % failure times (when it keeps failing, meaning there is no more particles in the remaining mesh, exit the loop)
    
    tic
    % Iteratively segment particles from mesh
    % for x = 1 : 1 % Debug mode (one particle at a time)
    %   seed = XXX; % manual seeding
    while true % Release mode (all particles at a time)

    % Stop sign (when remaining faces is less than certain percentage of total faces)    
    if sum(face_remain(:)) < 0.1 * nface || failure > 20
        break;
    end

    % Seeding
    % Option 1: random seeding
    % rng('shuffle'); % random seeding
    face_remain_idx = find(face_remain == 1);
    seed = face_remain_idx(randi(length(face_remain_idx))); % randomly start from one of the remain faces
    % Option 2: start from face with large convexity (thus more likely on
    % object surface), Caveat: some incomplete particle may not have convex
    % faces. It's evil (missing particles), and it's good (pre-filter
    % incomplete particles for us, and more stable)
    [~, idx] = maxk(face_angles(face_remain_idx), round(0.1 * length(face_remain_idx)));
    seed = face_remain_idx(idx(randi(length(idx))));

    % Initialization
    % Label all faces into one of the following categories in array 'state':
    % Cat 0: Unvisited face (black-colored)
    % Cat 1: Object face    (red-colored)
    % Cat 2: Boundary face  (white-colored)
    state = zeros(nface, 1); % array for recording face status: 0-unvisited face(black), 1-object face(red), 2-boundary face(white)
    state(face_segmented) = 1; % label already segmented faces as boundary
    q = CQueue();            % queue for BFS, storing the face index, i.e. ith face. Any face that entered into this queue is labelled 'object face'.
    q.push(seed); q.push(placeholder);
    state(seed) = 1; % label starting face as 'object'
    update_times = 0; % record the times of a ring BFS

    % BFS
    % tic
    while q.isempty() ~= 1
        % Pop a face from the queue
        curr_id = q.pop();

        % If a placeholder is met, increment the BFS loop count
        if curr_id == placeholder
            update_times = update_times + 1;
            if q.isempty() ~= 1
                q.push(placeholder); % if there are unvisited faces in queue, push the placeholder
            else
                break; % otherwise, this placeholder is the last element, exit (otherwise you will infinite loop)
            end
            continue; % for placeholder step, just skip to the next loop after updating
        end

        % Label neighboring faces into object face (1) OR boundary face (2) 
        if neighboring == 1
            neighbors = face_rings{curr_id};
        else
            neighbors = vertex_rings{curr_id};
        end
        for f = 1 : length(neighbors)
            adj_id = neighbors(f);
            if state(adj_id) == 0 % skip all visited faces
                if face_angles(adj_id) >= threshold
                    % object face (push to queue)
                    state(adj_id) = 1;
                    q.push(adj_id); 
                else
                    % boundary face (do not push to queue, BFS ends at this face)
                    state(adj_id) = 2;
                    % label its adjacency as well (stricter criterion, may or may not need this)
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

    % fprintf('BFS: %f seconds\n', toc);
    % fprintf('BFS completes at %d loops\n', update_times);

    % Display raw mesh
    if PLOT && plot_mesh_raw
        % Per-face coloring: Non-visited face-0 (black), object vertex-0.3 (red), boundary vertex-1.0 (white)
        face_colors = zeros(nface, 1);
        objects = state == 1; face_colors(objects) = 0.3;
        boundaries = state == 2; face_colors(boundaries) = 1.0;
        plot_face_color('Raw Mesh', 1);
        % Label the seed location by plotting normal vector
        hold on;
        id = seed;
        quiver3(face_centers(1,id), face_centers(2,id), face_centers(3,id), face_normals(1,id)/4, face_normals(2,id)/4, face_normals(3,id)/4, 'r');
    end

    % Mesh cleaning by Detecting connected components in unvisited faces (BFS)
    % tic
    state2 = state ~= 0 & face_segmented == 0; % state ~= 0 collects all current visited faces, face_segmented == 0 collects all unvisited faces from last loop, so the intersection is the objects faces identified in the current loop
    % previous BUG here (fixed): if I remove the face_segmented part here,
    % the complement step below may result in an overlapped portion
    % including the already-segmented part. Should take the intersect when
    % this clean mesh step is completed.
    state_new = state ~= 0 & face_segmented == 0; % the newly traversed faces (before cleaning), this is what we should apply the 0.01 criterion (tiny regions) to
    if sum(state_new(:)) < 0.01 * nface
        % if only tiny regions are visited (most likely means an useless BFS),
        % ignore them and skip the following mesh cleaning process
        failure = failure + 1;
        % fprintf('BFS skip\n');
        continue;
    end

    connect = 0; % No. of connected components
    component = {};
    while true
        unvisited = find(state2 == 0); % index of all unvisited faces
        if (length(unvisited) == 0) 
            break; % if no more unvisited faces, exit
        end
        starter = unvisited(round(max(length(unvisited))/2+min(length(unvisited))/2));
        q.empty();
        q.push(starter);
        state2(starter) = 1;
        face_set = zeros(nface, 1);
        face_set(starter) = 1; % don't forget!
        count = 0;
        while q.isempty() ~= 1
            curr_id = q.pop();
            count = count + 1;
            if neighboring == 1
                neighbors = face_rings{curr_id};
            else
                neighbors = vertex_rings{curr_id};
            end
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
        % fprintf('Connected Component %d: %d faces, starter %d\n', connect, count, starter);
    end
    % fprintf('Total connected components: %d\n', connect); 
    % There are usually one big component (the unvisited west), several medium
    % components (they're unvisited regions on the object surface) and many tiny
    % components (with < 10 faces, they're small unvisited faces that are
    % trapped by a surrounding of boundary faces). What we need is the
    % complement of the big component.

    % Find the connected component that contains the max No. of faces
    max_set = 1;
    max_face = 0;
    for c = 1 : connect
        if length(component{c}) > max_face
            max_set = c;
            max_face = length(component{c});
        end
    end

    % Clean the mesh by taking the complement of the largest component
    objects = ~component{max_set} & face_segmented == 0; % bug fix
    
    % Extract boundary faces from the cleaned mesh
    boundaries = boundary_extract(objects);

    % fprintf('Mesh cleaning: %f seconds\n', toc);

    % Display cleaned mesh
    if PLOT && plot_mesh_clean
        face_colors = zeros(nface, 1);
        face_colors(objects) = 0.3;
        face_colors(boundaries) = 1.0;
        plot_face_color('Cleaned Mesh', 1);
    end

    % Record successfully segmented particle faces
    face_remain(objects == 1) = 0;
    face_segmented(objects == 1) = 1;
    object_no = object_no + 1;
    object_set(:, object_no) = objects;
    % fprintf('Segment out particle No.%d, %d faces\n', object_no, sum(objects(:)));
    failure = 0; % if a particle is successfully segmented, reset 'failure' count

    end

    % fprintf('Segmentation: %f seconds\n', toc);

    if INCREMENT == false
        break; % non incremental mode
    end

    if sum(object_set(:)) > max_faces
        max_faces = sum(object_set(:));
    end
    if object_no > max_objects
        max_objects = object_no;
        optimal_threshold = threshold;
        optimal_object_no = object_no;
        optimal_object_set = object_set;
    end
    
    fprintf('Total objects: %d, Total faces: %d, Threshold: %f\n', object_no, sum(object_set(:)), threshold);
    
    if sum(object_set(:)) < 0.9 * max_faces
        break;
    end
    
    threshold = threshold + increment;

end

if INCREMENT
    fprintf('Best threshold found at %f.Task completed!\n', optimal_threshold);
    object_no = optimal_object_no;
    object_set = optimal_object_set;
else
    fprintf('Fix threshold at %f. Task completed!\n', threshold);
end

%% Display segmented particle(s)
% Particles on different figures
% if PLOT && plot_particle 
%     for i = 1 : object_no
%         face_colors = zeros(nface, 1);
%         face_colors(object_set(:, i)) = 0.3;
%         plot_face_color(strcat('Particle ', num2str(i)), 1);
%     end
% end

% Particles on one figure
if PLOT && plot_particle
    face_colors = zeros(nface, 1);
    for i = 1 : object_no
        face_colors(object_set(:, i)) = i / object_no; % better color option here using linspecer.m
    end
    plot_face_color(strcat('Particle ', num2str(i)), 0);
    colormap jet;
end

%% Compute volume(s)
volumes_raw = zeros(object_no, 1);
particle_points{object_no} = [];
% Particles on different figures
% for i = 1 : object_no
%     object_faces = faces(:, object_set(:,i));
%     object_vertices = vertex(:, unique(object_faces(:)))';
%     particle_points{i} = object_vertices;
%     [B, volumes_raw(i)] = boundary(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
%     % plot
%     if PLOT && plot_volume
%         figure(PLOT_FIG); PLOT_FIG = PLOT_FIG + 1;
%         title(strcat('Particle ', num2str(i)));
%         subplot(1,2,1);
%         scatter3(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
%         axis equal off;
% 
%         subplot(1,2,2);
%         trisurf(B,object_vertices(:,1),object_vertices(:,2),object_vertices(:,3),'Facecolor','red','FaceAlpha',0.1)
%         axis equal off;
%     end
% end

% Particles subplotted on one figure
% if PLOT && plot_volume
%     figure(PLOT_FIG); PLOT_FIG = PLOT_FIG + 1;
%     % Arrange subplots
% %     fig_row = 2 * round(sqrt(object_no));
% %     fig_col = 2 * ceil(object_no/(fig_row/2));
%     fig_row = ceil(object_no/2);
%     fig_col = 2 * 2;
% end
% for i = 1 : object_no
%     object_faces = faces(:, object_set(:,i));
%     object_vertices = vertex(:, unique(object_faces(:)))';
%     particle_points{i} = object_vertices;
%     
%     % Calculate volume encompassed by a set of points
%     % volume_raw(i) = volumeFromPoints(object_vertices);
%     [B, volumes_raw(i)] = boundary(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
%     
%     % plot
%     if PLOT && plot_volume
%         subplot(fig_row,fig_col,2*i-1);
%         scatter3(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
%         axis equal off;
% 
%         subplot(fig_row,fig_col,2*i);
%         trisurf(B,object_vertices(:,1),object_vertices(:,2),object_vertices(:,3),'Facecolor','red','FaceAlpha',0.1)
%         axis equal off;
%     end
% end

% Cache
% save('particles.mat', 'particle_points');

%% Optimize boundary (shortest path)
STOP = false;
if STOP
tic
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

fprintf('Optimize boundary: %f seconds\n', toc);

% Display optimized mesh
if PLOT && plot_mesh_optimized
    face_colors = zeros(nface, 1);
    face_colors(objects) = 0.3;
    face_colors(boundaries) = 1.0;
    plot_face_color('Optimized Mesh', 1);
end
end
