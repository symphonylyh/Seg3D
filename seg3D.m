%% Control panel
% Libraries related:
% Graph & Mesh: https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph 
% Queue/Stack (Cqueue.m:96 bug fix): https://www.mathworks.com/matlabcentral/fileexchange/28922-list-queue-stack
addpath(genpath('toolbox_graph'));
addpath(genpath('datastructure'));

close all;
% name = 'clean_mesh';
name = 'mesh/multiple_02';
%name = 'mesh/01';
global PLOT PLOT_FIG plot_mesh_original plot_mesh_curvature plot_mesh_region plot_mesh_raw plot_mesh_clean plot_mesh_optimized plot_particle plot_volume;
PLOT = true;
PLOT_FIG = 1;
    plot_mesh_original = 0;
    plot_mesh_curvature = 0;
    plot_mesh_region = 1;
    plot_mesh_raw = 0;
    plot_mesh_clean = 0;
    plot_mesh_optimized = 0;
    plot_particle = 1;
    plot_volume = 1;

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
    face_colors = ones(nface, 1);
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
face_angles = (face_angles - min(face_angles)) ./ (max(face_angles) - min(face_angles));
fprintf('Compute face angle: %f seconds\n', toc);

% Display curvature mesh
if PLOT && plot_mesh_curvature
    face_colors = face_angles;
    plot_face_color('Curvature Mesh', 0);
end

%% Region separation
tic

% Backup all faces parameters
faces_all = faces;
nface_all = nface;
face_angles_all = face_angles;
face_centers_all = face_centers;
face_normals_all = face_normals;
face_rings_all = face_rings;
vertex_rings_all = vertex_rings;

% Separate unconnected regions by different labels (1 to N)
region_visited = zeros(nface_all, 1);
region_label = 0; 

neighboring = 2;
q = CQueue();
while true
    region_unvisited = find(region_visited == 0);
    if length(region_unvisited) == 0
        break;
    end
    region_label = region_label + 1;
    start_id = region_unvisited(1);
    q.push(start_id); region_visited(start_id) = region_label;
    % BFS
    while q.isempty() ~= 1
        curr_id = q.pop();
        if neighboring == 1
            neighbors = face_rings{curr_id};
        else
            neighbors = vertex_rings{curr_id};
        end
        for f = 1 : length(neighbors)
            adj_id = neighbors(f);
            if region_visited(adj_id) == 0
                region_visited(adj_id) = region_label;
                q.push(adj_id);
            end
        end
    end
       
end

fprintf('Region separation: %f seconds\n', toc);

% Display region separated mesh
if PLOT && plot_mesh_region
    face_colors = zeros(nface, 1);
    region_colors = linspace(0, 1, region_label);
    for region_no = 1 : region_label
        face_colors(region_visited == region_no) = region_colors(region_no);
    end
    plot_face_color('Regional Mesh', 0);
end

%% Breadth-First Search Traversal (Curvature criterion) on each region
if region_label == 1
    % Single region BFS
    [object_no,object_set] = BFS_universal();
else
    % Multi-region BFS
    global face_subset_idx; % face subset indices
    global M; % index map
    object_no = 0; % total No. of particles segmented
    object_set = logical(zeros(nface_all, [])); % object_set(:,i) is all faces that belong to ith particle
    for region_no = 1 : region_label % Max No. of region label = total No. of unconnected regions
        tic
        % Prepare subset region info
        face_subset = region_visited == region_no;
        face_subset_idx = find(face_subset == 1);
        faces = faces_all(:, face_subset);
        nface = sum(face_subset);
        face_angles = face_angles_all(face_subset);
        face_centers = face_centers_all(:, face_subset);
        face_normals = face_normals_all(:, face_subset);
        face_rings = face_rings_all(face_subset);
        vertex_rings = vertex_rings_all(face_subset);
        M = zeros(nface_all, 1);
        M(face_subset) = 1:nface;
        % Regional BFS
        [sub_object_no,sub_object_set] = BFS_regional();
        % Convert subset indices to global face set indices
        for obj = 1 : sub_object_no
            global_set = logical(zeros(nface_all, 1));
            global_set(face_subset_idx(sub_object_set(:,obj))) = 1;
            object_no = object_no + 1;
            object_set(:, object_no) = global_set;
        end
        fprintf('Region %d completed, %d particles: %f seconds\n', region_no, sub_object_no, toc);
    end
end

% Restore all faces parameters
faces = faces_all;
nface = nface_all;
face_angles = face_angles_all;
face_centers = face_centers_all;
face_normals = face_normals_all;
face_rings = face_rings_all;
vertex_rings = vertex_rings_all;

%% Display segmented particle(s)
% Particles on different figures
% if PLOT && plot_particle 
%     for i = 1 : object_no
%         % All particles in figure, highlight one
%         % face_colors = zeros(nface, 1);
%         % face_colors(object_set(:, i)) = 0.3;
%         % plot_face_color(strcat('Particle ', num2str(i)), 1);
%         
%         % Single particle in figure
%         faces_all = faces;
%         faces = faces(:, object_set(:,i));
%         face_colors = 0.3 * ones(length(faces), 1);
%         plot_face_color(strcat('Particle ', num2str(i)), 1);
%         faces = faces_all;
%     end
% end

% Particles on one figure
if PLOT && plot_particle
    face_colors = zeros(nface, 1);
    for i = 1 : object_no
        face_colors(object_set(:, i)) = i / object_no; % better color option here using linspecer.m
    end
    plot_face_color('Segmented Particles', 0);
    colormap jet;
end
    
%% Compute volume(s)
volumes_raw = zeros(object_no, 1);
particle_points{object_no} = [];
particle_faces{object_no} = [];
% Particles on different figures
for i = 1 : object_no
    object_faces = faces(:, object_set(:,i));
    object_vertices = vertex(:, unique(object_faces(:)))';
    particle_points{i} = object_vertices;
    particle_faces{i} = object_faces;
    [B, volumes_raw(i)] = boundary(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3), 0.1);
    % plot
    if PLOT && plot_volume
        figure(PLOT_FIG); PLOT_FIG = PLOT_FIG + 1;
        title(strcat('Particle ', num2str(i)));
        subplot(1,2,1);
        scatter3(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
        axis equal off;

        subplot(1,2,2);
        trisurf(B,object_vertices(:,1),object_vertices(:,2),object_vertices(:,3),'Facecolor','red','FaceAlpha',0.1)
        axis equal off;
    end
end

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
%     particle_faces{i} = object_faces;
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
% save('particles.mat', 'particle_points', 'particle_faces', 'vertex', 'faces');

test = false;
if test
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

end