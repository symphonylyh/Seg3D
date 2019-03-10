function [object_no,object_set] = BFS_regional(threshold)
% Regional version of Breadth First Search
%   1. "Regional" means BFS applys to connected/disconnected mesh constructed 
%   from 'vertex' & 'faces'
%   2. This version is slower than BFS_universal() due to index mapping
%   3. This version is versatile in that it can apply to both types of mesh
global PLOT plot_mesh_raw plot_mesh_clean;
global vertex faces nvertex nface face_rings vertex_rings face_normals face_centers face_colors face_angles;
% Double Map for the face subset indices
global face_subset_idx; % linear indices (0,1,...) to global face indices (5876,1293,...)
global M; % global face indices (5876,1293,...) to linear indices (0,1,...)
% Map convertion (everything pushed to queue is the global face index! everyting queried in the global variables is the mapped index!)

% Dynamic thresholding
% threshold = 0.7; % curvature for convexity/concavity (adjustable, and since both vectors are normalized, this threshold value is actually a critical angle that you can specify)
increment = 0.01;
INCREMENT = false;
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
    state(seed) = 1; % label starting face as 'object'

    % Map convertion (everything pushed to queue is the global face index! everyting queried in the global variables is the mapped index!)
    seed = face_subset_idx(seed);
    q.push(seed); q.push(placeholder);
    
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
            neighbors = face_rings{M(curr_id)};
        else
            neighbors = vertex_rings{M(curr_id)};
        end
        for f = 1 : length(neighbors)
            adj_id = neighbors(f);
            if state(M(adj_id)) == 0 % skip all visited faces
                if face_angles(M(adj_id)) >= threshold
                    % object face (push to queue)
                    state(M(adj_id)) = 1;
                    q.push(adj_id); 
                else
                    % boundary face (do not push to queue, BFS ends at this face)
                    state(M(adj_id)) = 2;
                    % label its adjacency as well (stricter criterion, may or may not need this)
                    if neighboring == 1
                        boundary_neighbors = face_rings{M(adj_id)};
                    else
                        boundary_neighbors = vertex_rings{M(adj_id)};
                    end
                    for b = 1 : length(boundary_neighbors)
                        boundary_id = boundary_neighbors(b);
                        if state(M(boundary_id)) == 0
                            state(M(boundary_id)) = 2;
                        end
                    end
                end
            end
        end
    end

    % fprintf('BFS: %f seconds\n', toc);
    % fprintf('BFS completes at %d loops\n', update_times);

    % Display raw mesh (this will accumulatively display the newly identified faces along with already segmented parts)
    if PLOT && plot_mesh_raw
        % Per-face coloring: Non-visited face-0 (black), object vertex-0.3 (red), boundary vertex-1.0 (white)
        face_colors = zeros(nface, 1);
        objects = state == 1; face_colors(objects) = 0.3;
        boundaries = state == 2; face_colors(boundaries) = 1.0;
        plot_face_color('Raw Mesh', 1);
        % Label the current seed location by plotting normal vector
        hold on;
        id = M(seed); % since seed is converted by face_subset_idx, here we should convert back to query face_centers and face_normals
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
    starters = {}; % for plot
    while true
        unvisited = find(state2 == 0); % index of all unvisited faces ('unvisited' in the context is the complement of the newly-identified faces i.e. complement(state2))
        if (length(unvisited) == 0) 
            break; % if no more unvisited faces, exit
        end
        starter = unvisited(round(max(length(unvisited))/2+min(length(unvisited))/2));
        q.empty();
        state2(starter) = 1;
        face_set = zeros(nface, 1);
        face_set(starter) = 1; % don't forget!
        
        % Map convertion (everything on queue is the global face index!)
        starter = face_subset_idx(starter);
        q.push(starter);

        count = 0;
        while q.isempty() ~= 1
            curr_id = q.pop();
            count = count + 1;
            if neighboring == 1
                neighbors = face_rings{M(curr_id)};
            else
                neighbors = vertex_rings{M(curr_id)};
            end
            for f = 1 : length(neighbors)
                adj_id = neighbors(f);
                if state2(M(adj_id)) == 0
                    state2(M(adj_id)) = 1;
                    face_set(M(adj_id)) = 1;
                    q.push(adj_id);
                end
            end
        end
        connect = connect + 1;
        component{connect} = face_set; % component{i} is the logical face index array of the ith component
        starters{connect} = starter; % for plot
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
        if sum(component{c}) > max_face
            max_set = c;
            max_face = sum(component{c});
            % Bug fix March 9:
            % if length(component{c}) > max_face
            % max_face = length(component{c}); 
            % 'length' is always = nface...should both be 'sum'!
        end
    end

    % Clean the mesh by taking the complement of the largest component
    if connect == 0 % special case: if there is nothing left (the whole mesh is traversed), component{} is empty, connect = 0, we should skip
        objects = face_segmented == 0; % bug fix March_6
    else
        objects = ~component{max_set} & face_segmented == 0; % bug fix
    end
    
    % Extract boundary faces from the cleaned mesh
    boundaries = boundary_extract(objects);

    % fprintf('Mesh cleaning: %f seconds\n', toc);

    % Display cleaned mesh (this will only display the newly identified faces in current iteration)
    if PLOT && plot_mesh_clean
        face_colors = zeros(nface, 1);
        face_colors(objects) = 0.3;
        face_colors(boundaries) = 1.0;
        plot_face_color('Cleaned Mesh', 1);
        % Label the current starter location by plotting normal vector
        hold on;
        id = M(starters{max_set}); % since starter is converted by face_subset_idx, here we should convert back to query face_centers and face_normals
        quiver3(face_centers(1,id), face_centers(2,id), face_centers(3,id), face_normals(1,id)/4, face_normals(2,id)/4, face_normals(3,id)/4, 'b');
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
        break; % non incremental mode stops here
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

end

