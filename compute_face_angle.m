function face_angles = compute_face_angle(K)
global face_normals face_centers vertex_rings;
% Compute face angle (curvature indicator) for each face. This can be
% integrated into pre-compute mesh step.
% Input:
%   K: order of neighbor-search ring to compute the average face angle (at least 1)
% Output: 
%   face_angles: 1 x F array, the curvature measure at each face 

% Method 1: average up to K rings of neighbors 
face_angles = zeros(length(face_centers), 1); 
for i = 1: length(face_centers)
    curr_center = face_centers(:, i);
    curr_normal = - face_normals(:, i);
    
    % Get up to K rings of neighbors
    neighbor_list = i;
    for k = 1 : K
        neighbor_list = [neighbor_list vertex_rings{neighbor_list}];
    end
    neighbor_list = unique(neighbor_list);
    neighbors = neighbor_list(neighbor_list ~= i);

    % Compute and average neighbor curvature
    % new version (accelerated)
    adj_center = face_centers(:, neighbors);
    adj_normal = - face_normals(:, neighbors);
    local_center_diff = adj_center - curr_center;
    local_center_diff = local_center_diff ./ vecnorm(local_center_diff); % normalize vector along dimension = 1
    angle_set = curr_normal' * (local_center_diff + adj_normal);
    % old version (for-loop)
%     angle_set = zeros(length(neighbors), 1);
%     for j = 1: length(neighbors)
%         adj_center = face_centers(:, neighbors(j));
%         adj_normal = - face_normals(:, neighbors(j));
%         local_center_diff = adj_center - curr_center;
%         local_center_diff = local_center_diff / norm(local_center_diff); % normalize
%         angle_set(j) = dot(curr_normal, local_center_diff + adj_normal);
%     end
    
    % Average the angle_set and store it in face_angles
    face_angles(i) = mean(angle_set);
end

% Method 2 (old): compute 1st-order curvature at each face, and compute curvature based on top of 1st neighbors' 1st-order curvature  
% face_angles = zeros(1, length(face_centers)); 
% % Start algorithm
% for i = 1: length(face_centers)
%     curr_center = face_centers(:, i);
%     curr_normal = - face_normals(:, i);
%     % Loop over the neighbors of current face
%     neighbors = vertex_rings{i};
%     angle_set = [];
%     for j = 1: length(neighbors)
%         adj_center = face_centers(:, neighbors(j));
%         adj_normal = - face_normals(:, neighbors(j));
%         local_center_diff = adj_center - curr_center;
%         local_center_diff = local_center_diff / norm(local_center_diff); % normalize
%         angle = dot(curr_normal, local_center_diff + adj_normal);
%         angle_set = [angle_set angle];
%     end
% 
%     % Average the angle_set and store it in face_angles
%     face_angles(i) = mean(angle_set);
% end
% 
% % Smooth
% face_angles_temp = zeros(size(face_angles));
% for i = 1: length(face_centers)
%     curr_angle = face_angles(i);
%     neighbors = vertex_rings{i};
%     neighbor_angle = [];
%     for j = 1: length(neighbors)
%         neighbor_angle = [neighbor_angle face_angles(neighbors(j))];        
%     end
% 
%     % Sum up & Averaging
%     face_angles_temp(i) = curr_angle * 0.8 + sum(neighbor_angle)/length(neighbor_angle) * 0.2;
% end
% face_angles = face_angles_temp;

end