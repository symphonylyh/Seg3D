function boundary_face_idx = boundary_extract(faces, objects)
% Extract the boundary faces of a given set of faces
% Input:
%   face_rings: connectivity of each face
%   objects: a label array of the subset of target faces among all faces
% Return:
%   boundary: indices of a close boundary ring
object_face_idx = find(objects == 1);
object_faces = faces(:, object_face_idx);
object_face_rings = compute_face_ring(object_faces);

boundary_face = zeros(size(object_faces, 2), 1);
for f = 1 : length(object_face_rings)
    if length(object_face_rings{f}) < 3
        boundary_face(f) = 1;
    end
end
boundary_face_idx = object_face_idx(boundary_face == 1);

% DFS to reach the edge face
% s = CStack();
% starter = randi(size(object_faces, 2)); % random start
% visited = zeros(size(object_faces, 2), 1);
% s.push(starter); visited(starter) = 1;
% found = 0;
% while true
%     curr = s.pop();
%     neighbors = object_face_rings{curr};
%     if length(neighbors) < 3
%         found = curr;
%         break;
%     end
%     for f = 1 : length(neighbors)
%         if visited(neighbors(f)) == 0
%             visited(neighbors(f)) = 1;
%             s.push(neighbors(f));
%         end
%     end
% end
% edge_face_idx = object_face_idx(found)

end



