function boundary_faces = boundary_extract(objects)
global faces;
% Extract the boundary faces of a given set of faces
% Input:
%   objects: 1 x F logical array of the target faces among all faces
% Output:
%   boundary_faces: an array of indices that forms a closed boundary ring

object_faces = find(objects == 1);
object_face_rings = compute_face_ring(faces(:, object_faces)); % compute face ring on the subset object faces

boundary_face = zeros(length(object_faces), 1);
for f = 1 : length(object_faces)
    if length(object_face_rings{f}) < 3
        boundary_face(f) = 1;
    end
end
boundary_faces = object_faces(boundary_face == 1);

% old
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



