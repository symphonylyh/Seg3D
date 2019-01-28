function [component_no,component_set] = componentCount()
global vertex faces nvertex nface face_rings vertex_rings;
% Check the connected components in a graph/mesh. This is a sanity check
% step: if there is only one graph/mesh including all particles, the BFS is the
% same as before; if the mesh is comprised of multiple components, each
% component needs a separate BFS (by updating the global 'faces').

component_no = 0; % total No. of connected components
component_set = logical(zeros(nface, [])); % component_set(:, i) is all faces that belong to ith component

visited = zeros(nface, 1);
q = CQueue();
while true
    unvisited = find(visited == 0); % index of all unvisited faces
    if (length(unvisited) == 0) 
        break; % if no more unvisited faces, exit
    end
    starter = unvisited(rand(length(unvisited)));
    q.empty();
    q.push(starter);
    visited(starter) = 1;
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
            if visited(adj_id) == 0
                visited(adj_id) = 1;
                face_set(adj_id) = 1;
                q.push(adj_id);
            end
        end
    end
    component_no = component_no + 1;
    component_set(:, component_no) = objects;
end

end

