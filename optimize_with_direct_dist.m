function [new_boundary] = optimize_with_direct_dist(E_local, E_global, V_all, boundary_points, num_change) 
% count: how many local paths have been changed
    start_point = boundary_points;
    [local_dist, ~] = dijkstra(V_all,E_local, start_point, boundary_points);
%    [c_g_all_all, p_g_all_all] = dijkstra(V_all, E_global, start_point, boundary_points);
    global_dist = zeros(size(boundary_points, 2), size(boundary_points, 2));
    for i = 1: size(boundary_points, 2)
        for j = 1: size(boundary_points, 2)
            global_dist(i, j) = norm(V_all(boundary_points(i),:) - V_all(boundary_points(j),:));
        end
    end
    % Find the largest distance change with respected to local distance
    A = (local_dist - global_dist)./local_dist;
    [maximum, ind] = maxk(A(:), 2*num_change);
    new_boundary = boundary_points;
    for i = 1: num_change
        [x,y] = ind2sub(size(A),ind(2*i));
        start = boundary_points(x);
        last = boundary_points(y);
        % Change the path
        [c_g, p_g] = dijkstra(V_all, E_global, start, last);
        % Add the path
        new_boundary = union(new_boundary, p_g);    
    end  
   
end