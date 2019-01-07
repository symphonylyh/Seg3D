function v = volumeFromPoints(points)
% Calculate the volume encompassed by a point cloud
% Input:
%   points: N x 3 matrix, points(i, :) is the 3D coordinates of ith point
% Return:
%   volume: volume of the entity encompassed by the set of points

% Two methods are optional:
% 1: 3D Denaulay triangulation + Convex hull (convex hull is not we want)
% 2: 3D Boundary face (not necessarily convex, exactly what we want)
% For each method, I list two ways of calculation, one is how the built-in
% MATLAB function does for us, the other is how to manually calculate it.
% For method 1: Built-in convexHull() OR summing up volume of tetrahedrons.
% For method 2: Built-in boundary()'s alphaShape OR Gauss's divergence theorem (for closed surface mesh ONLY).
method = 1;

% Method 1
if method == 1
% Built-in
    DT = delaunayTriangulation(points);
    [C, v] = convexHull(DT); 
    % C: connectivity list (vertex indices) of the triangles on the convex hull
    % v: volume of the convex hull
    trisurf(C,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),'Facecolor','red','FaceAlpha',0.1);
    axis equal;
% Manual calculation (very slow)
%     vol = 0; area = 0;
%     for t = 1 : length(DT.ConnectivityList) % each one in the connectivity list is a tetrahedron
%         indices = DT.ConnectivityList(t, :);
%         V1 = DT.Points(indices(1), :);
%         V2 = DT.Points(indices(2), :);
%         V3 = DT.Points(indices(3), :);
%         V4 = DT.Points(indices(4), :);
%         % Calculate volume and area of a tetrahedron
%         V = 1/6 * abs(dot(cross(V2-V1,V3-V1),V4-V1));
%         A123 = 1/2 * norm(cross(V2-V1,V3-V1));
%         A124 = 1/2 * norm(cross(V2-V1,V4-V1));
%         A134 = 1/2 * norm(cross(V3-V1,V4-V1));
%         A234 = 1/2 * norm(cross(V3-V2,V4-V2));
%         vol = vol + V;
%         area = area + A123 + A124 + A134 + A234;
%     end
% v == vol can be verified
end

% Method 2
if method == 2
% Built-in
    [B, v] = boundary(points(:,1), points(:,2), points(:,3));
    trisurf(B,points(:,1),points(:,2),points(:,3),'Facecolor','red','FaceAlpha',0.1)
% Manual calculation (boundary face mesh returned from boundary() is ensured close)
    [vol, area] = triangulationVolume(B,points(:,1),points(:,2),points(:,3)); % triangulationVolume is a function based on Gauss's divergence theorem
% v == vol can be verified
end

end

% [Ref] Volume calculation:
% https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
% [Ref] Triangulation methods:
% http://www.csie.ntnu.edu.tw/~u91029/Triangulation.html
