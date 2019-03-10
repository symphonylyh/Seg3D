function [scale, centroid, dist_map, dist_mask] = shapeDeflate(cloudVertex, objectFace)
% Deflate a 3D shape onto normalized 2D distance map.
% Input:
%   cloudVertex: 3 x V matrix, vertex coordinates (for simplification, this is the vertex of the whole point cloud)
%   objectFace:  3 x F matrix, indexed-face (for simplification, this is the indexed-face among the whole point cloud that belongs to one particle)
% Output:
%   scale: scalar, used for later recovery the real distance = dist * scale, picked as the maximum value of distance
%   centroid: 3 x 1 vector, used for later recovery the real coordinates
%   dist_map: N_theta x N_phi 2D image, with pixel values being the normalized distance from centroid
%   dist_mask: N_theta x N_phi mask labelling the missing region

%% Libraries related
% 1. Centroid calculation (cloudCentroid.m, ignore the error checking): https://www.mathworks.com/matlabcentral/fileexchange/8514-centroid-of-a-convex-n-dimensional-polyhedron
% 2. Uniform sampling points on sphere surface: 
% sphereSampling: https://www.mathworks.com/matlabcentral/fileexchange/37004-suite-of-functions-to-perform-uniform-sampling-of-a-sphere
% OR
% spheretri: https://www.mathworks.com/matlabcentral/fileexchange/58453-spheretri
addpath(genpath('sphereSampling'));
addpath(genpath('spheretri'));
% 3. Ray-Mesh collision (first run src_matlab/mexall.m): https://www.mathworks.com/matlabcentral/fileexchange/41504-ray-casting-for-deformable-triangular-3d-meshes
addpath(genpath('opcodemesh'));
% 4. Other libraries that might be useful:
% Geom3D: https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d
% Mesh-Mesh intersection: https://www.mathworks.com/matlabcentral/fileexchange/49160-fast-mesh-mesh-intersection-using-ray-tri-intersection-with-octree-spatial-partitioning
% Simpler Ray-Mesh intersection using Möller and Trumbore (1997): https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection
% Lecture notes: http://www.cs.cornell.edu/courses/cs417/2003sp/Lectures/Lecture33/33rayintersection.pdf

% Key idea
% Based on the point cloud of a segmented particle:
% Step 1: Estimate the centroid of the point cloud
    % Old attempt: average the point coordinates. This is flawed for incomplete
    % particles b/c there are missing portions of the particle that cannot
    % contribute to the centroid
    % New attempt: use centroid of the convex hull. This is better and may be
    % the best-effort solution. The convex hull reconstruction to some extent
    % compensate the effect of missing portions (although not perfect)
% Step 2: Evaluate all-around shape faces using ray-tracing
    % Step 2.1: Generate rays that cover the all-around space
    % Staring from the centroid, shoot rays in all directions. 
    % Old attempt: convert to phi-theta spherical coordinates. This is not
    % uniformly distributed on a sphere.
    % New attempt: generate uniform distributed points on a unit sphere by
    % subdivision
    % Step 2.2: Detect ray-mesh collision
    % Check how many rays will intersect with the object's mesh faces, and the 
    % percentage of intersected directions indicates the particle shape completeness
% Step 3: Deflate the 3D directions onto 2D spherical coordinates, and
% set the pixel value at the coordinate as normalized face-centroid distance

% -------------------------------------------------------------------------
% Step 1 Compute centroid
% -------------------------------------------------------------------------   
object_vertices = cloudVertex(:, unique(objectFace(:)))';  
% Old attempt
centroid = mean(object_vertices, 1);
% New attempt
% centroid = cloudCentroid(object_vertices); 
% very slow, try: https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d
% try use MATLAB builtin alphaShape centroid functions

% -------------------------------------------------------------------------
% Step 2.1 Initialize directional vectors
% -------------------------------------------------------------------------
nSamples = 1000; % No. of directions to be generated on a unit sphere
[sample_vertex, sample_face] = spheretri(nSamples); % Note: this library will not generate EXACT number of samples
nSamples = length(sample_vertex);
% Visualize
% figure('color','w') 
% fv=struct('faces',sample_face,'vertices',sample_vertex); 
% h=patch(fv); 
% set(h,'EdgeColor','b','FaceColor','w') 
% axis equal	
% set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1]) 
% view(3) 
% grid on 

% -------------------------------------------------------------------------
% Step 2.2 Ray-Mesh Intersection
% -------------------------------------------------------------------------
% Initialize the AABB tree for collision detection (such data structure is flexible, i.e. mesh is deformable)
tree = opcodemesh(cloudVertex,objectFace);
% Check collision
orig = repmat(centroid', [1 length(sample_vertex)]);
dir = sample_vertex';
[hit,d,hit_face,hit_bary,hit_point] = tree.intersect(orig,dir); 
% Input:
%   orig: 3 x N matrix, origins of vector
%   dir:  3 x N matrix, directions of vector
% Output:
%   hit: N x 1 logical array, 0 if no intersection, 1 if intersected
%   d: N x 1 array, distance from origin to intersection point (NaN if the ray misses)
%   hit_face: N x 1 array, index of the intersected face (0 if the ray misses)
%   hit_bary: 2 x N matrix, barycentric coordinates of the intersected
%   points (NaN if the ray misses)
%   hit_point: 3 x N matrix, 3D coordinates of the intersection points (NaN if the ray misses)

% -------------------------------------------------------------------------
% Step 3 Deflate onto 2D distance map
% ------------------------------------------------------------------------- 
% Normalize face distance to 0~1
distance = d / max(d); % Note: distance for no-intersection points will be NaN
% distance = (d - min(d))/(max(d) - min(d));
% distance = abs(d - mean(d)) / mean(d);
scale = max(d);

N_theta = 32; % image height
N_phi = 32; % image width
dist_map = zeros(N_theta, N_phi);
dist_mask = zeros(N_theta, N_phi);
cumulate_count = zeros(N_theta, N_phi);

% Cartesian coordinates (x, y, z) to spherical coordinates (r, phi, theta)
% r = sqrt(x^2 + y^2 + z^2);
% phi = atan2(y / x); range [-pi, pi]
% theta = acos(z / r); range [0, pi], 0-north pole, pi-south pole
% Spherical coordinates (r, phi, theta) to Cartesian coordinates (x, y, z)
% x = r * sin(theta) * cos(phi);
% y = r * sin(theta) * sin(phi);
% z = r * cos(theta);
% In our case, the normal vector are already normalized so r = 1
phi = atan2(sample_vertex(:,2), sample_vertex(:,1)); % azimuthal angle
theta = acos(sample_vertex(:,3)); % polar angle

% Grid bins, something like a histogram voting
d_theta = (pi - 0) / N_theta;
d_phi = (pi - (-pi)) / N_phi;

% First create the mask
for n = 1 : nSamples
    if hit(n) == 0
        i = max(1, ceil( (theta(n) - 0) / d_theta ));
        j = max(1, ceil( (phi(n) - (-pi)) / d_phi ));
        dist_mask(i,j) = 1;
    end
end

% Then fill in the distance value filtered by the above mask
for n = 1 : nSamples
    i = max(1, ceil( (theta(n) - 0) / d_theta ));
    j = max(1, ceil( (phi(n) - (-pi)) / d_phi ));
    if dist_mask(i,j) ~= 1
        dist_map(i,j) = dist_map(i,j) + distance(n);
        cumulate_count(i,j) = cumulate_count(i,j) + 1;
    end
end
cumulate_count(cumulate_count == 0) = 1; % to avoid divide by 0
dist_map = dist_map ./ cumulate_count; % averaging
% dist_map(dist_mask == 1) = -1; % label no-intersection regions as -1

PLOT = false;
if PLOT
% Plot histogram
figure(1)
imagesc(dist_map), colormap gray, axis image;
end

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
visual = false;
if visual
    figure('color','w') 
    % Plot spherical directions
    centroid_dist = object_vertices - centroid;
    centroid_dist = vecnorm(centroid_dist, 2, 2);
    fv=struct('faces',sample_face,'vertices',centroid + 0.3*max(centroid_dist)/1*sample_vertex);
    % per-vertex coloring
    vertex_color = zeros(length(sample_vertex), 1);
    vertex_color(hit) = 1;
    h=patch(fv,'FaceVertexCData',vertex_color,'FaceColor','interp'); 
    % h.LineStyle = 'none'; % don't plot edges
    %set(h,'EdgeColor','b','FaceColor','w') 
    axis equal	
    %set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1]) 
    view(3) 
    grid on 
    % Plot particle point cloud
    hold on
    scatter3(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
end

end

