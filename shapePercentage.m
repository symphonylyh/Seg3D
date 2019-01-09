function percent_shape = shapePercentage(cloudVertex, objectFace)
% Determine the percentage of shape completeness of a segmented particle in the point cloud
% Input:
%   cloudVertex: 3 x V matrix, vertex coordinates (for simplification, this is the vertex of the whole point cloud)
%   objectFace:  3 x F matrix, indexed-face (for simplification, this is the indexed-face among the whole point cloud that belongs to one particle)
% Output:
%   percent_shape: percentage of shape completeness

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

% Key idea
% Based on the point cloud of a segmented particle:
% Step 1: Estimate the centroid of the point cloud
    % Old attempt: average the point coordinates. This is flawed for incomplete
    % particles b/c there are missing portions of the particle that cannot
    % contribute to the centroid
    % New attempt: use centroid of the convex hull. This is better and may be
    % the best-effort solution. The convex hull reconstruction to some extent
    % compensate the effect of missing portions (although not perfect)
% Step 2: Evaluate spatial shape completeness using ray-tracing
    % Step 2.1: Generate rays that cover the all-around space
    % Staring from the centroid, shoot rays in all directions. 
    % Old attempt: convert to phi-theta spherical coordinates. This is not
    % uniformly distributed on a sphere.
    % New attempt: generate uniform distributed points on a unit sphere by
    % subdivision
    % Step 2.2: Detect ray-mesh collision
    % Check how many rays will intersect with the object's mesh faces, and the 
    % percentage of intersected directions indicates the particle shape completeness

% -------------------------------------------------------------------------
% Step 1 Compute centroid
% -------------------------------------------------------------------------   
object_vertices = cloudVertex(:, unique(objectFace(:)))';  
% Old attempt
% centroid = mean(object_vertices, 1);
% New attempt
centroid = cloudCentroid(object_vertices); 
% very slow, try: https://www.mathworks.com/matlabcentral/fileexchange/24484-geom3d
% try use MATLAB builtin alphaShape centroid functions

% -------------------------------------------------------------------------
% Step 2.1 Initialize directional vectors
% -------------------------------------------------------------------------
method = 1;
if method == 1
    % Option 1
    nSamples = 1000; % No. of directions to be generated on a unit sphere
    [sample_vertex, sample_face] = spheretri(nSamples); % Note: this library will not generate EXACT number of samples
    % Visualize
    % figure('color','w') 
    % fv=struct('faces',sample_face,'vertices',sample_vertex); 
    % h=patch(fv); 
    % set(h,'EdgeColor','b','FaceColor','w') 
    % axis equal	
    % set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1]) 
    % view(3) 
    % grid on 
else
    % Option 2
    [sample_vertex,sample_face,~,~]=ParticleSampleSphere('N',200); % Note: this library will generate EXACT number, but slower; need subdivision
    fv = struct('faces',sample_face,'vertices',sample_vertex);
    fv = SubdivideSphericalMesh(fv,2);
    % Visualize
    % figure('color','w') 
    % subplot(1,2,1) 
    % fv=struct('faces',sample_face,'vertices',sample_vertex); 
    % h=patch(fv); 
    % set(h,'EdgeColor','b','FaceColor','w') 
    % axis equal	
    % set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1]) 
    % view(3) 
    % grid on 
    % set(get(gca,'Title'),'String','N=200 (base mesh)','FontSize',30)
    % fv_new=SubdivideSphericalMesh(fv,2); 
    % subplot(1,2,2) 
    % h=patch(fv_new); 
    % set(h,'EdgeColor','b','FaceColor','w') 
    % axis equal	
    % set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1]) 
    % view(3) 
    % grid on 
    % set(get(gca,'Title'),'String','N=3170 (after 2 subdivisions)','FontSize',30)
end

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

% Calculate percentage
percent_shape = sum(hit(:))/length(hit);

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

