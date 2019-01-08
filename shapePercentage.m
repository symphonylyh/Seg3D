function percent_shape = shapePercentage(cloudVertex, cloudFace)
% Determine the percentage of shape completion of a point cloud
% Input:
%   cloudVertex: 3 x V matrix, vertex coordinates
%   cloudFace  : 3 x F matrix, indexed-face
% Output:
%   percent_shape: percentage of shape completion

% Libraries related:
% Ray-Mesh collision (first run src_matlab/mexall.m): https://www.mathworks.com/matlabcentral/fileexchange/41504-ray-casting-for-deformable-triangular-3d-meshes
addpath(genpath('opcodemesh'));
% Uniform sampling on sphere surface: 
% sphereSampling: https://www.mathworks.com/matlabcentral/fileexchange/37004-suite-of-functions-to-perform-uniform-sampling-of-a-sphere
% OR
% spheretri: https://www.mathworks.com/matlabcentral/fileexchange/58453-spheretri
addpath(genpath('sphereSampling'));
addpath(genpath('spheretri'));

%% Initialize directional vectors
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

%% Ray-Mesh Intersection
% Initialize the AABB tree (such data structure is flexible, i.e. mesh is deformable)
tree = opcodemesh(cloudVertex,cloudFace);
% Check collision
object_vertices = cloudVertex(:, unique(cloudFace(:)))';
centroid = mean(cloudVertex(:, unique(cloudFace(:))), 2);
orig = repmat(centroid, [1 length(sample_vertex)]);
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
% Plot
figure('color','w') 
fv=struct('faces',sample_face,'vertices',centroid' + 0.03*sample_vertex);
vertex_color = zeros(length(sample_vertex), 1);
vertex_color(hit) = 1;
h=patch(fv,'FaceVertexCData',vertex_color,'FaceColor','interp'); 
%set(h,'EdgeColor','b','FaceColor','w') 
axis equal	
%set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1]) 
view(3) 
grid on 
hold on
scatter3(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
percent_shape = sum(hit(:))/length(hit);
end

