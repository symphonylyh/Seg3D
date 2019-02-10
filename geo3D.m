close all;

COMPLETNESS = false;

%% Load particle information
load('particles.mat');
% 'faces' and 'vertex' are still global variables

%% Determine shape completeness
if COMPLETNESS
tic

object_no = length(particle_faces);
completion = zeros(object_no, 1);
for i = 1 : length(particle_faces)
    object_faces = particle_faces{i};
    completion(i) = shapePercentage(vertex, object_faces);
end

fprintf('Shape completeness: %f seconds\n', toc);
end

%% Generate incomplete distance map
folder = 'fooling_set';
object_no = length(particle_faces);
scale_all = zeros(length(particle_faces), 1); 
centroid_all = zeros(3, length(particle_faces)); 
for i = 1 : 1 %length(particle_faces)
    object_faces = particle_faces{i};
    % Deflate
    [scale_all(i), centroid_all(:,i), dist_map, dist_mask] = shapeDeflate(vertex, object_faces); 
    % imwrite(dist_map, fullfile(folder, strcat(num2str(i, '%04.f'), '.png')));
    % imwrite(dist_mask, fullfile(folder, strcat(num2str(i, '%04.f'), '_mask.png')));
    % GAN generate inpainted map
    % ...
    dist_map(dist_mask == 1) = 0.5 + 0.5 * rand(sum(dist_mask(:)), 1);
    % Inflate
    vertex_inpaint = shapeInflate(scale_all(i), centroid_all(:,i), dist_map, dist_mask);
    % Plot
    figure(1); hold on;
    object_vertices = vertex(:, unique(object_faces(:)));  
    scatter3(object_vertices(1,:), object_vertices(2,:), object_vertices(3,:), 'filled', 'b');
    scatter3(vertex_inpaint(1,:), vertex_inpaint(2,:), vertex_inpaint(3,:), 'filled', 'r');
    axis equal
end
% save(fullfile(folder, 'scale.mat'), 'scale_all', 'centroid_all');

