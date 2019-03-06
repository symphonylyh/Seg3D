close all;

COMPLETNESS = false;

%% Load particle information
clear all;
% load('particles.mat');
load('mesh/01_06_2019/01.mat');
% 'faces' and 'vertex' are still global variables

volumes_raw = zeros(object_no, 1);
particle_points{object_no} = [];
particle_faces{object_no} = [];
% Particles on different figures
for i = 1 : object_no
    object_faces = faces(:, object_set(:,i));
    object_vertices = vertex(:, unique(object_faces(:)))';
    particle_points{i} = object_vertices;
    particle_faces{i} = object_faces;
    [B, volumes_raw(i)] = boundary(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3), 0.1);
    % plot
    if PLOT && plot_volume
        figure(PLOT_FIG); PLOT_FIG = PLOT_FIG + 1;
        title(strcat('Particle ', num2str(i)));
        subplot(1,2,1);
        scatter3(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
        axis equal off;

        subplot(1,2,2);
        trisurf(B,object_vertices(:,1),object_vertices(:,2),object_vertices(:,3),'Facecolor','red','FaceAlpha',0.1)
        axis equal off;
    end
end

% Particles subplotted on one figure
% if PLOT && plot_volume
%     figure(PLOT_FIG); PLOT_FIG = PLOT_FIG + 1;
%     % Arrange subplots
% %     fig_row = 2 * round(sqrt(object_no));
% %     fig_col = 2 * ceil(object_no/(fig_row/2));
%     fig_row = ceil(object_no/2);
%     fig_col = 2 * 2;
% end
% for i = 1 : object_no
%     object_faces = faces(:, object_set(:,i));
%     object_vertices = vertex(:, unique(object_faces(:)))';
%     particle_points{i} = object_vertices;
%     particle_faces{i} = object_faces;
%
%     % Calculate volume encompassed by a set of points
%     % volume_raw(i) = volumeFromPoints(object_vertices);
%     [B, volumes_raw(i)] = boundary(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
%     
%     % plot
%     if PLOT && plot_volume
%         subplot(fig_row,fig_col,2*i-1);
%         scatter3(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
%         axis equal off;
% 
%         subplot(fig_row,fig_col,2*i);
%         trisurf(B,object_vertices(:,1),object_vertices(:,2),object_vertices(:,3),'Facecolor','red','FaceAlpha',0.1)
%         axis equal off;
%     end
% end

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

