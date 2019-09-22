close all;

global SINGLE COMPLETENESS MAP;
SINGLE = true;
COMPLETENESS = false;
MAP = true;

global PLOT PLOT_FIG plot_volume;
PLOT = true;
PLOT_FIG = 1;
    plot_volume = 1;

%% Load particle and scale information
if SINGLE
    % Single particle mode (with manually cleaned .off mesh)
    NAME = 'mesh/03_20_2019/RR6_keystone_2';
    addpath(genpath('toolbox_graph'));
    [vertex,faces] = read_mesh(NAME);
    particle_points{1} = vertex';
    particle_faces{1} = faces';
    object_no = 1;
else
    % Normal mode (with mesh info from seg3D.m steps)
    NAME = 'mesh/01_06_2019/multiple_02'; 
    load(strcat(NAME, '.mat'));
    % 'faces' and 'vertex' are still global variables
    object_no = length(particle_faces);
end

[real_scale, sfm_scale] = textread(strcat(NAME, '.txt'), '%f %f');
scale = mean(real_scale ./ sfm_scale); % convert sfm scale to cm scale

%% Display segmented particle(s) and volume calculation
plot_volume_separate = 0; % plot particles on different figures
plot_volume_allinone = 1; % plot particles on one figure
shrink = 0.4; % shrink factor for boundary(). 0-convex hull, 1-compact, 0.5-default

if PLOT && plot_volume && plot_volume_allinone
    % Particles subplotted on one figure
    figure(PLOT_FIG); PLOT_FIG = PLOT_FIG + 1;
    fig_row = ceil(object_no/2);
    fig_col = 2 * 2;
end

volume_raw = zeros(object_no, 1); % volume of incomplete particle
completion = zeros(object_no, 1); % percentage of shape completeness
for i = 1 : object_no
    object_vertices = particle_points{i};
    object_faces = particle_faces{i};
    
    % Calculate volume encompassed by a set of points
    % volume_raw(i) = volumeFromPoints(object_vertices);
    [B, volume_raw(i)] = boundary(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3), shrink); 
    
    % Calculate percentage of shape completeness
    if COMPLETENESS
        completion(i) = shapePercentage(vertex, object_faces);
    end
    
    % Visualize
    if PLOT && plot_volume
        if plot_volume_separate
            figure(PLOT_FIG); PLOT_FIG = PLOT_FIG + 1;
            title(strcat('Particle ', num2str(i)));
            % point cloud
            subplot(1,2,1);
            scatter3(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3));
            axis equal off;
            % volume hull
            subplot(1,2,2);
            trisurf(B,object_vertices(:,1),object_vertices(:,2),object_vertices(:,3),'Facecolor','red','FaceAlpha',1)
            % shading
            view(3);
            light; lighting phong;
            camlight('left');
            shading interp;
            axis equal off;
        end
        if plot_volume_allinone
            % point cloud
            subplot(fig_row,fig_col,2*i-1);
            scatter3(object_vertices(:,1), object_vertices(:,2), object_vertices(:,3), 2);
            axis equal off;
            % volume hull
            subplot(fig_row,fig_col,2*i);
            trisurf(B,object_vertices(:,1),object_vertices(:,2),object_vertices(:,3),'Facecolor','red','FaceAlpha',0.1)
            axis equal off;
        end
    end
end

% Apply scale
volume_raw = volume_raw * scale^3; % in cm^3
mass_raw = volume_raw * 2.65; % in g (assume density = 2.65g/cm^3)

%% Generate incomplete distance map
if MAP
    load('particles.mat');
    folder = 'fooling_set';
    object_no = length(particle_faces);
    scale_all = zeros(length(particle_faces), 1); 
    centroid_all = zeros(3, length(particle_faces)); 
    for i = 1 : length(particle_faces)
    %for i = 9 : 9
        object_faces = particle_faces{i};
        % Deflate
        [scale_all(i), centroid_all(:,i), dist_map, dist_mask] = shapeDeflate(vertex, object_faces); 
        % imwrite(dist_map, fullfile(folder, strcat(num2str(i, '%04.f'), '.png')));
        % imwrite(dist_mask, fullfile(folder, strcat(num2str(i, '%04.f'), '_mask.png')));
        % GAN generate inpainted map
        % ...
        dist_map(dist_mask == 1) = 0.5 + 0.5 * rand(sum(dist_mask(:)), 1);
        % test
        dist_map = imread('output.png');
        dist_map = dist_map(:,:,1);
        dist_map = imresize(dist_map, [32 32]);
        dist_mask = logical(dist_mask);
        dist_map(dist_mask) = dist_map(dist_mask) * 1.1;
        dist_map = mat2gray(dist_map);
        % Inflate
        vertex_inpaint = shapeInflate(scale_all(i), centroid_all(:,i), dist_map, dist_mask);
        % Plot point cloud
        figure(1); hold on;
        object_vertices = vertex(:, unique(object_faces(:)));  
        point_size = 5; % small size (5) for group view, bigger (30) for single view
        scatter3(object_vertices(1,:), object_vertices(2,:), object_vertices(3,:), point_size, 'filled', 'b');
        scatter3(vertex_inpaint(1,:), vertex_inpaint(2,:), vertex_inpaint(3,:), point_size, 'filled', 'r');
        axis equal
        % Plot mesh
        figure(2); hold on;
        merge_vertices = [object_vertices vertex_inpaint];
        [B, ~] = boundary(merge_vertices(1, :)', merge_vertices(2, :)', merge_vertices(3, :)', 0.5); 
        trisurf(B,merge_vertices(1,:), merge_vertices(2,:), merge_vertices(3,:),'Facecolor','red','FaceAlpha',1)
        view(3);
        light; lighting phong;
        camlight('left');
        shading interp;
        colormap gray;
        axis equal
    end
    
    % save(fullfile(folder, 'scale.mat'), 'scale_all', 'centroid_all');
end