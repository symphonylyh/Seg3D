close all;

SHAPE_COMPLETE = false;
SHAPE_FLATTEN = true;

%% Load particle information
% load('particles.mat');
% 'faces' and 'vertex' are still global variables
addpath(genpath('toolbox_graph'));

name = 'rockGenerator/samples_5/rock_001';
[vertex,faces] = read_mesh(name);

%% Determine shape completeness
if SHAPE_COMPLETE
    tic

    object_no = length(particle_faces);
    completion = zeros(object_no, 1);
    for i = 1 : length(particle_faces)
        object_faces = particle_faces{i};
        completion(i) = shapePercentage(vertex, object_faces);
    end

    fprintf('Shape completeness: %f seconds\n', toc);
end

%% Generate flatten shape for training
if SHAPE_FLATTEN
    completion = shapeFlatten(vertex, faces);
end

% Issue:
% Scale the distance into grayscale 0~1, but for some uniform particle
% (e.g. a sphere), the grayscale value will differ a lot although the
% absolute distance is tiny.