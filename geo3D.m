close all;

%% Load particle information
load('particles.mat');
% 'faces' and 'vertex' are still global variables

%% Determine shape completeness
tic

object_no = length(particle_faces);
completion = zeros(object_no, 1);
for i = 1 : length(particle_faces)
    object_faces = particle_faces{i};
    completion(i) = shapePercentage(vertex, object_faces);
end

fprintf('Shape completeness: %f seconds\n', toc);