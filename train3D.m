% To compile MATLAB script to executable, you need to first install two large toolbox:
% 1. MATLAB Compiler SDK // to compile to a standalone
% executable/application package
% 2. MATLAB Runtime // to run the executable (this is free, can work on any computer without MATLAB license)
% Ref: https://www.nas.nasa.gov/hecc/support/kb/compiling-matlab-scripts-into-executables-to-reduce-the-use-of-licenses_527.html
% A little too complicated...

close all;

DEBUG = false;
BATCH = true;

addpath(genpath('toolbox_graph'));

%% Load rock information
folder = 'rockGenerator/samples_50';
files = [];

if DEBUG
    files(1).name = 'rock.obj';
end

if BATCH
    files = dir(folder);
    files = files(~ismember({files.name},{'.','..','.DS_Store'}));
end

scale_all = zeros(length(files), 1);
centroid_all = zeros(3, length(files));
for f = 1 : length(files)
    % For Blender .obj file, here I slightly modified the read_obj.m in toolbox_graph
    [vertex,faces, ~] = read_obj(fullfile(folder, files(f).name));
    [scale_all(f), centroid_all(:,f), dist_map, dist_mask] = shapeDeflate(vertex, faces);
    dist_map = imresize(dist_map, [256 256]); % for DeepFill, at least 64x64 images are required, I don't know why
    dist_map(dist_map < 0) = 0;
    dist_map = repmat(dist_map, 1, 1, 3); % for DeepFill 3-channel images are required
    imwrite(dist_map, fullfile(folder, strcat(num2str(f, '%04.f'), '.png')));
end
% save(fullfile(folder, 'scale.mat'), 'scale_all');

% Issue:
% Scale the distance into grayscale 0~1, but for some uniform particle
% (e.g. a sphere), the grayscale value will differ a lot although the
% absolute distance is tiny.
