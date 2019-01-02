% Blog: http://www.voidcn.com/article/p-mitplaal-mw.html
% https://blog.csdn.net/moneyhoney123/article/details/78454837

%% Control panel
clc;
close all;

inFolderName = '09_26_2018'; % user-define

% Turn on/off the following options in sequence to complete the analysis
PLOT = false;           % figure windows
POINT_READ = true;     % read point cloud data
GROUND_REMOVAL = true; % remove ground points based on principal direction
CAMERA_READ = true;     % read camera parameters
REPROJECTION = true;    % 3D reprojection point cloud onto different 2D image planes

%% Read point cloud
if POINT_READ % takes around 5 secs
    
inFileName = strcat('./Results/', inFolderName, '/', strcat(inFolderName, '.nvm.cmvs/00/models/option-0000.ply')); 
outFolderName = strcat('./Results/', inFolderName, '/Point Cloud');
if ~exist(outFolderName, 'dir')
	mkdir(outFolderName);
end

% Total number of points
[~ ,~ , N] = textread(inFileName, '%s %s %u', 1, 'headerlines', 2);
% Read in all points
[x, y, z, nx, ny, nz, R, G, B] = textread(inFileName, '%f %f %f %f %f %f %u %u %u', N, 'headerlines', 13);
points = [x y z]; % n x 3 matrix
normals = [nx ny nz]; % n x 3 matrix
colors = [R G B]; % n x 3 matrix
% Convert RGB to grayscale
I = 0.2989 * R + 0.5870 * G + 0.1140 * B;

end

%% Ground removal
if GROUND_REMOVAL
% -------------------------------------------------------------------------
% Step 1: Find principal plane (floor/ground) from normal vectors
% -------------------------------------------------------------------------
% Cartesian coordinates (x, y, z) to spherical coordinates (r, phi, theta)
% r = sqrt(x^2 + y^2 + z^2);
% phi = atan2(y / x);
% theta = acos(z / r);
% Spherical coordinates (r, phi, theta) to Cartesian coordinates (x, y, z)
% x = r * sin(theta) * cos(phi);
% y = r * sin(theta) * sin(phi);
% z = r * cos(theta);
% In our case, the normal vector are already normalized so r = 1
theta = acos(nz);
phi = atan2(ny, nx);

interval = 1; % 1 deg

votes = hist3([phi theta], [360/interval 180/interval]); 
votes = votes';
temp = padarray(votes, [1 1], 'post'); % pad the matrix by 1 to allow correct labels (matrix's index start at 1 instead of 0)

if PLOT
% Plot histogram
figure(1)
hist3([phi theta], [360/interval 180/interval], 'EdgeColor', 'none'); % theta range [0, pi], phi range [-pi, pi]
xlabel('\phi (rad)'); ylabel('\theta (rad)'); zlabel('Bin Count'); title('Distribution of normal vectors in spherical coordinates');
set(gca,'XTick',-pi:pi/2:pi); set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
set(gca, 'YDir', 'reverse'); set(gca,'YTick',0:pi/4:pi); set(gca,'YTickLabel',{'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
colormap hot;
daspect([1 1 250]);

figure(2)
imagesc(temp), colormap hot, axis image;
set(gca,'XTick',1:90:361); set(gca,'XTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
set(gca, 'YDir', 'reverse'); set(gca,'YTick',1:45:181); set(gca,'YTickLabel',{'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
xlabel('\phi (rad)'); ylabel('\theta (rad)'); title('Normal vectors in \phi-\theta space');
end

[~, idx] = max(votes(:));
[idx_i, idx_j] = ind2sub(size(votes), idx);
% Principal direction of the normal vectors (floor/ground normal)
phi_p = (idx_j-180)/180*pi;
theta_p = idx_i/180*pi;
if PLOT
figure(1), hold on, scatter3(phi_p, theta_p, votes(idx_i, idx_j), 'sb', 'LineWidth', 2);
figure(2), hold on, plot(idx_j, idx_i, 'sb', 'LineWidth', 2);
end
% Convert to Cartesian coordinates
np_x = 1 * sin(theta_p) * cos(phi_p);
np_y = 1 * sin(theta_p) * sin(phi_p);
np_z = 1 * cos(theta_p);
fprintf('Principal Direction: phi = %3.2f, theta = %3.2f\n', phi_p, theta_p);
fprintf('Normal Vector: (%3.2f, %3.2f, %3.2f)\n', np_x, np_y, np_z);

% -------------------------------------------------------------------------
% Step 2: Remove floor/ground points from cloud
% -------------------------------------------------------------------------
% General idea:
% After knowing the principal normal direction, we can calculate for each
% point P their distance to origin O along the normal direction, as the
% projection to the principal normal vector.
% Suppose point & normal vector form are p(x,y,z), n(nx, ny, nz), then
% d = |PO|cos<p,n> = dot(p,n)/|n| = dot(p,n)
% Use direction cosines of normal vector to define a plane
% x*cosa + y*cosb + z*cosc = d (d is distance from origin to the plane)
points = [x y z]; % n x 3 matrix
normal_p = [np_x np_y np_z]'; % principal normal vector
distances = points * normal_p; % (n x 3) * (3 x 1) = (n x 1) dot products
% Origin is usually at above the objects, so distances are negative
dist_min = min(distances(:)); % farest (floor)
dist_max = max(distances(:)); % nearest
% Based on the prior that floor pixels will always be the farest ones with
% repsect to origin, we can simply remove those
dist_percentage = 0.2;
dist_threshold = dist_min + (dist_max - dist_min) * dist_percentage;
survivors = points(distances > dist_threshold, :);

if PLOT
figure(3), hold on; % raw point cloud
scatter3(points(:,1), points(:,2), points(:,3), 'b');
scatter3(np_x, np_y, np_z, 'r','filled');
scatter3(0, 0, 0, 'b','filled');
quiver3(0, 0, 0, np_x, np_y, np_z, 'r'); % draw principal normal
% [plane_x, plane_y] = meshgrid(-2:0.5:2);
% plane_z = (dist_min - np_x * plane_x - np_y * plane_y) / np_z;
% surf(plane_x, plane_y, plane_z); zlim([-5 5]); % draw ground plane
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Raw point cloud');
daspect([1 1 1]);
view(-180.3, -89.2); % no plane view
% view(86.5, -50.8); % plane view

figure(4), hold on; % ground-removal point cloud
scatter3(survivors(:,1), survivors(:,2), survivors(:,3), 'b');
scatter3(np_x, np_y, np_z, 'r','filled');
scatter3(0, 0, 0, 'b','filled');
quiver3(0, 0, 0, np_x, np_y, np_z, 'r'); % draw principal normal
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Ground-removal point cloud');
daspect([1 1 1]);
view(-180.3, -89.2); % no plane view
end

% Save the ground-removed points cloud
points = survivors;
normals = normals(distances > dist_threshold, :);
% I_survivors = I(distances > dist_threshold, :) / 255;
% result = [survivors I_survivors];
% dlmwrite(fullfile(outFolderName, 'ground_removal.txt'), result, 'delimiter', ' ', 'newline', 'pc');
end

%% Read camera parameters
if CAMERA_READ 
    
inFileName = strcat('./Results/', inFolderName, '/', strcat(inFolderName, '.nvm')); 
outFolderName = strcat('./Results/', inFolderName, '/Camera Parameter');
if ~exist(outFolderName, 'dir')
	mkdir(outFolderName);
end

% Total number of cameras
N_cam = textread(inFileName, '%u', 1, 'headerlines', 2);
% Read in all cameras <File name> <focal length> <quaternion WXYZ> <camera center> <radial distortion> 0
[~, f, Q_w, Q_x, Q_y, Q_z, C_x, C_y, C_z, D, ~] = textread(inFileName, '%s %f %f %f %f %f %f %f %f %f %f', N_cam, 'headerlines', 3);
% Wrap parameters into struct
for i = 1 : N_cam
    cameras(i).focal = f(i); % focal length
    cameras(i).quat = quatnormalize([Q_w(i) Q_x(i) Q_y(i) Q_z(i)]); % w + x*i + y*j + z*k, w-real part, x/y/z-imaginary part. Don't forget to normalize!
    cameras(i).center = [C_x(i) C_y(i) C_z(i)];
end

end

%% 3D Reprojection problem
if REPROJECTION
% Idea: 
% Step 1: Recover camera parameters (extrinsic [R|t] & intrinsic [K] matrices)
% Step 2: Re-project 3D point cloud onto camera/image coordinates (with back-face culling)
% Step 3: Apply detected boundary restriction on the re-projected pixels
% Step 4: Partition and label the points
% Step 5: Loop Step 1 ~ 4 for all cameras
% Note: re-projection looks good. Task next is to remove some points based
% on Z depth hidden surface removal
rawImages = dir(fullfile(strcat('./ImageSet/', inFolderName, '/', inFolderName), '*.JPG'));
for i = 1 : 1%N_cam
% Step 1 ------------------------------------------------------------------
    % Quaternion Q(w + x*i + y*j + z*k) to rotation matrix R
    % Advantage of quaternion is to represent rotation in a straightforward
    % way. Rotation with quaternion is defined by an axis of rotation (x0, y0, z0)
    % and an angle of rotation about the axis (theta). Then the quaternion Q:
    % Q = cos(theta/2) + sin(theta/2)*(x0*i + y0*j + z0*k) = w + x*i + y*j + z*k
    % Note: when read from file, the quaternion is normalized first.
    % Ref: http://run.usc.edu/cs520-s12/quaternions/quaternions-cs520.pdf
    % http://work.thaslwanter.at/Kinematics/html/04_Quaternions.html
    % Extrinsic parameter (translation, rotation)
    q = cameras(i).quat;
    R = quat2rotm(q); % Rotation matrix, from quaternion
    t = - R * cameras(i).center'; % translation matrix, - R * center. Don't forget R and '-'!
    Extrinsic = [R t]; % extrinsic matrix [R t] = [R -RC] where C is camera center
    % Intrinsic parameter (focal length, image center)
    imageName = strcat('./Results/', inFolderName, '/', strcat(inFolderName, '.nvm.cmvs/00/visualize/', num2str(i - 1, '%08u'), '.jpg')); % use the renamed and rotated image file duplicated by PMVS. % Previous BUG here! PMVS folder use 0-based index...so (i-1)!
    info_raw = imfinfo(fullfile(rawImages(i).folder, rawImages(i).name));
    info = imfinfo(imageName);
    c_x = info.Width / 2; % principal point of image (usually image center)
    c_y = info.Height / 2;
    focal = cameras(i).focal;
    horizontal_FOV = 57.716; % in degree % https://developer.apple.com/library/archive/documentation/DeviceInformation/Reference/iOSDeviceCompatibility/Cameras/Cameras.html#//apple_ref/doc/uid/TP40013599-CH107-SW22
    vertical_FOV = horizontal_FOV * c_y/c_x;
    focal_x = c_x/(tan(deg2rad(horizontal_FOV/2)));%cameras(i).focal; % * 29/35; % info_raw.DigitalCamera.FocalLengthIn35mmFilm * 100; % we should search more about equivalent 35mm film focal length and why *100 is needed
    focal_y = c_y/(tan(deg2rad(vertical_FOV/2)));%cameras(i).focal;
    %Intrinsic = [focal_x 0 c_x; 0 focal_y c_y; 0 0 1];
    Intrinsic = [focal 0 c_x; 0 focal c_y; 0 0 1]; % after the bug fix, using SfM focal length has perfect results
    % Back-face culling
    cam_view = points - repmat(cameras(i).center, size(points,1), 1); % view vector of camera
    culling = dot(normals, cam_view, 2); % n x 3 and n x 3, treat row (dim = 2) as vectors
    survivors = points(culling < 0, :); % dot product < 0 means the point normal doesn't face the camera
    survivor_colors = colors(culling < 0, :);
    fprintf('Camera %u:\n%u points removed during back-face culling.\n', i, size(points,1) - size(survivors, 1));
    % Re-projection
    survivors = [survivors ones(size(survivors,1),1)]; % add homogeneous coordinates
    % Hidden surface/point removal based on depth (just like Z buffer method in Computer Graphics)
    camera_space = (Extrinsic * survivors')'; % points in camera coordinates
    survivor_depths = camera_space(:, 3);
    projections = (Intrinsic * Extrinsic * survivors')'; % n x 3 = (3 x 3 * 3 x 4 * 4 x n)'
    projections(:, 1) =  round(projections(:, 1) ./  projections(:, 3)); % divide by homogeneous coordinates
    projections(:, 2) =  round(projections(:, 2) ./  projections(:, 3)); % [x/z y/z]
    % Remove re-projected points that are out camera FoV (Field of View)
    inPlane = projections(:, 1) > 0 & projections(:, 1) < info.Width & projections(:, 2) > 0 & projections(:, 2) < info.Height; 
    fprintf('%u points removed due to out-of-plane projection.\n', size(projections,1) - sum(inPlane(:)));
    projections = projections(inPlane, :);
    % Remove hidden surface/point
    survivor_depths = survivor_depths(inPlane);
    survivor_colors = survivor_colors(inPlane);
    hidden_points = ones(length(survivor_depths), 1); % 1-seen; 0-hidden
    pixel_location = zeros(info.Height, info.Width);
    % Display re-projected 2D image
    canvas = ones(info.Height, info.Width, 3);
    removal = ones(info.Height, info.Width, 3);
    for p = 1 : length(survivor_depths)
        % Without hidden surface removal
        % canvas(projections(p, 2), projections(p, 1), :) = survivor_colors(p, :) / 255;
        % With hidden surface removal
        previous = pixel_location(projections(p, 2), projections(p, 1));
        if previous == 0 % first visitor (default value 0)
            pixel_location(projections(p, 2), projections(p, 1)) = p;
            canvas(projections(p, 2), projections(p, 1), :) = survivor_colors(p, :) / 255;
        else % a previous point visits and its index 'previous' was recorded in pixel_location, compare two depths
            if survivor_depths(previous) > survivor_depths(p) % remove previous (previous is deeper thus covered)
                hidden_points(previous) = 0;
                canvas(projections(p, 2), projections(p, 1), :) = survivor_colors(p, :) / 255; % update pixel color
            else % remove current
                hidden_points(p) = 0;
            end
            removal(projections(previous, 2), projections(previous, 1), :) = survivor_colors(p, :) / 255;
        end
    end
    fprintf('%u points removed due to hidden surface removal.\n', size(projections,1) - sum(hidden_points(:)));
    imshow(removal);
    % Closure: hidden surface removal attempt here is not a perfect
    % solution. Maybe I should try to rely on multiple images to correct
    % this. For example, if I apply the boundary restriction here. Due to
    % the overlapping effect which hasn't fully resolved by the depth
    % buffer method yet, I may wrongly label some points from the rock
    % beneath into this surface rock. I can keep it for now. But later when
    % an image from another angle indicates that they should be labelled
    % differently.
    % Notice that: the inclusion/mistaken of overlapping points is always
    % an overestimation (includes some points that don't belong to you).
    % So, when we observe from other views that this point should not
    % belong to this group anymore, it's ALWAYS true.
    % If following such logic, the hidden surface removal SHOULDN'T be
    % applied.
    % Another thought is to abandone the controversial points, and leave
    % them for object completion step.
    % Another thought is to store the depths of each point to the camera
    % in each image, and do a 3D verification among all images.
    %img_overlap = imfuse(canvas, imread(imageName));
    %imshow(img_overlap);
    %imwrite(img_overlap, strcat(outFolderName, '/', num2str(i, '%03u'), '.jpg'));
    imshowpair(canvas, imread(imageName), 'montage');
    fprintf('Complete\n');
end
   
end

%% K-mean
KMEAN = false;
if KMEAN
close all;
k = 9;
[idx, centers] = kmeans(survivors, k);
% scatter3(survivors(:,1), survivors(:,2), survivors(:,3), 'b');
%scatter3(centers(:,1), centers(:,2), centers(:,3), 'r');
[uniqueGroups, uga, ugc] = unique(idx); 
colors = brewermap(length(uniqueGroups),'Set1');  %or any other way of creating the colormap
markersize = 20;   %change to suit taste
figure(5), hold on;
scatter3(survivors(:,1), survivors(:,2), survivors(:,3), markersize, colors(ugc,:));
daspect([1 1 1]);
view(-180.3, -89.2); % no plane view
end

% Finding
% shadow or edges can't be obtained from sfm, because SIFT features can't
% be effectively extractly from a shadow region
% Maybe we can try segment/boundary detection each image, and use the
% N-view match matrix to create a 3D skeleton for boundaries. Then we apply
% these boundary skeleton as a restriction for the k-means algorithm.

