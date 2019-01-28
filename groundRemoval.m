%% Control panel
clc;
close all;

% Turn on/off the following options in sequence to complete the analysis
PLOT = true;           % figure windows
POINT_READ = true;     % read point cloud data
GROUND_REMOVAL = true; % remove ground points based on principal direction

%% Read point cloud
if POINT_READ % takes around 5 secs

inFileName = 'option-0000.ply';

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
dist_percentage = 0.3;
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

pcwrite(pointCloud(points), 'option-0000-removal.ply');
end
