function vertex_inpaint = shapeInflate(scale, centroid, dist_map, dist_mask)
% Inflate 2D distance map to a 3D point cloud.
% Input:
%   scale: scalar, used for later recovery the real distance = dist * scale, picked as the maximum value of distance
%   centroid: 3 x 1 vector, used for later recovery the real coordinates
%   dist_map: N_theta x N_phi inpainted 2D distance map by GAN
%   dist_mask: N_theta x N_phi mask labelling the missing region
% Output:
%   vertex_inpaint: 3 x V matrix, newly inpainted vertex coordinates

% -------------------------------------------------------------------------
% Step 1 Recover 3D directions from 2D mask
% ------------------------------------------------------------------------- 
N_theta = size(dist_map, 1); % image height
N_phi = size(dist_map, 2); % image width
d_theta = (pi - 0) / N_theta;
d_phi = (pi - (-pi)) / N_phi;

% Cartesian coordinates (x, y, z) to spherical coordinates (r, phi, theta)
% r = sqrt(x^2 + y^2 + z^2);
% phi = atan2(y / x); range [-pi, pi]
% theta = acos(z / r); range [0, pi], 0-north pole, pi-south pole
% Spherical coordinates (r, phi, theta) to Cartesian coordinates (x, y, z)
% x = r * sin(theta) * cos(phi);
% y = r * sin(theta) * sin(phi);
% z = r * cos(theta);
% In our case, the normal vector are already normalized so r = 1
[i,j] = find(dist_mask == 1);
theta = 0 + (i - 0.5) * d_theta;
phi = -pi + (j - 0.5) * d_phi;

% -------------------------------------------------------------------------
% Step 2 Recover 3D coordinates from distance map
% ------------------------------------------------------------------------- 
distance = scale * dist_map(sub2ind(size(dist_map),i,j)); % linear indexing
vertex_inpaint = [distance .* sin(theta) .* cos(phi) distance .* sin(theta) .* sin(phi) distance .* cos(theta)];
vertex_inpaint = centroid + vertex_inpaint'; % centroid offset

end

