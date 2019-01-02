function plot_face_color(figure_title, edge_visibility)
global PLOT_FIG vertex faces face_colors;
% Plot mesh with per face coloring
% Input:
%   figure_title: figure title
%   edge_visibility: 0-face only, 1-display edge
figure(PLOT_FIG); PLOT_FIG = PLOT_FIG + 1;
p = patch('vertices',vertex','faces',faces','FaceVertexCData', face_colors, 'FaceColor', 'flat'); % 'flat' for per-face, 'interp' for per-vertex
lighting phong, shading faceted; % 'faceted' for per-face, 'interp' for per-vertex
colormap hot, caxis([0 1]); % fix colormap range
axis equal, axis off;
camlight, cameratoolbar;
if ~isempty(figure_title)
    title(figure_title);
end
if edge_visibility == 0
    p.LineStyle = 'none';
end
end

