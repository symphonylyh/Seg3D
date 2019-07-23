function obj = discrete_boundary(cc)
s  = regionprops(cc, 'centroid', 'Area','Perimeter','MajorAxisLength','MinorAxisLength', 'Orientation');
sz=[cc.ImageSize];

% boundaries = cell(cc.NumObjects, 1);
% obj.cetroid = zeros(cc.NumObjects, 2);
for k = 1:cc.NumObjects
   BW = zeros(sz);
    ctd = round([s(k).Centroid]);
    BW(cc.PixelIdxList{k}) = 1;
%     Area = bwarea(BW);
    ContourImage = bwmorph(BW,'remove');
    %     find the beginning point
    [~,x0]=find(ContourImage(ctd(2), :),1, 'first');
    
    contour = bwtraceboundary(ContourImage, [ctd(2), x0], 'W', 8,...
        Inf,'counterclockwise');
    X=contour(:,2);
    Y=contour(:,1);
    
    % conver the cartesian coordinates to polar coordinates
    
    [ro, phi] = mycarte2polar(X, Y, ctd, 2);
    obj.objects(k).polar = [phi, ro];
    obj.objects(k).centroid = s(k).Centroid;
    obj.objects(k).orientation = s(k).Orientation;
    obj.objects(k).area = s(k).Area;
    obj.objects(k).perimeter = s(k).Perimeter;
    obj.objects(k).rawXY = [X, Y];
    obj.objects(k).d1d2 = [s(k).MajorAxisLength, s(k).MinorAxisLength];
end
    obj.NumObjects = cc.NumObjects;
    obj.ImageSize = [sz];
end
