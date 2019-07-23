function obj = nonparametric_fitting(obj, span)


% roughness: 

for k = 1:obj.NumObjects
phi = obj.objects(k).polar(:, 1);
% phi = [phi(1); phi(1); phi; phi(end); phi(end)];
ro =  obj.objects(k).polar(:, 2);
% ro = [ro(1); ro(1); ro; ro(end); ro(end)];
ctd = obj.objects(k).centroid;
%% Use "cvpartition" and "crossval"

% opt = statset('UseParallel','always');

% Npts = length(phi);

% span = 0.022;
% linspace(.001,.5,num);
% 
% sse = zeros(size(spans));
% cp = cvpartition(Npts,'k',3);
% % mylowess([phi(train), ro(train)],phi(test),0.01)
% for j=1:length(spans)
%     f = @(train,test) norm(test(:,2) - mylowess(train,test(:,1),spans(j)))^2;
%     sse(j) = sum(crossval(f,[phi,ro],'partition',cp, 'Options', opt));
% %         sse(j) = sum(crossval(f,[phi,ro],'partition',cp));
% end
% 
% [minsse,minj] = min(sse);
% span = spans(minj);
% % span = 0.0091;

ro_mean = smooth(phi,ro,span,'loess');
%% the figure 7 of the geotechnique paper
% figure
% plot(phi*180/pi, ro, 'k-')
% 
% plot(phi*180/pi,ro,'k','linestyle','-','linewidth',1);
% hold on
% plot(phi*180/pi,ro_mean,'k','linestyle','- -','linewidth',1);
% xlim([0,360]); ylim([150, 450]);
% % axis equal
% xlabel('\theta (degrees)', 'FontSize', 15);
% ylabel('\rho (pixels)','FontSize', 15)
% ax = gca;
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
% ax.XTick = [0 45 90 135 180 225 270 315 360];
% 
% 
% ax.PlotBoxAspectRatio = [1, 1, 1];
% 
% ax.FontSize = 12;
% text(2.53, 3.7,'(b)', 'FontSize', 15);
% L2 = legend('Measurements', 'LOESS: \alpha = 0.0150');


% length_loess = arclength(phi,ro_mean,  'spline');

% Roughness = sqrt(mean((ro-ro_mean).^2))/length_loess;


% ro_mean = mylowess([phi,ro],phi,span);
% x = linspace(min(phi),max(phi));
% line(x,mylowess([phi,ro],x,span),'color','k','linestyle','-', 'linewidth',2)
% 
% legend('Noisy Data Sample','LOWESS: Optimal Span',2)


%%
% opt = statset('UseParallel','always');
% hold off
% scatter(phi, ro)
% 
% f = @(xy) mylowess(xy, phi, span);
% [yboot2, bootsam] = bootstrp(100, f, [phi,ro], 'Options', opt);
% yboot2 = yboot2';
% % yboot2 = bootstrp(10, f, [phi,ro], 'Options', opt)';
% meanloess = mean(yboot2,2);
% 
% h1 = line(phi, meanloess,'color','b','linestyle','-','linewidth',2);
% 
% stdloess = std(yboot2,0,2);
% mean(stdloess)
% h2 = line(phi, meanloess+2*stdloess,'color','r','linestyle','--','linewidth',2);
% h3 = line(phi, meanloess-2*stdloess,'color','r','linestyle','--','linewidth',2);
% 
% L5 = legend('Localized Regression','Confidence Intervals');

% figure
% hold on
% % for i =1:size(bootsam, 2);
% %     plot(phi,...
% %         mylowess([phi(bootsam(:,i)), ro(bootsam(:,i))],phi,span));
% % end

% [X2, Y2] = mypolar2carte(meanloess, phi, ctd);
[X, Y] = mypolar2carte(ro_mean, phi, ctd, 2);

X(end) = X(1);
% X(1) = (X(end) + X(1))/2;
Y(end) = Y(1);
% Y(1) =(Y(end)+ Y(1))/2;
obj.objects(k).cartesian = [X, Y];
% obj.objects(k).roughness = Roughness;
end

end