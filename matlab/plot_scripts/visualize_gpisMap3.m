% test 3D GPIS using depth image

% for visualization
gcf_size=[360 450];
gcf_color = 'k';
fig_visible = 'on';
param.bias = 0.2; % This must be the same to xxxxx

res  = mexGPisMap3('test',xtest1);
fval = res(1,:);
val = reshape(fval,size(xg))+param.bias;

% mesh by isosurface
[f,v] = isosurface(xg,yg,zg,val,0,'noshare');
va = mexGPisMap3('test',single(v'));
va = va(5,:);
vis_thre = 0.02;
falpha = va;
min_falpha = min(falpha);
falpha = max((vis_thre - falpha)/(vis_thre-(min_falpha+0.001)),0);

figure(3); hold off;
set(gcf,'visible',fig_visible,'renderer','opengl');
gcf_pos = get(gcf,'Position');
gcf_pose(3:4) = gcf_size;
set(gcf,'Color',gcf_color,'Position',gcf_pose);
pcshow([-100,-100,-100],'k','VerticalAxis','Z','VerticalAxisDir','Up'); hold on; % dummy point for better axis setting
p = patch('Faces',f,'Vertices',v,'FaceVertexCData',repmat([0.5 0.5 0.5],size(v,1),1),'FaceVertexAlphaData',falpha','FaceAlpha','interp'); hold on;
axis off;
isonormals(xg,yg,zg,reshape(fval,size(xg)),p)
p.EdgeColor = 'none';
clt = camlight(0,20);
lighting gouraud
axis([-0.09 0.17 -0.13 0.17 0 0.3])
shading interp

ang = 20*pi/180;
[xg2, yg2] = meshgrid(-0.05:0.01:0.13, -0.1:0.01:0.14);
xtest2 = [xg2(:)'; yg2(:)'; 0.12*ones(1,numel(xg2))];
R_ = [cos(ang) sin(ang) 0; -sin(ang) cos(ang) 0; 0 0 1];
xtest2 = R_'*(xtest2-[0.04 0.02 0]') + [0.04 0.02 0]';
xtest2 = single(xtest2);
res2 = mexGPisMap3('test',xtest2);
a2 = res2(1,:);

[yg3, zg3] = meshgrid(-0.1:0.01:0.14, 0:0.01:0.30);
xtest3 = [zeros(1,numel(yg3)); yg3(:)'; zg3(:)'];
R_ = [cos(ang) sin(ang) 0; -sin(ang) cos(ang) 0; 0 0 1];
xtest3 = R_'*xtest3 + [0.04 0 0]';
xtest3 = single(xtest3);
res3 = mexGPisMap3('test',xtest3);
a3 = res3(1,:);

h = surf(reshape(xtest2(1,:),size(xg2)),reshape(xtest2(2,:),size(xg2)),reshape(xtest2(3,:),size(xg2)));
h.CData = reshape(a2+param.bias,size(yg2));
h.EdgeColor = 'none';
h.FaceAlpha =0.6;
colormap('jet')
set(gca,'CLim',[-0.05 0.2])
h2 = surf(reshape(xtest3(1,:),size(yg3)),reshape(xtest3(2,:),size(yg3)),zg3);
h2.CData = reshape(a3+param.bias,size(yg3));
h2.EdgeColor = 'none';
h2.FaceAlpha = 0.6;
colormap('jet')
set(gca,'CLim',[-0.05 0.2])
shading interp

view(-30,30);
set(gca,'units','pixels','Position',[-50 -60 gcf_size+[100 100]]);
drawnow;
set(gcf,'Position',gcf_pose);
