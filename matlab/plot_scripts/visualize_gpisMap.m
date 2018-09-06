disp('drawing...');
% Note: Displaying plots sometimes causes an (uknown) error.
%       If you want to save the whole sequence,
%       then start Matlab with the '-nodisplay' option 
%       and add codes that save the images in a folder.

gp_bias = 0.2;
sensor_offset = [0.08;0];
head = 0.5*[0.25 0 -0.25 0.25; 0 1 0 0 ];

res = mexGPisMap('test',xtest);
arr = [0 1; -1 0]*Rot*head+tr;
fval = res(1,:) + gp_bias;
va = reshape(res(4,:),size(xg));
% res(1,:) : fval
% res(1:2,:) : grad
% res(4,:) : var

figure(1), hold off;
gcf_position = get(gcf,'Position');
gcf_position(3:4) = [(xmax-xmin)*20+10 (ymax-ymin)*20+10];
set(gcf,'Color','w','Position',gcf_position);
ax1 = gca;  hold off;

% sdf
h = pcolor(xg,yg,reshape(fval,size(xg)));  hold on;
set(h,'EdgeColor','none'); axis equal;
set(ax1,'CLim',[-0.4 0.4]);
axis([xmin xmax ymin ymax])

% surface extraction
[~,temp] = isocontour(reshape(fval,size(xg)),0);
Ind = (round(temp(:,2))-1)*size(xg,1) + round(temp(:,1));
valid= find(Ind>0 & Ind<=prod(size(xg)));
valid2 = find(va(Ind(valid)) < 0.4);
Ind = valid(valid2);
tx_draw = ((xmax-xmin)*(temp(Ind,2)-0.5)/size(xg,2)) + xmin;
ty_draw = ((ymax-ymin)*(temp(Ind,1)-0.5)/size(yg,1)) + ymin;
plot(tx_draw,ty_draw,'r.','MarkerSize',5);

% variance
ax2 = axes; hold off;
h = pcolor(ax2,xg,yg,ones(size(va))); hold on;
set(h,'EdgeColor','none');
alpha(ax2,va);
linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
axis(ax2,'equal');
colormap(ax1,parula)
colormap(ax2,[1 1 1])
set([ax1 ax2],'Units','pixels');
set([ax1,ax2],'Position',[5 5 (xmax-xmin)*20 (ymax-ymin)*20]);
set([ax1,ax2],'XLim',[xmin xmax]);
set([ax1,ax2],'YLim',[ymin ymax]);
set([ax1,ax2],'XTick',{});
set([ax1,ax2],'YTick',{});
drawnow;

% measurement
valid = find((ranges(nframe,:)'<3e1) & (ranges(nframe,:)'>2e-1) & (~isinf(ranges(nframe,:)')));
XY = polar2xy(thetas(valid),ranges(nframe,(valid))');
XY(1,:) = XY(1,:) + sensor_offset(1);
XY_ref = Rot*XY + tr;
for t=1:numel(valid)
    plot([tr(1) XY_ref(1,t)],[tr(2) XY_ref(2,t)],'g-');
end

% robot pose
headpatch = patch( arr(1,:),arr(2,:),'w');
set(headpatch, 'EdgeColor','r','LineWidth',2);
plot(poses(initframe:skip:nframe,1),poses(initframe:skip:nframe,2),'k-','LineWidth',1);
drawnow;
set(gcf,'Position',gcf_position);
