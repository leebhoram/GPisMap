% GPisMap demo 2D 
%
% by Bhoram Lee 
% Aug. 28. 2018
clearvars
close all

addpath(genpath('./util'));
addpath('./plot_scripts');
addpath('../mex');
load('../data/2D/gazebo1.mat');

% test points (regular samples for visualization)
xmin = -5;
xmax = 20;
ymin = -15;
ymax = 5;
test_intv = 0.1;
[xg, yg] = meshgrid((xmin+test_intv):test_intv:(xmax-test_intv), (ymin+test_intv):test_intv:(ymax-test_intv));
xtest = single([xg(:)'; yg(:)']);

skip = 10;
initframe = 101; % first 100 frames are almost static... 
lastframe = (floor((size(poses,1)-initframe)/skip))*skip+initframe;
for nframe = initframe:skip:lastframe

    % pose
    tr = poses(nframe,1:2)';
    phi = poses(nframe,3);
    Rot = [cos(phi) -sin(phi); sin(phi) cos(phi)];

    % update
    mexGPisMap('update',single(thetas),...
                        single(ranges(nframe,:)'),...
                        single([tr; reshape(Rot,[],1)]));
  
    % test visualization
    if  1 % nframe == lastframe % set the condition to 1 to visualize every update 
        visualize_gpisMap        
    end
    
    % % pause if needed
    % disp('Press a button to continue')
    % pause
    
end

% delete all resources
mexGPisMap('reset');


