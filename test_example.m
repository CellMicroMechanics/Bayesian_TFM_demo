%                        Introduction
% This script runs an example to demonstrate the Bayesian Traction Force Microscopy
% code BL2 described in Huang et al. https://arxiv.org/abs/1810.05848
%
% This artificial test data involves 3 circular traction spots with having constant
% traction magnitudes in the range of [0-250] Pa. The test displacements used for traction 
% reconstruction are perturbed with Gaussian noise having a variance 0.0025
% pix^2. Using these displacements, the traction data is reconstructed.
%
% This program requires the free Matlab package "Regularization tools" by P.C.
% Hansen Math. Model. Comput.(1998) that can be found at
% https://www.mathworks.com/matlabcentral/fileexchange/52-regtools.
% To use this tool, download the package and add it to the Matlab path
% ,e.g, by entering the command "path(path,'Hansen')".
%
% Yunfei Huang and Benedikt Sabass, Forschungszentrum J�lich, 2018
% Reference: Huang et al. https://arxiv.org/abs/1810.05848.
tic 
clear all; 
disp('Running Bayesian TFM demo. The calculations can take a few seconds.');

load('test_X.mat');
load('test_displacement.mat');
load('test_noise.mat');
load('test_traction.mat')

u(1:2:size(displacement.vec,1)*2,1) = displacement.vec(:,1);
u(2:2:size(displacement.vec,1)*2,1) = displacement.vec(:,2);

noise_u(1:2:size(noise,1)*2,1) = noise(:,1);
noise_u(2:2:size(noise,1)*2,1) = noise(:,2);
beta = 400; % For the test data, the variance of noise is 0.0025. 

F = reconstruction_BL2(X, u,noise_u,beta);

f(1:length(F)/2,1) = F(1:2:end);
f(1:length(F)/2,2) = F(2:2:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot reconstructed traction and original traction for comparison  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for n=1:21*21
 ft(n,1) = sqrt(f(n,1)^2+f(n,2)^2);
 end
 x = reshape(test_traction.r_pos(:,1),[21,21]);
 y = reshape(test_traction.r_pos(:,2),[21,21]);
 ff = reshape(ft(:,1),[21,21]);
   
figure(2);
surf(x,y,ff,'LineStyle','none','FaceColor','interp');
colorbar;
view(2);
hold on;
quiver(test_traction.r_pos(:,1),test_traction.r_pos(:,2),f(:,1),f(:,2),'r');
view(90,-90);
colormap('jet');
pbaspect([1 1 1]);
caxis([0 250]);
title('Reconstruction using BL2');
hold on;
for k = 1:3
radius=20;  
pos = [test_traction.center(k,1)-radius test_traction.center(k,2)-radius...
       2*radius 2*radius];
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','r','linewidth',2)
hold on; 
end
hold off;
 

figure(3);
surf(test_traction.pos_x, test_traction.pos_y, test_traction.vec, 'LineStyle', 'none', 'FaceColor', 'interp');
colorbar;
caxis([0 250]);
view(2);
view(90,-90);
pbaspect([1 1 1]);
colormap('jet');
title('Artificial input traction');
hold on;
h1 = quiver(test_traction.center(:,1),test_traction.center(:,2),...
     test_traction.direction(:,1),test_traction.direction(:,2),'r');
set(h1,'AutoScale','on', 'AutoScaleFactor', 0.3,'linewidth',2);
 
 toc
