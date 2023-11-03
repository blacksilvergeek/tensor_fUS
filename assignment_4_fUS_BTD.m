
% Clear previous work
close all; 
clear; 
clc; 
addpath(genpath('data')); 
addpath(genpath('given'));
%%
%%%%%%%%%%%%%%%%% Load experimental data %%%%%%%%%%%%%%%%%

load('data/params.mat'); 
x_axis = params.x_axis; % Pixel locations along width [mm]
z_axis = params.z_axis; % Pixel locations along depth [mm]
Fs = params.Fs; % Sampling rate of the experiment

% Load the power-Doppler images
load('data/pdi.mat'); 
% PDI = PDI;

% Load the binary stimulus vector
load('data/stim.mat');
% stim = stim;

Nz = size(PDI,1); % Number of pixels along depth
Nx = size(PDI,2); % Number of width pixels
Nt = size(PDI,3); % Number of timestamps
t_axis = 0:1/Fs:(Nt-1)/Fs; % Time-axis of the experiment
%%
%%%%%%%%%%%%%%%%%% Get to know the data %%%%%%%%%%%%%%%%%%

figure; imagesc(x_axis,z_axis,PDI(:,:,100)); 
ylabel('Depth [mm]'); xlabel('Width [mm]');
title(['PDI at ' num2str(round(100/Fs)) 'st second']);

figure; plot(t_axis,squeeze(PDI(10,10,:))); 
hold on; s = plot(t_axis,stim*3*10^6); % the stimulus is
  % multiplied with a factor to visualize it at the same 
                   % y-axis scale as the fUS time-series
xlabel('Time (s)'); ylabel('Power Doppler Amplitude');
title('Time-Series of the Pixel at (x,z)=(10,10)'); 
legend(s,'Stimulus');

% Calculate the mean PDI
mean_PDI = mean(abs(PDI), 3);
mean_PDI = mean_PDI./(max(mean_PDI(:)));
% Display log of mean_PDI to enhance the contrast
figure; imagesc(x_axis,z_axis,log(mean_PDI));  
title('Mean PDI');
ylabel('Depth [mm]'); xlabel('Width [mm]');
%%
%%%%%%%%%%%%%%%%%%%%% Pre-processing %%%%%%%%%%%%%%%%%%%%%

% Standardize the PDI along time
P = (PDI - mean(PDI, 3))./std(PDI,[],3); 

% Spatial Gaussian smoothing
ht = fspecial('gaussian',[4 4],2);
Pg = double(convn(P,ht,'same'));

% Temporal low pass filter at 0.3 Hz per pixel time-series
f1 = 0.3;
[b, a] = butter(5,f1/(Fs/2),'low');
PDIlinear = reshape(Pg,Nz*Nx,Nt);
Pgf = reshape(filtfilt(b,a,PDIlinear')',size(PDI));
PDI = Pgf;
%%
%%%%%%%%%%%%%%% Show the correlation image %%%%%%%%%%%%%%%
delay = 0:10;
delay_num = round(delay * Fs);
PDI_newlinear = reshape(PDI,Nz*Nx,Nt);
for i = 0:round(10*Fs)
    if i ~= 0
        %stim_new = [zeros(delay_num(i),1)',stim'];
        stim_new = [zeros(i,1)',stim'];
        stim_new = stim_new(1:length(stim));
    end
    if i == 0
        stim_new = stim';
    end
    PCC{i+1} = corr(PDI_newlinear',stim_new','Type','Pearson');
    PCC_aver(i+1) = mean(abs(PCC{i+1}));
end
delay_stim = find(PCC_aver==max(PCC_aver)); 
PCC_delay = PCC{delay_stim};
PCC_delay(abs(PCC_delay)<0.3)=0;
pc_image = reshape(abs(PCC_delay), [Nz Nx]);
delay = (delay_stim-1)/Fs
% Two ways to visualize the correlation image -->
plot_version = 1;
display_brain_img(pc_image,log(mean_PDI),z_axis,x_axis,...
    'Significantly Correlated Regions',plot_version);
%plot_version = 2;
%display_brain_img(pc_image,log(mean_PDI),z_axis,x_axis,...
 %   'Significantly Correlated Regions',plot_version);

%%%%%%%%%%%%%%%%%%%%%%%%%% BTD %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Fill in btd_ll1_als_3d.m
% Include plots for all your claims (you can use display_brain_img.m to 
% help with the visualization of spatial maps)
% Fixing a set of initial value

%% L = 3
% Fixing a set of initial value

%R2 = 12;
%L1= 3;
%A_ini = randn(size(PDI, 1), R2*L1);
%B_ini = randn(size(PDI, 2), R2*L1);
%C_ini = randn(size(PDI, 3), R2);
%Init_1{1} = A_ini;
%Init_1{2} = B_ini;
%Init_1{3} = C_ini;
%save('Init1.mat', 'Init_1');

%% BTD when L = 3
%Init_1 = load('Init1');
%Init_1 = Init_1.Init_1;
%R2 = 12;
%L1 = 3;
%options.maxiter = 30; 
%options.th_relerr = 0.6;
%[A_1, B_1, C_1, const, output] = btd_ll1_als_3d(PDI, R2, L1, options,Init_1);
%% Plot the Spatial maps
%stim_new = [zeros(delay_stim,1)',stim'];
%stim_new = stim_new(1:length(stim));
%figure;
%for i = 1:R2
%    sp_map = A_1(:,L1*(i-1)+1:L1*i)*B_1(:,L1*(i-1)+1:L1*i)';
%    subplot(4,3,i);
%    imagesc(sp_map);
%    title(['corr(',num2str(i),')=',num2str(corr(C_1(:,i), stim_new', 'Type','Pearson'))]);
%    colorbar;
%    c = colorbar;
%    c.Label.String = 'Amplitude';
%end
%sgtitle('Spatial Maps of BTD when L = 3');

%% Plot the BTD Spatial maps when L = 3
%figure
%corb3 = zeros(1, R2);
%for i = 1:R2
%    sc1 = A_1(:,L1*(i-1)+1:L1*i)*B_1(:,L1*(i-1)+1:L1*i)';
%    corb3(i) = corr(C_1(:, i), stim_new', 'Type', 'Pearson');  
%    ax_ = subplot(4, 3, i);
%     %imagesc(sc1);
%    title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
%    display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
%    size(sc1);
%    title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
%    colorbar;
%end
%% L = 2 
% Fixing a set of initial value
R2 = 12;
L2= 2;
A_ini = randn(size(PDI, 1), R2*L2);
B_ini = randn(size(PDI, 2), R2*L2);
C_ini = randn(size(PDI, 3), R2);
Init_2{1} = A_ini;
Init_2{2} = B_ini;
Init_2{3} = C_ini;
%save('Init2.mat', 'Init_2');
%% BTD when L = 2
R2 = 12;
L2 = 2;
Init_2 = load('Init2');
Init_2 = Init_2.Init_2;
options.maxiter = 300; 
options.th_relerr = 0.6;
[A_2, B_2, C_2, const, output] = btd_ll1_als_3d(PDI, R2, L2, options,Init_2);
%% Plot the Spatial maps when L = 2
figure
corb3 = zeros(1, R2);
for i = 1:R2
    sc1 = A_2(:,L2*(i-1)+1:L2*i)*B_2(:,L2*(i-1)+1:L2*i)';
    corb3(i) = corr(C_2(:, i), stim_new', 'Type', 'Pearson');  
    ax_ = subplot(5, 4, i);
    % imagesc(sc1);
    title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
    display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
    size(sc1);
    title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
    colorbar;
end
%% figure
%n=[1,7,8];
%corb3 = zeros(1, R2);
 %   i = 1;
 %   sc1 = A_2(:,L2*(i-1)+1:L2*i)*B_2(:,L2*(i-1)+1:L2*i)';
 %   corb3(i) = corr(C_2(:, i), stim_new', 'Type', 'Pearson');  
 %   ax_ = subplot(1, 3, 1);
 %   % imagesc(sc1);
 %   title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
 %   display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
 %   size(sc1);
 %   title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
 %   colorbar;
 %   i = 7
 %   sc1 = A_2(:,L2*(i-1)+1:L2*i)*B_2(:,L2*(i-1)+1:L2*i)';
 %   corb3(i) = corr(C_2(:, i), stim_new', 'Type', 'Pearson');  
 %   ax_ = subplot(1, 3, 2);
 %   % imagesc(sc1);
 %   title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
 %   display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
 %   size(sc1);
 %   title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
 %   colorbar;
 %   i = 8
 %   sc1 = A_2(:,L2*(i-1)+1:L2*i)*B_2(:,L2*(i-1)+1:L2*i)';
 %   corb3(i) = corr(C_2(:, i), stim_new', 'Type', 'Pearson');  
 %   ax_ = subplot(1, 3, 3);
 %   % imagesc(sc1);
 %   title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
 %   display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
 %   size(sc1);
 %   title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
 %   colorbar;

%% Plot the Spatial maps
%figure;
%for i = 1:R2
 %   sp_map = A_2(:,L2*(i-1)+1:L2*i)*B_2(:,L2*(i-1)+1:L2*i)';
  %  subplot(4,3,i);
   % imagesc(sp_map);
    %title(['corr(',num2str(i),')=',num2str(corr(C_2(:,i), stim_new', 'Type','Pearson'))]);
    %colorbar;
    %c = colorbar;
    %c.Label.String = 'Amplitude';
%end
%sgtitle('Spatial Maps of BTD when L = 2');





%% Function
function corb3 = plot_spatial_maps_cpd(B1, B2, B3, stim, opt_lag)  % use opt_lag instead of lagindex
    % load mean_PDI , z_axis, x_axis, plot_version from workspace
    mean_PDI = evalin('base', 'mean_PDI');
    z_axis = evalin('base', 'z_axis');
    x_axis = evalin('base', 'x_axis');
    plot_version = evalin('base', 'plot_version');

    figure;
    corb3 = zeros(1, size(B1, 2));
    for i = 1:size(B1, 2)
        sc1 = B1(:, i) * B2(:, i)';
        corb3(i) = corr(B3(:, i), circshift(stim, opt_lag), 'Type', 'Pearson');  % use opt_lag instead of lagindex
        subplot(5, 3, i);
        % imagesc(sc1);
        display_brain_img(sc1,log(mean_PDI),z_axis,x_axis,...
    'Significantly Correlated Regions',plot_version);
        size(sc1);
        title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
        colorbar;
    end
    sgtitle('Spatial Maps of CP Decomposed Components');
end