
% Clear previous work
close all; 
clear; 
clc; 
addpath(genpath('data')); 
addpath(genpath('given'));
addpath(genpath('tensorlab'))
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

%% ll1 function 
nn = 10;
L = [1,1,1,1,1,...
    2,2,2,2,2];
% L = ones(1,nn-1);
% % L = [L*2 4];
% 
% L = ones(1,nn/2);
% L = [L ones(1,nn/2)*2];

L = ones(1,nn-2);
L = [L 1 1];

L = [ones(1,nn-2)...
    ones(1,1)*10 ones(1,1)*5];
% L = [1,1,1,1,...
%      2,2,3,3];

% L = ones(1,nn)*2;

refinementOptions = struct;
refinementOptions.MaxIter = 300;
init = @ll1_gevd;
method_main = @ll1_minf;
compression_method = @mlsvd;
Refinement = false;

corb3 = zeros(1, nn);
close_figure = false;


% while last 2 smaller than 0.3
% while abs(corb3(end)) < 0.25 && abs(corb3(end-1)) < 0.25
% while abs(corb3(end)) < 0.25
for i = 1
    % close last 2 figures
    if close_figure
        close(gcf);
        close(gcf);
    end
    close_figure = true;
    
   
    Uhat = ll1(PDI, L,'Display', 1, 'Initialization', ...
                init, 'Algorithm', method_main, ...
                'Compression',compression_method);
    % Uhat = ll1_gevd(PDI, L,'Display', 1, 'Initialization', init, 'Algorithm', method, 'Compression',compression_method);

    %%
    figure

    corb3 = zeros(1, nn);
    for i = 1:nn
        sc1 = normc(Uhat{i}{1})*normc(Uhat{i}{2})';
        corb3(i) = corr(Uhat{i}{3}, stim_new', 'Type', 'Pearson');  
        ax_ = subplot(3, 4, i);
        %imagesc(sc1);
        title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
        display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
        size(sc1);
        title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
        colorbar;
    end

    % full screen of the figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    figure
    corb3 = zeros(1, nn);
    for i = 1:nn
        sc1 = normc(Uhat{i}{1})*normc(Uhat{i}{2})';
        corb3(i) = corr(Uhat{i}{3}, stim_new', 'Type', 'Pearson');  
        ax_ = subplot(5, 6, i);
        imagesc(sc1);
        % title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
        % display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
        size(sc1);
        title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
        colorbar;
    end
    % full screen of the figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

    % save Uhat.mat Uhat
    save('Uhat.mat', 'Uhat');

    
end


% T2 = ll1gen(Uhat);
% nn2 = 5;
% L = ones(1,nn-2);
% L = [L 1 1];
% Uhat2 = ll1(T2, L,'Refinement', @ll1, 'RefinementOptions', refinementOptions,'Display', 1);
% %%
% figure

% corb3 = zeros(1, nn);
% for i = 1:nn2
%     sc1 = normc(Uhat2{i}{1})*normc(Uhat2{i}{2})';
%     corb3(i) = corr(Uhat2{i}{3}, stim_new', 'Type', 'Pearson');  
%     ax_ = subplot(5, 6, i);
%     %imagesc(sc1);
%     title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
%     display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
%     size(sc1);
%     title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
%     colorbar;
% end

% % full screen of the figure
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% figure
% corb3 = zeros(1, nn);
% for i = 1:nn2
%     sc1 = normc(Uhat2{i}{1})*normc(Uhat2{i}{2})';
%     corb3(i) = corr(Uhat2{i}{3}, stim_new', 'Type', 'Pearson');  
%     ax_ = subplot(5, 6, i);
%     imagesc(sc1);
%     % title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
%     % display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
%     size(sc1);
%     title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
%     colorbar;
% end
% % full screen of the figure
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% % save Uhat2.mat Uhat2
% save('Uhat2.mat', 'Uhat2');






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