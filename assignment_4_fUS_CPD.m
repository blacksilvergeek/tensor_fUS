
% Clear previous work
close all; 
clear; 
clc; 

addpath(genpath('data')); 
addpath(genpath('given')); 

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

%%%%%%%%%%%%%%% Show the correlation image %%%%%%%%%%%%%%%

% pc_image = ...

max_lag = 10; % maximum lag in seconds
step = 1/Fs; % step for lag in seconds
max_corr = -Inf;
opt_lag = 0;

filename = 'opt_lag.mat';

if exist(filename, 'file')
    load(filename);
    disp('opt_lag loaded successfully');
else
    for lag = 0:step:max_lag
        shift_size = round(lag * Fs);
        shifted_stim = [zeros(1, shift_size), stim(1:end-shift_size)']';

        PDIlinear = reshape(PDI, Nz*Nx, Nt);
        corr_values = corr(PDIlinear', shifted_stim);
        avg_corr = mean(abs(corr_values));
        if avg_corr > max_corr
            max_corr = avg_corr;
            opt_lag = lag;
        end
    end
end


% save and load
save('opt_lag.mat', 'opt_lag');



% Calculate the correlation image at the optimal lag
shifted_stim = circshift(stim, round(opt_lag * Fs));
PDIlinear = reshape(PDI, Nz*Nx, Nt);
corr_values = corr(PDIlinear', shifted_stim);
corr_image = reshape(corr_values, Nz, Nx);

% Find significantly correlated pixels
% p_values = 1 - normcdf(abs(corr_values), 0, 1);
% p_sig = p_values < 0.3;

% to observe same figure in assignment
% threshold = 0.36;

threshold = 0.3;
p_sig = (corr_values) > threshold;

pc_vals = corr_values.*p_sig;

pc_image = reshape(corr_values.*p_sig, [Nz Nx]);

% Display the correlation image at the optimal lag
plot_version = 1;
display_brain_img(pc_image, log(mean_PDI), z_axis, x_axis, ...
    ['Significantly Correlated Regions at lag ' num2str(opt_lag) 's'], plot_version);


% Two ways to visualize the correlation image -->
plot_version = 1;
display_brain_img(pc_image,log(mean_PDI),z_axis,x_axis,...
    'Significantly Correlated Regions',plot_version);
% plot_version = 2;
% display_brain_img(pc_image,log(mean_PDI),z_axis,x_axis,...
%      'Significantly Correlated Regions',plot_version);

disp("First part done")

%%%%%%%%%%%%%%%%%%%%%%%%%% CPD %%%%%%%%%%%%%%%%%%%%%%%%%%%

% You may use hidden_cpd_als_3d.m
% Include plots for all your claims (you can use display_brain_img.m to 
% help with the visualization of spatial maps)

% Set parameters
R1 = 15;
options.maxiter = 500; 
options.th_relerr = 0.6;

filename = 'cpd_save.mat';


if_new_CPD = false;
if exist(filename, 'file') && ~if_new_CPD
    load(filename);
    B1 = cpd_save{1};
    B2 = cpd_save{2};
    B3 = cpd_save{3};
    c = cpd_save{4};
    output = cpd_save{5};
    disp('CPD results loaded successfully');
else
    disp('File does not exist. Performing CPD...');
    [B1, B2, B3, c, output] = cpd_als_3d(PDI, R1, options);
    
    % Save CPD results
    cpd_save = {B1, B2, B3, c, output};
    save(filename, 'cpd_save');
    disp('CPD results calculated and saved');
end


% Plot spatial maps of CP decomposed components
% corb = plot_spatial_maps_cpd(B1, B2, B3, stim, round(opt_lag*Fs));

n_shift = round(opt_lag*Fs);

figure
corb3 = zeros(1, size(B1, 2));
for i = 1:size(B1, 2)
    sc1 = B1(:, i) * B2(:, i)';
    corb3(i) = corr(B3(:, i), circshift(stim, n_shift), 'Type', 'Pearson');  
    ax_ = subplot(4, 5, i);
    % imagesc(sc1);
    title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
    display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
    size(sc1);
    title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
    colorbar;
end
% sgtitle('Spatial Maps of CP Decomposed Components');

% Set figure to full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, 'imgs/cpd_all.png');


    % Plot relative error per iteration
    plot_relative_error(output.relerr, 'CP decomposition');
    
    % Display strongest correlations
    disp_strongest_correlations(corb3, 0.3);
    
    
    % Compare spatial maps with correlation
    select_indices = [6, 9, 10];
    compare_spatial_maps_with_correlation(B1, B2, corb3, select_indices, pc_image);
    
    disp("Second part done")


%%% R1 = 12
% For R1 = 12
R1 = 12;
options.maxiter = 500; 
options.th_relerr = 0.6;

filename = 'cpd_save_R12.mat';
if_new_CPD = true;
if exist(filename, 'file') && ~if_new_CPD
    load(filename);
    B1 = cpd_save{1};
    B2 = cpd_save{2};
    B3 = cpd_save{3};
    c = cpd_save{4};
    output = cpd_save{5};
    disp('CPD results loaded successfully');
else
    disp('File does not exist. Performing CPD...');
    [B1, B2, B3, c, output] = cpd_als_3d(PDI, R1, options);
    
    % Save CPD results
    cpd_save = {B1, B2, B3, c, output};
    save(filename, 'cpd_save');
    disp('CPD results calculated and saved');
end

% Plot spatial maps of CP decomposed components
% corb = plot_spatial_maps_cpd(B1, B2, B3, stim, round(opt_lag*Fs));

n_shift = round(opt_lag*Fs);

figure
corb3 = zeros(1, size(B1, 2));
for i = 1:size(B1, 2)
    sc1 = B1(:, i) * B2(:, i)';
    corb3(i) = corr(B3(:, i), circshift(stim, n_shift), 'Type', 'Pearson');  
    ax_ = subplot(3, 4, i);
    % imagesc(sc1);
    title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
    display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
    size(sc1);
    title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
    colorbar;
end
% sgtitle('Spatial Maps of CP Decomposed Components');

% Set figure to full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, 'imgs/cpd_r12.png');



%%% R1 = 20
% For R1 = 20
R1 = 20;
options.maxiter = 500; 
options.th_relerr = 0.6;

filename = 'cpd_save_R20.mat';
if_new_CPD = true;
if exist(filename, 'file') && ~if_new_CPD
    load(filename);
    B1 = cpd_save{1};
    B2 = cpd_save{2};
    B3 = cpd_save{3};
    c = cpd_save{4};
    output = cpd_save{5};
    disp('CPD results loaded successfully');
else
    disp('File does not exist. Performing CPD...');
    [B1, B2, B3, c, output] = cpd_als_3d(PDI, R1, options);
    
    % Save CPD results
    cpd_save = {B1, B2, B3, c, output};
    save(filename, 'cpd_save');
    disp('CPD results calculated and saved');
end

% Plot spatial maps of CP decomposed components
% corb = plot_spatial_maps_cpd(B1, B2, B3, stim, round(opt_lag*Fs));

n_shift = round(opt_lag*Fs);

figure
corb3 = zeros(1, size(B1, 2));
for i = 1:size(B1, 2)
    sc1 = B1(:, i) * B2(:, i)';
    corb3(i) = corr(B3(:, i), circshift(stim, n_shift), 'Type', 'Pearson');  
    ax_ = subplot(4, 5, i);
    % imagesc(sc1);
    title_string = ['corr(', num2str(i), ')=', num2str(corb3(i))];
    display_brain_img_sub((sc1),log(mean_PDI),z_axis,x_axis,title_string,plot_version, ax_);
    size(sc1);
    title(['corr(', num2str(i), ')=', num2str(corb3(i))]);
    colorbar;
end
% sgtitle('Spatial Maps of CP Decomposed Components');

% Set figure to full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, 'imgs/cpd_r20.png');






%%%%%%%%%%%%%%%%%%%%%%%%%% BTD %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fill in btd_ll1_als_3d.m
% Include plots for all your claims (you can use display_brain_img.m to 
% help with the visualization of spatial maps)




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

function plot_relative_error(relerr, title_str)
    figure; 
    semilogy(relerr)
    grid on
    xlabel('Iteration number')
    ylabel('Relative error $\frac{\| T - T_{dec} \|_F}{\| T \|_F}$', 'interpreter', 'latex');
    title(title_str);
end

function disp_strongest_correlations(corb3, threshold)
    strongest_corr = find(abs(corb3) > threshold);
    disp(['The strongest correlation appears at ' num2str(strongest_corr)]);
end

function compare_spatial_maps_with_correlation(B1, B2, corb3, select_indices, pc_image)
    figure;
    subplot(2, 2, 1);
    imagesc(pc_image);
    title('(a) correlation image');
    colorbar;
    for i = 1:length(select_indices)
        idx = select_indices(i);
        subplot(2, 2, i + 1);
        imagesc(B1(:, idx) * B2(:, idx)');
        title({['(b) component ', num2str(idx)]; ['corr=' num2str(corb3(idx))]});
        colorbar;
    end
    sgtitle({'Spatial Maps'; '- comparison between correlation and CPD component'});
end