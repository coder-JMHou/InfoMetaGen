close all;
clear;
clc;
iter_kkk=36;
total_time_s = zeros(iter_kkk,1);
psnr_s = zeros(iter_kkk,1);
ssim_s = zeros(iter_kkk,1);
for kkk=1:iter_kkk
tStart=tic;
lambda = 1;
num_units = 64;
d = lambda/3;
p = d;
z0 = 20;
dis = 0;

num_x = num_units;
num_y = num_units;

A = digitTrain4DArrayData;
B = A(:,:,:,1);
B = double(B);

% 归一化到 [0,1]
B = (B - min(B(:))) / (max(B(:)) - min(B(:)) + eps);

C = zeros(num_units, num_units);

[hB, wB] = size(B);
rowStart = floor((num_units - hB)/2) + 1;
colStart = floor((num_units - wB)/2) + 1;
rowEnd   = rowStart + hB - 1;
colEnd   = colStart + wB - 1;

C(rowStart:rowEnd, colStart:colEnd) = B;
known_abs_spatial = C;

horn_imamp_single = ones(num_x*num_y, 1);

[metrix, ~] = Creat_metrix(num_x, num_y, p, lambda, z0, 1, dis);

problem = build_hologram_problem( ...
    known_abs_spatial, metrix, ...
    'horn_imamp_single', horn_imamp_single, ...
    'optimizer', 'ga', ...
    'designMode', 'integer8', ...
    'quantLevel', 8, ...
    'beta', 0.25, ...
    'gamma', 0.00, ...
    'maskThreshold', 0.05);

optsGA = optimoptions('ga', ...
    'Display', 'iter', ...
    'PopulationSize', 10, ...      % 4096维下建议先小一点试
    'MaxGenerations', 1000, ...
    'UseParallel', false);

resultGA = optimize_hologram_builtin(problem, 'ga', optsGA);

figure;
subplot(1,3,1); imagesc(resultGA.code_2d); axis image; colorbar; title('GA编码矩阵');
subplot(1,3,2); imagesc(known_abs_spatial); axis image; colorbar; title('目标图像');
subplot(1,3,3); imagesc(resultGA.target_amp_norm_2d); axis image; colorbar; title('GA重建结果');

fprintf('GA: score=%.6f, RMSE=%.6f, eta=%.6f, PCC=%.6f\n', ...
    resultGA.score, resultGA.rmse, resultGA.eta_target, resultGA.pcc);
%%
load('xianzai1.mat')
% yuanlai=B;
xianzai=resultGA.target_amp_norm_2d;
% xianzai=mapminmax(xianzai,0,1);
xianzai2=xianzai((num_units-28)/2+1:(num_units-28)/2+27+1,(num_units-28)/2+1:(num_units-28)/2+27+1);
% save('xianzai1','xianzai1')
P=psnr(xianzai1,xianzai2);
S=ssim(xianzai1,xianzai2);
total_time_s(kkk,1) = toc(tStart);
psnr_s(kkk,1) = P;
ssim_s(kkk,1) = S;
end
%%
GA_matrix(:,1)=psnr_s;
GA_matrix(:,2)=ssim_s;
GA_matrix(:,3)=total_time_s;
save('GA.mat','GA_matrix')