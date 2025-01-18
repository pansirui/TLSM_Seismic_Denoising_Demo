clc;clear;close all;

load('C:\Users\WangYan\Desktop\MSSA\MSSA\Dip_steered_median_filter_gauss_02.mat'); %643 943 463 x y t
d = ori_(1:48,1:48,1:450);
ori_ = d;
save('clean.mat','ori_');
%figure;imagesc(reshape(d(:,50,:),[48,450]));
% decimate
load('C:\Users\WangYan\Desktop\MSSA\MSSA\mask_0.5.mat');
mask = ori_;
for i=1:450
    d0(:,:,i) = d(:,:,i) .* mask(1:48,1:48);
end
ori_=d0;
save('noise.mat','ori_'); %x y t
%figure;imagesc(reshape(d0(:,50,:),[100,450]));
%% TNN
T = permute(d0,[2 1 3]); % y x t
normalize = max(T(:));
Xn = T/normalize;
[n1,n2,n3] = size(Xn);
p = 0.5;
Omega = zeros(size(Xn));
chosen = randperm(n1*n2*n3,round(p*n1*n2*n3));
Omega(chosen) = 1;
Omega = Xn ~= 0;

alpha  = 1;
maxItr = 200; % maximum iteration
rho = 0.08;
myNorm = 'tSVD_1'; % dont change for now

A = diag(sparse(double(Omega(:)))); % sampling operator
b  = A * Xn(:); % available data
bb = reshape(b,[n1,n2,n3]);
X = tensor_cpl_admm(A,b,rho,alpha,[n1,n2,n3],maxItr,myNorm,0);
X = X * normalize;
X = reshape(X,[n1,n2,n3]);    
recover = permute(X,[2 1 3]);% x y t
%figure;imagesc(Xn(:,:,150));
ori_=recover;
save('recover.mat','ori_'); %x y t
%figure;imagesc(reshape(recover(:,50,:),[48,450]));