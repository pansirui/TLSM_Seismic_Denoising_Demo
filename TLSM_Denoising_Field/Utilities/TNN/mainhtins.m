clc;clear;


load('F:\wangyan\matlabpython_workspace\mat\htins_vsp.mat')
d=ori_;
d = d(1:218,1:31,1:205,1:205);
d = d(:,1,1:48,1:48);
d = permute(d,[4 3 2 1]);
% ori_ = d;
% save('clean_htins.mat','ori_');
d=reshape(d,48,48,218);

%% add noise
randn('state',201314);  % 加随机噪声
var=0.01;
d0=d+var*randn(size(d));
%% TNN
T = d0; % y x t
normalize = max(T(:));
Xn = T/normalize;   % 归一化
[n1,n2,n3] = size(Xn);
p = 0.5;
Omega = zeros(size(Xn));
chosen = randperm(n1*n2*n3,round(p*n1*n2*n3)); % 502272个数中随机选择251136个数[round(四舍五入)] 
Omega(chosen) = 1;  %随机取一半的数置1
Omega = Xn ~= 0;   % Omega 采样矩阵

alpha  = 1.5;
maxItr = 1000; % maximum iteration
rho = 0.08;
myNorm = 'tSVD_1'; % dont change for now

A = diag(sparse(double(Omega(:)))); % sampling operator
b  = A * Xn(:); % available data
bb = reshape(b,[n1,n2,n3]);
X = tensor_cpl_admm(A,b,rho,alpha,[n1,n2,n3],maxItr,myNorm,0);
X = X * normalize;
X = reshape(X,[n1,n2,n3]);    
recover = X;% x y t
%figure;imagesc(Xn(:,:,150));
% ori_=recover;
% save('recover.mat','ori_'); %x y t

d = reshape(d,48,48,1,218);
d0 = reshape(d0,48,48,1,218);
recover = reshape(recover,48,48,1,218);
figure;imagesc(reshape(d(:,:,1,200),[48,48]));
figure;imagesc(reshape(d0(:,:,1,200),[48,48]));
figure;imagesc(reshape(recover(:,:,1,200),[48,48]));
figure;imagesc(reshape(X(:,:,200),[48,48]));

A = d0;
B = d;
fprintf('***********************Quality factor = %d, RSE = %d ***********\n', -20*log10(norm(A(:) - B(:)) / norm(B(:))), norm(A(:) - B(:)) / norm(B(:)));



