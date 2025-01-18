clc;clear;close all;

%% generate 3D synthetic data
a1=zeros(300,20);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
  t1(i)=round(140);
  t3(i)=round(-6*i+180);
  t4(i)=round(6*i+10);
  a1(t1(i):t1(i)+k-1,i)=b1; 
  a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b1;
end

temp=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
for j=1:20
    a4=zeros(300,20);
    for i=1:m
  t4(i)=round(6*i+10+3*j); 
  a4(t4(i):t4(i)+k-1,i)=b1;
  
  t1(i)=round(140-2*j);
  a1(t1(i):t1(i)+k-1,i)=b1;
    end
    shot(:,:,j)=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
end
plane3d=shot;
d=plane3d/max(max(max(plane3d)));

% without noise
dn=d;
size(d)  %300 20 20
% decimate
[nt,nx,ny]=size(d);
ratio=0.5;
mask=genmask(reshape(d,nt,nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;

% adding noise
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));
d1=dn.*mask;
figure;imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9)]);
d = permute(dn,[3,2,1]);
in = d(:,1:19,110);
%% TNN
T = permute(dn,[3,2,1]);
normalize = max(T(:));
Xn = T/normalize;
[n1,n2,n3] = size(Xn);
p = 0.5;
Omega = zeros(size(Xn));
chosen = randperm(n1*n2*n3,round(p*n1*n2*n3));
Omega(chosen) = 1;
Omega = Xn ~= 0;

alpha  = 1;
maxItr = 1000; % maximum iteration
rho = 0.08                          ;
myNorm = 'tSVD_1'; % dont change for now

A = diag(sparse(double(Omega(:)))); % sampling operator
b  = A * Xn(:); % available data
bb = reshape(b,[n1,n2,n3]);
X = tensor_cpl_admm(A,b,rho,alpha,[n1,n2,n3],maxItr,myNorm,0);
X = X * normalize                 ;
X = reshape(X,[n1,n2,n3]);    

figure;imagesc(Xn(:,:,150));
figure;imagesc(X(:,:,150));