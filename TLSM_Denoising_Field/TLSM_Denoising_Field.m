clear all
close all
clc

% noisy seismic data
Y = load('peno_250_150_40.mat'); 
Y = Y.peno_250_150_40;
Y = permute(Y, [2 1 3]);

% normalization
Y = ((Y - min(Y,[],'all')) / (max(Y,[],'all') - min(Y,[],'all')));
Y = 2 * Y - 1;
 
randn ('seed',0);

% parameters setting
tau = 0.1;
alpha = 1;
beta = 1; 
gamma = 0.05;
lambda1 = 10;
lambda2 = 1;

[Nx,Ny,Nt] = size(Y);

% Initialization
Bx = zeros(Nx,Ny,Nt); 
By = zeros(Nx,Ny,Nt); 
B = zeros(Nx,Ny,Nt);
Dx = zeros(Nx,Ny,Nt); 
Dy = zeros(Nx,Ny,Nt); 
Z = Y; 

Iter = 20;
% Iteration process
for k = 1: Iter

    % Update X
    Xnew = UpdateXTen(Y, Z, Dx, Dy, B, Bx, By, alpha, beta, gamma);
    X = Xnew;    
    
    % Update Z
    [Znew, objV] = proxF_tSVD_LSM(X+B, tau/alpha, 1);
    Z = Znew;

    % Update D1
    Dxx = zeros(Nx,Ny,Nt);
    for i = 1:Nt 
        temp = diff(X(:,:,i)-Y(:,:,i),1,2);
        dx = [temp temp(:,Ny-1)];  
        Dxx(:,:,i) = dx;
    end
    Dx_new = proxF_LSM_sparsity_3D_W(Dxx+Bx, lambda1/beta);   
    Dx = Dx_new; 

    % Update D2
    Dyy = zeros(Nx,Ny,Nt);
    for i = 1:Nt
        temp1 = diff(X(:,:,i),1,1);
        dy = [temp1;temp1(Nx-1,:)];
        Dyy(:,:,i) = dy;
    end
    Dy_new = proxF_LSM_sparsity_3D_W(Dyy+By, lambda2/gamma);    
    Dy = Dy_new;

    % Update B, Bx, and By
    B_new  = B+(X-Z);
    Bx_new = zeros(Nx,Ny,Nt);
    By_new = zeros(Nx,Ny,Nt);
    for i = 1:Nt
        Bx_new(:,:,i) = Bx(:,:,i) + (X(:,[2:Ny,Ny],i) - X(:,:,i)) - (Y(:,[2:Ny,Ny],i)- Y(:,:,i)) - Dx(:,:,i);
        By_new(:,:,i) = By(:,:,i) + X([2:Nx,Nx],:,i)- X(:,:,i) - Dy(:,:,i);
    end
    B = B_new;    
    Bx = Bx_new;     
    By = By_new;

end

% Results of seismic data noise suppression
recover = X; 
removed_noise = Y - X;

% plot
% inline, crossline
s_cplot(squeeze(Y(:,:,30))); title('Noisy data'); xlabel('Inline'); ylabel('Crossline');
s_cplot(squeeze(recover(:,:,30))); title('TLSM'); xlabel('Inline'); ylabel('Crossline');
s_cplot(squeeze(removed_noise(:,:,30))); title('Removed noise'); xlabel('Inline'); ylabel('Crossline');

% time, crossline
s_cplot(squeeze(Y(:,100,:))'); title('Noisy data'); xlabel('Crossline'); ylabel('Time');
s_cplot(squeeze(recover(:,100,:))'); title('TLSM'); xlabel('Crossline'); ylabel('Time');
s_cplot(squeeze(removed_noise(:,100,:))'); title('Removed noise'); xlabel('Crossline'); ylabel('Time');

% time, inline
s_cplot(squeeze(Y(60,:,:))'); title('Noisy data'); xlabel('Inline'); ylabel('Time');
s_cplot(squeeze(recover(60,:,:))'); title('TLSM'); xlabel('Inline'); ylabel('Time');
s_cplot(squeeze(removed_noise(60,:,:))'); title('Removed noise'); xlabel('Inline'); ylabel('Time');


% No-reference quality metrics
[nx, ny, nt] = size(recover);

BRISQUE_ALL = zeros(1,nt);
NIQE_ALL = zeros(1,nt);
PIQE_ALL = zeros(1,nt);

for index = 1:1:nt

    img = squeeze(recover(:,:,index));

    min_value = min(img(:));
    max_value = max(img(:));
    
    img = uint8(255 * (img - min_value) / (max_value - min_value));
    
    fprintf('Time slice: %d\n', index)

    brisqueScore = brisque(img);
    disp(['BRISQUE Score: ', num2str(brisqueScore)]);
    
    niqeScore = niqe(img);
    disp(['NIQE Score: ', num2str(niqeScore)]);
    
    piqeScore = piqe(img);
    disp(['PIQE Score: ', num2str(piqeScore)]);

    BRISQUE_ALL(index) = brisqueScore;
    NIQE_ALL(index) = niqeScore;
    PIQE_ALL(index) = piqeScore;

end

allScores = [BRISQUE_ALL(:), NIQE_ALL(:), PIQE_ALL(:)];
T = table(allScores(:,1), allScores(:,2), allScores(:,3),...
          'VariableNames', {'BRISQUE', 'NIQE', 'PIQE'});

excelFileName = 'scores_TLSM_Peno.xlsx';
writetable(T, excelFileName);
