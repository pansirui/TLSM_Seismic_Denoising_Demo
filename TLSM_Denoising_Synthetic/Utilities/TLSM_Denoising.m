function [Denoising, k]  = TLSM_Denoising(ori, Y, par)

[Nx,Ny,Nt] = size(Y);

tau = par.tau;
lambda1 = par.lambda1;
lambda2 = par.lambda2;
alpha = par.alpha; 
beta = par.beta; 
gamma = par.gamma; 

% Initialization
Bx = zeros(Nx,Ny,Nt); 
By = zeros(Nx,Ny,Nt); 
B = zeros(Nx,Ny,Nt);
Dx = zeros(Nx,Ny,Nt); 
Dy = zeros(Nx,Ny,Nt); 
Z = Y; 

Denoising = cell(1, par.Iter);

% Iteration process
for k = 1: par.Iter

    % Update X
    Xnew = UpdateXTen(Y,Z,Dx,Dy,B,Bx,By,alpha,beta,gamma);
    X = Xnew;    
    
    % Update Z
    [Znew, objV] = proxF_tSVD_LSM(X+B,tau/alpha,1);
    Z = Znew;

    % Update D1
    Dxx = zeros(Nx,Ny,Nt);
    for i = 1:Nt 
        temp = diff(X(:,:,i)-Y(:,:,i),1,2);
        dx = [temp temp(:,Ny-1)];  
        Dxx(:,:,i) = dx;
    end
    Dx_new = proxF_LSM_sparsity_3D_W(Dxx+Bx,lambda1/beta);   
    Dx = Dx_new;

    % Update D2
    Dyy = zeros(Nx,Ny,Nt);
    for i = 1:Nt
        temp1 = diff(X(:,:,i),1,1);
        dy = [temp1;temp1(Nx-1,:)];
        Dyy(:,:,i) = dy;
    end
    Dy_new = proxF_LSM_sparsity_3D_W(Dyy+By,lambda2/gamma);      
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
 
    Denoising{k} = Xnew;

end

    % results of seismic data noise suppression
    recover = Denoising{k};
    removed_noise = Y - recover;
    
    % plot
    % inline, crossline
    s_cplot(squeeze(Y(:,300,:))'); title('Noisy data'); xlabel('Crossline'); ylabel('Inline');
    s_cplot(squeeze(recover(:,300,:))'); title('TLSM'); xlabel('Crossline'); ylabel('Inline');
    s_cplot(squeeze(removed_noise(:,300,:))'); title('Removed noise'); xlabel('Crossline'); ylabel('Inline');

    % time, crossline
    s_cplot(squeeze(Y(:,:,30))'); title('Noisy data'); xlabel('Crossline'); ylabel('Time');
    s_cplot(squeeze(recover(:,:,30))'); title('TLSM'); xlabel('Crossline'); ylabel('Time');
    s_cplot(squeeze(removed_noise(:,:,30))'); title('Removed noise'); xlabel('Crossline'); ylabel('Time');

    % time, inline
    s_cplot(squeeze(Y(10,:,:))); title('Noisy data'); xlabel('Inline'); ylabel('Time');
    s_cplot(squeeze(recover(10,:,:))); title('TLSM'); xlabel('Inline'); ylabel('Time');
    s_cplot(squeeze(removed_noise(10,:,:))); title('Removed noise'); xlabel('Inline'); ylabel('Time');

end

