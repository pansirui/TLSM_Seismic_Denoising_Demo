function [Denoising, k]  = TLSM_Denoising(ori, Y, par)

[Nx,Ny,Nt] = size(Y);

tau = par.tau;
lambda1 = par.lambda1;
lambda2 = par.lambda2;
alpha = par.alpha; 
beta = par.beta; 
gamma = par.gamma; 

Bx = zeros(Nx,Ny,Nt); 
By = zeros(Nx,Ny,Nt); 
B = zeros(Nx,Ny,Nt);
Dx = zeros(Nx,Ny,Nt); 
Dy = zeros(Nx,Ny,Nt); 
Z = Y; 

Denoising = cell(1, par.Iter);

for k = 1: par.Iter

    % Update X
    Xnew = UpdateXTen(Y,Z,Dx,Dy,B,Bx,By,alpha,beta,gamma);
    X = Xnew;    
    
    % Update Z
    [Znew, objV] = proxF_tSVD_LSM(X+B,tau/alpha,1);
    Z = Znew;

    % Update Dx
    Dxx = zeros(Nx,Ny,Nt);
    for i = 1:Nt 
        temp = diff(X(:,:,i)-Y(:,:,i),1,2);
        dx = [temp temp(:,Ny-1)];  
        Dxx(:,:,i) = dx;
    end
    Dx_new = proxF_GSM_sparsity_3D_W(Dxx+Bx,lambda1/beta);   
    Dx = Dx_new; 

    % Update Dy
    Dyy = zeros(Nx,Ny,Nt);
    for i = 1:Nt
        temp1 = diff(X(:,:,i),1,1);
        dy = [temp1;temp1(Nx-1,:)];
        Dyy(:,:,i) = dy;
    end
    Dy_new = proxF_GSM_sparsity_3D_W(Dyy+By,lambda2/gamma);   
    Dy = Dy_new;

    % Update B, Bx, and By
    B_new  = B+(X-Z);
    Bx_new = zeros(Nx,Ny,Nt);
    By_new = zeros(Nx,Ny,Nt);
    for i = 1:Nt
        Bx_new(:,:,i) = Bx(:,:,i) + (X(:,[2:Ny,Ny],i) - X(:,:,i)) - (Y(:,[2:Ny,Ny],i)- Y(:,:,i)) - Dx(:,:,i);
        By_new(:,:,i) = By(:,:,i) +   X([2:Nx,Nx],:,i)- X(:,:,i) - Dy(:,:,i);
    end
    B = B_new;    
    Bx = Bx_new;     
    By = By_new;
 
    Denoising{k} = Xnew;

end

end

