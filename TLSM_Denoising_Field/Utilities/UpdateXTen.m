function Xnew = UpdateXTen(Y,Z,Dx,Dy,B,Bx,By,alpha,beta,gamma)

    Xnew = zeros(size(Y));
    
    for i = 1:size(Y,3)
        
        y = Y(:,:,i); 
        z = Z(:,:,i);
        
        dx = Dx(:,:,i); 
        
        dy = Dy(:,:,i);
        
        b  = B(:,:,i); 
        
        bx = Bx(:,:,i); 
        
        by = By(:,:,i);
        
        xnew = Updatex(y,z,b,bx,by,dx,dy,alpha,beta,gamma);
        
        Xnew(:,:,i) = xnew;
        
    end

end