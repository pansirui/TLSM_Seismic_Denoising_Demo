function xnew = Updatex(y,z,b,bx,by,dx,dy,alpha,beta,gamma)
    
    ep = 1e-6;
    [conjoDx,conjoDy,num1,Denom1,Denom2] = getC(y);
    
    FPsix = fft2( dx - bx );
    FPsiy = fft2( dy - by );
    FGPx = conjoDx .* FPsix; 
    FGPy = conjoDy .* FPsiy; 
    Denom = alpha + beta*Denom1 + gamma*Denom2 + fft2(eye(1));   
    
    Fu = alpha*fft2(z-b) + beta*num1 + beta* FGPx + gamma*FGPy + fft2(y); 
    
    xnew = real( ifft2(Fu./(Denom+ep)) ); 

end