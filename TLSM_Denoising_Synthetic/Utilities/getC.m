function [conjoDx,conjoDy,num1,Denom1,Denom2] = getC(g)

sizeI = size(g);

otfDx = psf2otf([1,-1],sizeI);
otfDy = psf2otf([1;-1],sizeI);
conjoDx = conj(otfDx);
conjoDy = conj(otfDy);
num1 = abs(otfDx).^2.*fft2(g);

Denom1 = abs(otfDx).^2;
Denom2 = abs(otfDy).^2;

end