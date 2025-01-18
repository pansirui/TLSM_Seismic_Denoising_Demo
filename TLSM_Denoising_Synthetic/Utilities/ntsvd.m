%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [U,S,V] = ntsvd(A,fftOP,parOP)
% Nth order tension SVD
% Purpose:  Factors any higher order tensors into its component
%
% Inputs:   A - 3+ order tensor MxNx...
%
%           fftOP (optional) - boolean that determines if outputs have fft
%           Applied along each dimensions
%
%           1 ---- apply fft
%           0 ---- not apply fft
%
%           parOP (optional) - boolean, runs algorithm in parallel.
%
% Output:   U,S,V - (if fftOp true) where U*S*Vt = A.
%
%           or
%
%           U,S,V - (if fftOp false) ifft_T(U)*ifft_T(S)*ifft_T(V)^T = A
%
% Original author :  Misha Kilmer, Ning Hao
% Edited by       :  G. Ely, S. Aeron, Z. Zhang, ECE, Tufts Univ. 03/16/2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U,S,V] = ntsvd(A,fftOP,parOP)

sa = size(A);

la = length(sa);

[n1,n2,n3]=size(A);

fl=0;

if ~exist('parOP','var')
    parOP = false;
end

if ~exist('fftOP','var')
   fftOP = false; 
end

%% Perform FFT along 3 to P axeses of Tensor

if la == 3
    if n2 > n1
        transflag=1;
        A=tran(A);
        nn1=n1;
        n1=n2;
        n2=nn1;
    end
    U = zeros(n1,n1,n3);
    S = zeros(n1,n2,n3);
    V = zeros(n2,n2,n3);
else
    sU =sa;         
    sU(2) = sU(1);
    sV = sa;        
    sV(1) = sV(2);
    U = zeros(sU);  
    S = zeros(sa);
    V = zeros(sV);
end

for i = 3:la
    A = fft(A,[],i);
end

faces = prod(sa(3:la));    

if la == 3
    % Do the conjugate symetric trick here.
    if isinteger(n3/2)
        endValue = int16(n3/2 + 1);
        [U, S, V] = takeSVDs(U,S,V,A,endValue,parOP);
        
        for j =n3:-1:endValue+1
            U(:,:,j) = conj(U(:,:,n3-j+2));
            V(:,:,j) = conj(V(:,:,n3-j+2));
            S(:,:,j) = S(:,:,n3-j+2);
        end

    else 
        endValue = int16(n3/2 + 1);        
        [U,S,V] = takeSVDs(U,S,V,A,endValue,parOP);
        
        for j =n3:-1:endValue+1
            U(:,:,j) = conj(U(:,:,n3-j+2));
            V(:,:,j) = conj(V(:,:,n3-j+2));
            S(:,:,j) = S(:,:,n3-j+2);     
        end
    end 

%% for 4+ dimensional tensors do not perform the
% the conjugate trick.
else 
    [U, S, V] = takeSVDs(U,S,V,A,faces,parOP);
end

%%

if ~fftOP
    [U, S, V] = ifft_T(U,S,V);
end

if exist('transflag','var')
    Uold =U;
    U=V; 
    S=tran(S);
    V=Uold;  
end


end
