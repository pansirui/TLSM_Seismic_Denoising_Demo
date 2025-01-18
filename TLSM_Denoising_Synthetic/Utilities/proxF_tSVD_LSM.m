function [X, objV] = proxF_tSVD_LSM(Y,rho,opts)

% Input checking
parOP = false;
if ~exist('opts','var') && ~isempty(opts)
    nARG = length(opts);
    if nARG > 0
        parOP=opts{1};
    end    
end

sa = size(Y);

la = length(sa);

[n1,n2,n3] = size(Y);

% t-SVD
[U,S,V] = ntsvd(Y,1,parOP);

sTrueV = zeros(min(n1,n2),n3);

for i = 1:n3
    s = S(:,:,i);  
    s = diag(s);
    sTrueV(:,i) = s;
end

[sTrueV] = proxF_LSM_LR_W(sTrueV,rho);

objV = sum(sTrueV(:));

for i = 1:min(n1,n2) 
    for j = 1:n3
        S(i,i,j) = sTrueV(i,j);
    end
end
    
for i = la:-1:3
    U = ifft(U,[],i);
    S = ifft(S,[],i);
    V = ifft(V,[],i);
end

X = tprod(tprod(U,S), tran(V));

end
