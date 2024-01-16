function [x,w,D,G,Ghat] = GLLmodified(N,xnew)


x = 0;
eps = 1e-7;
Lmat = zeros(N+2,N+1);
D = zeros(N+1,N+1);
G = D;
Ghat = G;


while(max(abs(xnew-x))>eps)
    
    x = xnew;
    
    Lmat(1,:) = 1;
    Lmat(2,:) = x;
    
    for k = 2:N+1
    Lmat(k+1,:) = ((2*k-1).*x.*Lmat(k,:) - (k-1).*Lmat(k-1,:))/k; %from this link: https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
    end
    
    q = Lmat(end,:) - Lmat(end-2,:);
    qb = (2*N+1)*Lmat(end-1,:);
        
    xnew = x - q./qb;
  

end 

 w = 2./(N*(N+1)*Lmat(N+1,:).^2);
 
 
 
 for i = 1:N+1
     for j = 1:N+1
         if (i ~= j)
             D(i,j) = (Lmat(N+1,i)/Lmat(N+1,j))/(x(i)-x(j));
         end
     end
 end
 
 D(1,1) = N*(N+1)/4;
 D(end,end) = - D(1,1);
 
 for j = 2:N
    for n = 2:N
        G(j,n) = sum(D(:,n).*D(:,j).*w(:));
        Ghat(j,n) = G(j,n)/w(j);
    end
end
 
end