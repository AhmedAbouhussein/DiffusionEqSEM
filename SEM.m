clc 
clear all
close all

Nvec = [4,8,16,32,64];
L2mat = zeros(1,length(Nvec));
umat = zeros(max(Nvec),length(Nvec));
xmat = zeros(max(Nvec),length(Nvec));
color='rmbcy';
markers = '+o*.xd';

for k = 1:length(Nvec)
    
N = Nvec(k);
x0 = cos(pi*(N:-1:0)/N); 

u = zeros(1,N+1);
uold = zeros(1,N+1);
f = zeros(1,N+1); %RHS
fold = zeros(1,N+1);
deltaT = 1e-6;
tlimit = 0.1;
t = 0;

exactfunc = @(a,b) sin(pi.*(a+1)).*exp(-pi^2*b);

[x,w,D,G,Ghat] = GLLmodified(N,x0);



u = sin(pi.*(x+1));
exact = exactfunc(x,tlimit);

%[D,G,Ghat] = derivMat(x,w,N);

while (t < tlimit)
    
    uold = u;  
    
    for j = 2:N    
    f(j) = -1*sum(u.*Ghat(j,:));
    end
     
    u = deltaT*((3/2)*f-(1/2)*fold) + uold;    
    
    u(1) = 0;
    u(end) = 0;
    fold = f;
    t = t + deltaT;
    
end

L2 = sum((u-exact).^2);
L2mat(k) = L2;

umat(1:N+1,k) = u; 
xmat(1:N+1,k) = x;

end
figure
hold on
plot(xmat(:,1),umat(:,1),'-o',xmat(:,2),umat(:,2),'-x',xmat(:,3),umat(:,3),'-d',xmat(:,4),umat(:,4),'-*',xmat(:,5),umat(:,5),'-+',xmat(:,5),exact,'r')
legend('N=4','N=8','N=16','N=32','N=64','exact')
xlabel('x')
ylabel('u')
title(['Solution for t = ', num2str(tlimit)])
hold off

%% 

figure
loglog(Nvec,L2mat,'-x')
xlabel('N')
ylabel('L2 Error')
xlim([ 1 1e2])
title('L2 Error vs N')
%%


N = 16;
tvec = [0.05,0.1,0.3];
L2mat = zeros(1,length(tvec));
umat = zeros(N,length(tvec));
xmat = zeros(N,length(tvec));
exactmat = zeros(N,length(tvec));

for k = 1:length(tvec)
    
tlimit = tvec(k);
x0 = cos(pi*(N:-1:0)/N); 

u = zeros(1,N+1);
uold = zeros(1,N+1);
f = zeros(1,N+1); %RHS
fold = zeros(1,N+1);
deltaT = 1e-6;
t = 0;

exactfunc = @(a,b) sin(pi.*(a+1)).*exp(-pi^2*b);

[x,w,D,G,Ghat] = GLLmodified(N,x0);



u = sin(pi.*(x+1));
exact = exactfunc(x,tlimit);

%[D,G,Ghat] = derivMat(x,w,N);

while (t < tlimit)
    
    uold = u;  
    
    for j = 2:N    
    f(j) = -1*sum(u.*Ghat(j,:));
    end
     
    u = deltaT*((3/2)*f-(1/2)*fold) + uold;    
    
    u(1) = 0;
    u(end) = 0;
    fold = f;
    t = t + deltaT;
    
end

L2 = sum((u-exact).^2);
L2mat(k) = L2;

umat(1:N+1,k) = u; 
xmat(1:N+1,k) = x;
exactmat(1:N+1,k) = exact;

end
figure
hold on
plot(xmat(:,1),umat(:,1),'o',xmat(:,1),exactmat(:,1),'-',xmat(:,2),umat(:,2),'x',xmat(:,2),exactmat(:,2),'-',xmat(:,3),umat(:,3),'d',xmat(:,3),exactmat(:,3),'-')
legend('t=0.05 approx','t=0.05 exact','t=0.1 approx','t=0.1 exact','t=0.3 approx','t=0.3 exact')
xlabel('x')
ylabel('u')
title(['Solution for N = ', num2str(N)])
hold off
 
 



