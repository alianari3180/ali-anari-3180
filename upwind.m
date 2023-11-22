% upwind method for solving wave equation%
clc
clear
%N number of point &n number of time step point%
N=201;
xmin=0;
xmax=1;
h=xmax/(N-1);
a=0.2;
tmax=2.5;
T=0.02;
n=(tmax/T)+1;
c=a*T/h;
u=zeros(n,N);
x=zeros(N,1);
%initial condition%
for i=1:N
    x(i)=(i-1)*h;
    u(1,i)=sin(100*pi*x(i));
end
for k=1:n
    u(k,1)=u(k,N);
    for i=2:N
        u(k+1,i)=u(k,i)-(c*(u(k,i)-u(k,i-1)));
    end
end
ufinal=u(n,:);
finalenergy=0.5*sqrt(sum(ufinal.^2))
plot(x,ufinal)

    
    
    
    

