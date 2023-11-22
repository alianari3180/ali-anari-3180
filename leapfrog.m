% leap frog for solving wave equation%
clc
clear
%N number of point &n number of time step point%
N=201;
xmin=0;
xmax=1;
h=(xmax)/(2*(N-1));
a=0.2;
tmax=2.5;
T=0.01;
n=(tmax/T)+1;
c=a*T/h;
u=zeros(2*n,2*N+2);
x=zeros(2*N+2,1);
%initial condition%
for i=2:2*N+1
    x(i)=(i-2)*h;
    u(1,i)=sin(100*pi*x(i));
end
for k=2:2*n
    u(k,1)=u(k,2*N);
    u(k,2*N+2)=u(k,3);
    for i=2:2*N+1
        u(k+1,i)=u(k-1,i)-c*(u(k,i+1)-u(k,i-1));
      end
end
ufinal=u(n,[2:2:2*N+1]);
finalenergy=0.5*sqrt(sum(ufinal.^2))
xnew=x(2:2:2*N+1);
plot(xnew,ufinal)