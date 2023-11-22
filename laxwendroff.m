% lax-wendroff method for solving wave equation%
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
u=zeros(n,N+2);
x=zeros(N+2,1);
%initial condition%
for i=2:N+1
    x(i)=(i-2)*h;
    u(1,i)=sin(2*pi*x(i));
end
for k=1:n
    u(k,1)=u(k,N);
    u(k,N+2)=u(k,3);
    for i=2:N+1
     u(k+1,i)=u(k,i)-((c/2)*(u(k,i+1)-u(k,i-1)))+((c^2/2)*(u(k,i-1)-2*u(k,i)+u(k,i+1)));
    end
end
ufinal=u(n,[2:N+1]);
finalenergy=0.5*sqrt(sum(ufinal.^2))
xnew=x(2:N+1);
plot(xnew,ufinal)