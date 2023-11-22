%crank-nicolson%
clc
clear
%N number of point &n number of time step point%
N=203;%we calculate N+2 point for add i-1 ,i+1 point%
xmin=0;
xmax=1;
h=xmax/(N-3);
a=0.2;
tmax=2.5;
T=0.01;
n=(tmax/T)+1;
c=a*T/h;
U=zeros(N,1);
A=zeros(N,N);
B=zeros(N,1);
x=zeros(N,1);
%initial condition%
for k=1:n
   for i=2:N-1
        x(i)=(i-2)*h;
        U(i)= sin(100*pi*x(i));
        U(1)=U(N-2);
        U(N)=U(3);
        B(1)=(c/4)*U(N-2)+U(2)-(c/4)*U(3);
        B(N)=(c/4)*U(N-1)+U(N)-(c/4)*U(3);
        B(i)=(c/4)*U(i-1)+U(i)-(c/4)*U(i+1);
        for i=2:N-1
        A(i,i)=1;
        A(i,i+1)=c/4;
        A(i,i-1)=-c/4;
        A(N,N)=1;
         A(N,N-1)=-c/4;
          A(N,2)=c/4;
           A(1,1)=1;
            A(1,2)=c/4;
             A(1,N-1)=-c/4;
        
        end
      
        
       
       
   end
   U=A^-1*B;
end
 finalenergy=0.5*sqrt(sum(U(2:N-1).^2))
plot(x(2:N-1),(-U(2:N-1)))
