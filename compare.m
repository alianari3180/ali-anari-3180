%compare plots%
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
T=0.002;
n=(tmax/T)+1;
c=a*T/h;
u=zeros(n,N);
x=zeros(N,1);
%initial condition%
for i=1:N
    x(i)=(i-1)*h;
    if 0.1<=x(i)<=0.3
        u(1,i)=1;
    else 
        u(1,i)=0;
    end
    
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

hold on

    
    
    
    

% lax method for solving wave equation%
clc
clear
%N number of point &n number of time step point%
N=201;
xmin=0;
xmax=1;
h=xmax/(N-1);
a=0.2;
tmax=2.5;
T=0.002;
n=(tmax/T)+1;
c=a*T/h;
u=zeros(n,N+2);
x=zeros(N+2,1);
%initial condition%
for i=2:N+1
    x(i)=(i-2)*h;
     if 0.1<=x(i)<=0.3
        u(1,i)=1;
    else 
        u(1,i)=0;
    end
    
end
for k=1:n
    u(k,1)=u(k,N);
    u(k,N+2)=u(k,3);
    for i=2:N+1
     u(k+1,i)=((u(k,i-1)+u(k,i+1))/2)-((c/2)*(u(k,i+1)-u(k,i-1)));
    end
end
ufinal=u(n,[2:N+1]);
finalenergy=0.5*sqrt(sum(ufinal.^2))
xnew=x(2:N+1);
plot(xnew,ufinal);
hold on
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
T=0.002;
n=(tmax/T)+1;
c=a*T/h;
u=zeros(n,N+2);
x=zeros(N+2,1);
%initial condition%
for i=2:N+1
    x(i)=(i-2)*h;
     if 0.1<=x(i)<=0.3
        u(1,i)=1;
    else 
        u(1,i)=0;
    end
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
hold on
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
T=0.002;
n=(tmax/T)+1;
c=a*T/h;
u=zeros(2*n,2*N+2);
x=zeros(2*N+2,1);
%initial condition%
for i=2:2*N+1
    x(i)=(i-2)*h;
     if 0.1<=x(i)<=0.3
        u(1,i)=1;
    else 
        u(1,i)=0;
    end
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
hold on
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
T=0.002;
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
         if 0.1<=x(i)<=0.3
        u(1,i)=1;
    else 
        u(1,i)=0;
    end
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
plot(x(2:N-1),-U(2:N-1))
legend('upwind','lax','lax-wendroff','leap frog','crank-nicolson')
       
       
       


    
    
    
    



