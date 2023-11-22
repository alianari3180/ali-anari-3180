%initial energy%
N=201;
h=1/(N-1);
x=zeros(N,1);
u=zeros(N,1);
for i=1:N
    x(i,1)=(i-1)*h;
    u(i,1)=sin(2*pi*x(i));
end
initialenerg=0.5*sqrt(sum(u.^2))
    