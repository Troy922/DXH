function f=eqs(t,y,f1,f2,f3,f4,f5,f6,f7,f8)
%f=zeros(2,1001);
f(1,1)=(f1-f2*y(1))./f7;
f(2,1)=(f2*y(1)+f3-f4*y(2)-f5*(y(2)-f6))./f8;