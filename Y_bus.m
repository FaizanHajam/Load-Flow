function [nbus, nbranch, Ybus, Yb, theta, G, B]= Y_bus(linedata)
fb=linedata(:,1);
tb=linedata(:,2);
r=linedata(:,3);
x=linedata(:,4);
b2=linedata(:,5);
z=r+x.*1j;
Y=1./z;
nbus=max(max(fb), max(tb));
nbranch=length(fb);
Ybus=zeros(nbus,nbus);
%For off diagonal elements:
for i=1:nbranch
Ybus(fb(i),tb(i))=-Y(i);
Ybus(tb(i),fb(i))=Ybus(fb(i),tb(i));
end
%For diagonal elements:
b2=b2.*1j;
for m=1:nbus
for i=1:nbranch
if fb(i)==m || tb(i)==m
Ybus(m,m)=Ybus(m,m)+Y(i)+b2(i);
end
end
end
Ybus;
Yb=abs(Ybus);
theta=angle(Ybus);
G=real(Ybus);
B=imag(Ybus);