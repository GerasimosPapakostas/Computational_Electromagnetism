
syms r theta
%integral
fun=cart2pol(-ux,-uy);
first=int((norm(fun))^2*r,theta,0,2*pi);
we=1/2*8.854*int(first,r,0.76,1.75);
%-----------------------%
%real values of C and We%
%-----------------------%
We=double(we);
C=2*We*10^(-15);
%-----------------------------------------------------%
%discrete form of integral and computation of C and We%
%-----------------------------------------------------%
phi=zeros(3,1);
basis=zeros(3,1);
Wedis=0;
for i=1:Ne
    n(1:3)=t(1:3,i);%nodes of the element
    
    basis(1)=a(1)+b(1)*x(1)+c(1)*y(1);%æi of every element
    basis(2)=a(2)+b(2)*x(2)+c(2)*y(2);
    basis(3)=a(3)+b(3)*x(3)+c(3)*y(3);
     for j=1:3
       phi(j)=basis(j)*X0(t(j,i));
     end
    Wedis=1/2*8.854*transpose(phi)*Se*phi+Wedis;
 
end
 Cdis=2*Wedis*(10^(-12));

%relative error%
relerror=abs(1-Cdis/C)