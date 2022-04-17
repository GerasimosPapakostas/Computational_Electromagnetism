
syms r theta
%integral
fun=([-ux,-uy]);
z=int(norm(fun)^2,r,-7.5*0.01,7.5*0.01);
we=int(z,theta,-7.5*0.01,7.5*0.01);
We1=1/2*8.854*2.2*we*(10^(-12));
C1=2*We1/(100^2);
C=double(C1);
We=double(We1);
phi=zeros(3,1);
basis=zeros(3,1);
Wedis=0;
%ö of every element and discrete form
for i=1:Ne
    n(1:3)=t(1:3,i);%nodes of the element
    
    basis(1)=a(1)+b(1)*x(1)+c(1)*y(1); %æi of every element
    basis(2)=a(2)+b(2)*x(2)+c(2)*y(2);
    basis(3)=a(3)+b(3)*x(3)+c(3)*y(3);
     for j=1:3
       phi(j)=basis(j)*X0(t(j,i));%+basis(2)*X0(t(2,i))+basis(3)*X0(t(3,i));
     end
    Wedis=1/2*8.854*transpose(phi)*Se*phi+Wedis;
    
end
%discrete form


 
 Cdis=2/(100^2)*Wedis*(10^(-12));


