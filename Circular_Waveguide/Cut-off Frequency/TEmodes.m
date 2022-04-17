%--------%
%GEOMETRY%
%--------%

C=[1
    0
    0
    1*10^(-2)];
gd=[1
    0
    0
    1*10^(-2)];

d1=decsg(gd);
[p,e,t]=initmesh(d1);[p,e,t]=refinemesh(d1,p,e,t);[p,e,t]=refinemesh(d1,p,e,t);
[p,e,t]=refinemesh(d1,p,e,t);


%----------------------------------%
%Computing rigidity and mass tables%
%----------------------------------%
Nn=size(p,2);%number of nodes
Ne=size(t,2);%number of elements
S=spalloc(Nn,Nn,7*Nn);T=spalloc(Nn,Nn,7*Nn);
b=zeros(3,1);c=zeros(3,1);
for ie=1:Ne
    n(1:3)=t(1:3,ie);%nodes of the element
    x(1:3)=p(1,n(1:3));%x-coordinates of nodes
    y(1:3)=p(2,n(1:3));%y-coordinates of nodes
    De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
    Ae=abs(De/2);%element area
    %computing bi,ci%
    b(1)=(y(2)-y(3))/De;c(1)=(x(3)-x(2))/De;
    b(2)=(y(3)-y(1))/De;c(2)=(x(1)-x(3))/De;
    b(3)=(y(1)-y(2))/De;c(3)=(x(2)-x(1))/De;
    
    for i=1:3
      for j=1:3
        Se(i,j)=(b(i)*b(j)+c(i)*c(j))*Ae;
        S((n(i)),(n(j))) = S(n(i),n(j)) + Se(i,j);
        if i==j
          Te(i,j)=Ae/6;
        else
            Te(i,j)=Ae/12;
        end
        T(n(i),n(j))=T(n(i),n(j)) + Te(i,j);
           
       end
    end
end  
[V,D] = eigs(S,T,13,'smallestabs');


fc_theoritical=[8.790127507;14.58177589;18.29645226;20.05829748;25.38680497;25.45365005]*10^9;
fc=zeros(6,1); 
relativerror=zeros(6,1);

counter=1;
for i=2:1:13
    
    if i==3 || i==4 || i==6 || i==8 || i==9  || i==12 %desired modes met for these i's
      X0(:,counter)=V(:,i);
      figure;
      pdeplot(p,e,t,'XYData',X0(:,counter),'Colormap','jet');

      axis tight; 
      axis equal; 
      kc=sqrt(D(i,i));
      fc(counter,1)=3*10^(8)*kc/(2*pi);
      relativerror(counter,1)=abs((fc(counter,1)-fc_theoritical(counter,1)))/fc_theoritical(counter,1)*100;
      counter=counter+1;
    end
end

