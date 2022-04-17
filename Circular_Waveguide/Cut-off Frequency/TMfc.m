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

Nn=size(p,2);Ne=size(t,2);
node_id=ones(Nn,1);
X0=zeros(Nn,1);
%-----------------------%
%Finding the known nodes%
%-----------------------%
for id=1:Nn
  if p(1,id)^2+p(2,id)^2==C(4,1)^2
    node_id(id,1)=0;
  end
end
%---------------------------%
%Computing recounting vector%
%---------------------------%
counter=0;
for i=1:Nn
  if node_id(i,1)==1
    counter=counter+1;
    index(i,1)=counter;
  end
end
%------------------------------------------%
%computing rigidity and local mass matrixes%
%------------------------------------------%


S=spalloc(counter,counter,7*counter);
T=spalloc(counter,counter,7*counter);
b=zeros(3,1);
c=zeros(3,1);
for ie=1:Ne
    n(1:3)=t(1:3,ie);%nodes of the element
    x(1:3)=p(1,n(1:3));%x-coordinates of nodes
    y(1:3)=p(2,n(1:3));%y-coordinates of nodes
    De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
    Ae=abs(De/2);%element area
    %computing ai,bi,ci%
    b(1)=(y(2)-y(3))/De;c(1)=(x(3)-x(2))/De;
    b(2)=(y(3)-y(1))/De;c(2)=(x(1)-x(3))/De;
    b(3)=(y(1)-y(2))/De;c(3)=(x(2)-x(1))/De;
    %computing local rigidity and local mass matrixes
    for i=1:3
      for j=1:3
        Se(i,j)=(b(i)*b(j)+c(i)*c(j))*Ae;
        if i==j
          Te(i,j)=Ae/6;
        else
            Te(i,j)=Ae/12;
        end
          if node_id(n(i))~=0
              if node_id(n(j))~=0
                S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
                T(index(n(i)),index(n(j)))=T(index(n(i)),index(n(j))) + Te(i,j);
              end
          end
       
           
       end
    end
end  
%----------------------------------------%
%finding the eigenvectors and eigenvalues%
%----------------------------------------%
[V,D] = eigs(S,T,13,'smallestabs');


%---------------------------------------%
%Filling X0 vector with the eigenvectors%
%---------------------------------------&


for columns=1:13
    counter2=0;
  for i=1:Nn
    if p(1,i)^2+p(2,i)^2~=C(4,1)^2 %if node is not on the edge then pass the value
           counter2=counter2+1;

        X0(i,columns)=V(counter2,columns);
    end
  end
end



%---------------------------%
%Plotting and relative error%
%---------------------------%

fc_theoritical=[11.48302914;18.29645226;24.51781898;26.35605858;30.46225611;33.49893242]*10^9;
fc=zeros(6,1);
relativerror=zeros(6,1);
counter=1;
for i=1:13
   if i==1 || i==2 ||i==5 ||i==6 || i==8 ||i==9 %desired modes
      figure;
      pdeplot(p,e,t,'XYData',X0(:,i),'Colormap','jet'); 
      axis tight; 
      axis equal; 

      kc=sqrt(D(i,i));
      fc(counter,1)=3*10^(8)*kc/(2*pi);
      relativerror(counter,1)=abs((fc(counter,1)-fc_theoritical(counter,1)))/fc_theoritical(counter,1)*100;
      counter=counter+1;
  end
end

