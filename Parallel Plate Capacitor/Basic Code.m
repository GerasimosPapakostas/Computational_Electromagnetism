gd=[3  3  3;
    4  4  4;
    -7.5*0.01  -0.03  -0.03;
    -7.5*0.01  -0.03  -0.03;
    7.5*0.01 0.03 0.03;
    7.5*0.01 0.03 0.03;
    7.5*0.01 0.6*0.01 -0.4*0.01 ;
   -7.5*0.01 0.4*0.01  -0.6*0.01;
   -7.5*0.01 0.4*0.01  -0.6*0.01;
    7.5*0.01  0.6*0.01 -0.4*0.01];
ns=char('R1','R2','R3');
ns=ns';
sf='R1-R2-R3';
d1=decsg(gd,sf,ns);

[p,e,t]=initmesh(d1);
[p,e,t]=refinemesh(d1,p,e,t);
%[p,e,t]=refinemesh(d1,p,e,t);
[p,e,t]=refinemesh(d1,p,e,t);


%-------------------------------------------------------------%
%Construction of the node_id vector for known and unkown nodes%
%-------------------------------------------------------------%
Nn=size(p,2);%number of nodes
node_id=ones(Nn,1);
Nd=size(e,2);
Ne=size(t,2);%number of elements
for id=1:Nd
    if (e(6,id)==0 || e(7,id)==0)
        node_id(e(1,id),1)=0;
        node_id(e(2,id),1)=0;        
    end
    if p(1,e(1,id))>=-0.03 && p(1,e(1,id))<=0.03 && p(2,e(1,id))>=-0.01 && p(2,e(1,id))<=0.01
        node_id(e(1,id),1)=0;
    
    end
     if p(1,e(2,id))>=-0.03 && p(1,e(2,id))<=0.03 && p(2,e(2,id))>=-0.01 && p(2,e(2,id))<=0.01
        node_id(e(2,id),1)=0;
    
    end
  
   
   
end

%-----------------------------%
%Construstion of the X0 vector%
%-----------------------------%
X0=zeros(Nn,1);
V=100;
for id=1:Nn
  
  if p(1,id)>=-0.03 && p(1,id)<=0.03 
     
      if((p(2,id)>=0.4*0.01)&&p(2,id)<=0.6*0.01)
        X0(id,1)=V/2;
      elseif((p(2,id)<=-0.4*0.01&& p(2,id)>=-0.6*0.01))
          X0(id,1)=-V/2;
      end 
  end
end

%--------------------------------------------%
%Construction of the re-counting index vector%
%--------------------------------------------%
counter=0;
for i=1:Nn
  if node_id(i,1)==1
      counter=counter+1;
      index(i,1)=counter;
  end
end
%------------------%
%Parse matrix e.t.c%
%------------------%


S=spalloc(counter,counter,7*counter);
B=zeros(counter,1);
Se=zeros(3,3);
b=zeros(3,1);
c=zeros(3,1);
for ie=1:Ne
    n(1:3)=t(1:3,ie);%nodes of the element
    rj=t(4,ie);%region of element
    x(1:3)=p(1,n(1:3));%x-coordinates of nodes
    y(1:3)=p(2,n(1:3));%y-coordinates of nodes
    De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
    Ae=abs(De/2);%element area
    
    %computing bi,ci%
    a(1)=(x(2)*y(3)-x(3)*y(2))/De;
    b(1)=(y(2)-y(3))/De;
    c(1)=(x(3)-x(2))/De;
    a(2)=(x(3)*y(1)-x(1)*y(3))/De;
    b(2)=(y(3)-y(1))/De;
    c(2)=(x(1)-x(3))/De;
    a(3)=(x(1)*y(2)-x(2)*y(1))/De;
    b(3)=(y(1)-y(2))/De;
    c(3)=(x(2)-x(1))/De;
    %computing local rigidity table
    for i=1:3
      for j=1:3
        Se(i,j)=(b(i)*b(j)+c(i)*c(j))*Ae;
          if (node_id(n(i))~=0)
            if(node_id(n(j))~=0)
                S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
            else
                B(index(n(i))) = B(index(n(i))) - Se(i,j)*X0(n(j));
            end
          end
      end  
    end
end
X=S\B;
%------------------------------%
%Filling the X0 matrix%
%------------------------------%
counter=0;
for i=1:Nn
   
    if index(i)~=0 && X0(i)==0
          counter=counter+1;%each time the condition is met,place the element X(counter,1) to the corresponding place
          X0(i,1)=X(counter,1);
        
    end
    
end
%-----%
%Plots%
%-----%
pdeplot(p,e,t,'XYData',X0);%plot of φ0
[ux,uy] = pdegrad(p,t,X0);
figure;
u=-ux;
v=-uy;
pdeplot(p,e,t,'FlowData',[u;v]);
hold on;
pdegplot(d1);




