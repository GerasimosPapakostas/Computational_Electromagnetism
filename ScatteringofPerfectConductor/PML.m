E0=1;
R0=10^(-6);%reflection coefficient for vertical incidence
c=3*10^8;%light speed
f0=300*10^6;%frequency
l=c/f0;%wavelength
k0=2*pi*f0/c;
wfreq=2*pi*f0;%angular frequency
a=5*l/2;%scatterer's radiance
w=l;%distance between PML and scatterer
d=l;%PML's thickness
b0=-log(R0)/(2*k0*d);%β
e0=8.854*10^(-12);%dielectric constant
m0=4*pi*10^(-7);%magnetic permeability
a0=1-1j*b0;%a parameter


%---------------------%
%Geometry construction%
%---------------------%

gd=[3     3     3      3      3       3      3       3        3         1 
    4     4     4      4      4       4      4       4       4         0
   -a-w  -a-w   -a-w   a+w    a+w    a+w     a+w    -a-w     -a-w      0
    a+w  -a-w-d  a+w   a+w+d  a+w    a+w    -a-w    -a-w-d   -a-w      a
    a+w  -a-w-d  a+w   a+w+d  a+w+d  a+w+d  -a-w    -a-w-d   -a-w-d    0
   -a-w  -a-w   -a-w   a+w    a+w+d  a+w+d   a+w    -a-w     -a-w-d    0
    a+w   a+w    a+w   a+w    a+w   -a-w     -a-w   -a-w     -a-w      0
    a+w   a+w    a+w   a+w   -a-w   -a-w-d   -a-w   -a-w     a+w       0
   -a-w   a+w+d  a+w+d a+w+d -a-w   -a-w-d   -a-w-d -a-w-d   a+w       0
   -a-w   a+w+d  a+w+d a+w+d  a+w   -a-w     -a-w-d -a-w-d   -a-w      0    ];    
  
  ns=[82 82 82 82 82 82 82 82 82 67;
      49 50 51 52 53 54 55 56 57 49];
  sf='(R1+R2+R3+R4+R5+R6+R7+R8+R9)-C1';
  dl=decsg(gd,sf,ns);
  [p,e,t]=initmesh(dl);[p,e,t]=refinemesh(dl,p,e,t);[p,e,t]=refinemesh(dl,p,e,t);
  [p,e,t]=refinemesh(dl,p,e,t);
  
  
  %-----------------%
  %Computing tensors%
  %-----------------%
  for i=1:9
    if i==9 %air region
        mxx(i,1)=m0;
        myy(i,1)=m0;
        ezz(i,1)=e0;
    elseif i==2 || i==6 % Λx tensor's region
        mxx(i,1)=m0/a0;
        myy(i,1)=m0*a0;
        ezz(i,1)=e0*a0;
        
    elseif i==5 || i==4 %Λy tensor's region
       mxx(i,1)=m0*a0;
       myy(i,1)=m0/a0;
       ezz(i,1)=e0*a0;
    else %ΛxΛy tensor's region
        mxx(i,1)=m0;
        myy(i,1)=m0;
        ezz(i,1)=e0*a0^(2);
    end
  end
Nn=size(p,2);Ne=size(t,2);
node_id=ones(Nn,1);
X0=zeros(Nn,1);
%-----------------------%
%Finding the known nodes%
%-----------------------%
for id=1:Nn
  if p(1,id)^2+p(2,id)^2==a^2
    node_id(id,1)=0;
    X0(id,1)=-E0*exp(-1j*k0*(p(1,i)));
  end
end

%----------------------%
%Computing index vector%
%----------------------%
counter=0;
for i=1:Nn
  if node_id(i,1)==1
    counter=counter+1;
    index(i,1)=counter;
  end
end
%----------------------%
%Computing the matrixes%
%----------------------%
S=spalloc(counter,counter,7*counter);
T=spalloc(counter,counter,7*counter);
A1=spalloc(counter,counter,7*counter);%whole matrix of the system
B=zeros(counter,1);
b=zeros(3,1);
c=zeros(3,1);
for ie=1:Ne
    n(1:3)=t(1:3,ie);%nodes of the element
    rg=t(4,ie);
    x(1:3)=p(1,n(1:3));%x-coordinates of nodes
    y(1:3)=p(2,n(1:3));%y-coordinates of nodes
    De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);
    Ae=abs(De/2);%element area
     %computing ai,bi,ci%
    b(1)=(y(2)-y(3))/De;c(1)=(x(3)-x(2))/De;
    b(2)=(y(3)-y(1))/De;c(2)=(x(1)-x(3))/De;
    b(3)=(y(1)-y(2))/De;c(3)=(x(2)-x(1))/De;
  for i=1:3
      for j=1:3
         Se(i,j)=(b(i)*b(j)/myy(rg,1)+c(i)*c(j)/mxx(rg,1))*Ae;
        if i==j
          Te(i,j)=ezz(rg,1)*Ae/6;
        else
          Te(i,j)=ezz(rg,1)*Ae/12;
        end
          Ad(i,j)=Se(i,j)-wfreq^(2)*Te(i,j);%local matrix of the system
        if node_id(n(i))~=0
            if node_id(n(j))~=0
             S(index(n(i)),index(n(j))) = S(index(n(i)),index(n(j))) + Se(i,j);
             T(index(n(i)),index(n(j)))=T(index(n(i)),index(n(j))) + Te(i,j);
             A1(index(n(i)),index(n(j)))=A1(index(n(i)),index(n(j)))+Ad(i,j);
           
            else
               B(index(n(i)),1) = B(index(n(i)),1) - Ad(i,j)*X0(n(j));
                
            end
        end
      end
  end
  
end
X=A1\B;
%------------------------%
%Passing the known values%
%------------------------%
counter=0;
for i=1:Nn
   
    if index(i)~=0 && X0(i,1)==0
          counter=counter+1;%each time the condition is met,place the element X(counter,1) to the corresponding place
          X0(i,1)=X(counter,1);
        
    end
    
end

%---------------------------%
%Computing the overall field%
%---------------------------%
counter=0;
for i=1:Nn     
    X0(i,1)=E0*exp(-1j*k0*(p(1,i)))+X0(i,1);
    
      if (p(1,i)>a+w || p(1,i)<-a-w) || (p(2,i)<-a-w || p(2,i)>a+w)%excluding the nodes of the PML
          X0(i,1)=NaN;
      end
   
end
%----%
%Plot%
%----%

pdeplot(p,e,t,'XYData',abs(X0));
hold on;
pdegplot(dl);
colormap jet;
axis tight;
axis equal;
figure;
phase=atan2(imag(X0),real(X0));%phase of the wave
pdeplot(p,e,t,'XYData',phase);
hold on;
pdegplot(dl);
colormap jet;
axis tight;
axis equal;

  
  
