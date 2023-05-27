%% tapered helicopter blade problem
clear all 
clc
L=10;
N=100;
le=L/N;
Po=20000;
E=220*10^3;
%% area  vector
 for i=1:N+1;
        x(i)=(i-1)*le;
 end
    
 A(1)=100;

for i=1:N
   A(i+1)=A(1)*(exp(-(x(i+1))/(2*L)));
   
   
end  
%% global stiffness matrix    
GK=zeros(N+1,N+1);
for e=1:N
    alpha=(E*(A(e)+A(e+1)))/(2*le);
    EK=[1 -1;-1 1]*alpha;
    for i=1:2
      for  j=1:2
        GK(e+i-1,e+j-1)=GK(e+i-1,e+j-1)+EK(i,j);
      end
    end
end
GN=zeros(N+1,1);    
for e=1:N
     beta=Po/(3*A(1));
     EL=[((2*A(e))+(A(e+1)));(A(e)+2*A(e+1))]*beta;
     
     for i=1:2
             GN(e+i-1,1)=GN(e+i-1,1)+EL(i,1);
     end
  end
n=input('\n Enter the node at which deformation is zero');
GK(n,:)=[];
GK(:,n)=[];
GN(n)=[];
us=inv(GK)*GN;
U=[0;us];
plot(x,U,'b');
hold on;
y=[0:0.1:10];
uexact=(-((Po*le*le/2)*y)+((Po*(y.^3))/6))/(E*A(1));
plot(y,uexact,'r');
error=U-uexact;

z=input('\n Enter the value of x where u want to know deformation');

n1=fix(z/(L/N))+1;
n2=n1+1;
N1=(x(n2)-z)/(x(n2)-x(n1));
N2=(z-x(n1))/(x(n2)-x(n1));
Uz=(N1*U(n1))+(N2*U(n2))
