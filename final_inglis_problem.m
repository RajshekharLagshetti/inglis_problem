clc
clear all
%Connectivity Matrix
A=[1,2,3;2,4,5;2,6,3;2,5,6;3,6,7;5,7,6;5,8,7;8,9,10;7,8,10;3,7,11;11,12,13;7,12,11;7,14,12;7,10,14];
%Coordinate Matrix
C=[0,10;7.5,7.5;0,30;10,0;20,0;15,20;30,50;40,0;50,0;50,50;0,70;25,100;0,100;50,100]*0.01;
figure;
col1=(0:1/13:1)';
patch('Faces',A,'Vertices',C,'FaceVertexCData',col1,'FaceColor','flat');
D=((200*10^9)/(1-(0.3)^2))*[1,0.3,0;0.3,1,0;0,0,0.35];
G=zeros(28);
t=0.01;
for i=1:14
    for j=1:3
        x(j)=C(A(i,j),1);
        y(j)=C(A(i,j),2);
    end
    Ae=0.5*det([1,x(1),y(1);1,x(2),y(2);1,x(3),y(3)]);
    B=(0.5/Ae)*[y(2)-y(3),0,y(3)-y(1),0,y(1)-y(2),0;0,x(3)-x(2),0,x(1)-x(3),0,x(2)-x(1);x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)];
    k=0.5*t*Ae*B'*D*B;
    for m=1:3
        for n=1:3
    G(2*A(i,m)-1,2*A(i,n)-1)= G(2*A(i,m)-1,2*A(i,n)-1)+k(2*m-1,2*n-1);
    G(2*A(i,m)-1,2*A(i,n))= G(2*A(i,m)-1,2*A(i,n))+k(2*m-1,2*n);
    G(2*A(i,m),2*A(i,n))= G(2*A(i,m),2*A(i,n))+k(2*m,2*n);
    G(2*A(i,m),2*A(i,n)-1)= G(2*A(i,m),2*A(i,n)-1)+k(2*m,2*n-1);
        end
    end
    x=[];
    y=[];
end
 
f=zeros(28,1);
 
f(24,:)=100*10^6*t*0.25;
f(26,:)=100*10^6*t*0.25*0.5;
f(28,:)=100*10^6*t*0.25*0.5;
 
f(1,:)=[];
G(:,1)=[];
G(1,:)=[];
 
f(5-1,:)=[];
G(:,5-1)=[];
G(5-1,:)=[];
 
f(8-2,:)=[];
G(:,8-2)=[];
G(8-2,:)=[];
 
f(10-3,:)=[];
G(:,10-3)=[];
G(10-3,:)=[];
 
f(16-4,:)=[];
G(:,16-4)=[];
G(16-4,:)=[];
 
f(18-5,:)=[];
G(:,18-5)=[];
G(18-5,:)=[];
 
 
f(21-6,:)=[];
G(:,21-6)=[];
G(21-6,:)=[];
 
f(25-7,:)=[];
G(:,25-7)=[];
G(25-7,:)=[];
 
De=G\f;
GDe=[0;De(1:3);0;De(4:5);0;De(6);0;De(7:11);0;De(12);0;De(13:14);0;De(15:17);0;De(18:20)];
si=[];
for i=1:14
    for j=1:3
        x(j)=C(A(i,j),1);
        y(j)=C(A(i,j),2);
    end
    Ae=0.5*det([1,x(1),y(1);1,x(2),y(2);1,x(3),y(3)]);
    B=(0.5/Ae)*[y(2)-y(3),0,y(3)-y(1),0,y(1)-y(2),0;
        0,x(3)-x(2),0,x(1)-x(3),0,x(2)-x(1);
        x(3)-x(2),y(2)-y(3),x(1)-x(3),y(3)-y(1),x(2)-x(1),y(1)-y(2)];
    de=[];
    for k=1:3
    de=[de;GDe(((2*A(i,k))-1),1);GDe((2*A(i,k)),1)];
    end
    e=B*de;
    si=[si,D*e];
end
figure;
col2=(si(1,:))';
patch('Faces',A,'Vertices',C,'FaceVertexCData',col2,'FaceColor','flat');
colorbar;
figure;
col3=(si(2,:))';
patch('Faces',A,'Vertices',C,'FaceVertexCData',col3,'FaceColor','flat');
colorbar;
figure;
col4=(si(3,:))';
patch('Faces',A,'Vertices',C,'FaceVertexCData',col3,'FaceColor','flat');
colorbar;
 
