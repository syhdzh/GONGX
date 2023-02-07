function [Fx,Fy,F,Xa,Ya,Xv,Yv,X,Y,e,Dt,fai]=Fcalcu(deltsita,ps,PT,FN,DZ,R,L,c,Qx,Qy,M,g,Dt,X0,Y0,Xv0,Yv0,dyt,dxt,n,m)
dbstop if error
%   更新迭代
PT=PT*ps; %恢复成有量纲量
sitas=(0:n*DZ-1)*deltsita; %周向
for i=1:m
    A(:,i)=sin(sitas');
    B(:,i)=cos(sitas');
end
PA=PT.*A;PB=PT.*B;
Fx=-(sum(sum(PB))*deltsita*R*L/m);
Fy=(sum(sum(PA))*deltsita*R*L/m);
F=sqrt(Fx^2+Fy^2);
Xa=(-Fx+Qx+M*g)/M;
Ya=(-Fy+Qy)/M;
% if Xa<0
%     Xa=Qx/M;
%     Ya=Qy/M;
% end
Xv=Xv0+Xa*Dt;Yv=Yv0+Ya*Dt;
X=X0+Xv*Dt;Y=Y0+Yv*Dt;
e=sqrt(X^2+Y^2)/c;
kk=1;
while e>1
    disp 请减小Dt
    Dt=Dt/2;
    Xv=Xv0+Xa*Dt*dxt;Yv=Yv0+Ya*Dt*dyt;
    X=X0+Xv*Dt;Y=Y0+Yv*Dt;
    e=sqrt(X^2+Y^2)/c;
%     if kk>5
%         Y=Y*0.8;
%         X=X*0.8;
%     end
    kk=kk+1;
end
if Fy==0
    fai=0;
else
    fai=atan(Fy/Fx);
end
fai=abs(fai);
end