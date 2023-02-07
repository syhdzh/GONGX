function [Fp,fp,Fmc,fmc]=fu(HT,PT,c,deltsita,deltL,R,w,Qx,Qy,M,g,nd,DZ,FN,m,n)
%% Çó½âÄ¦²ÁÁ¦
th=zeros(n*DZ,m);
tmc=zeros(n*DZ,m);
for i=1:n*DZ
    for j=1:m        
th(i,j)=nd*w*R/(HT(i,j)*c)+H(i,j)*c*PT(i+1,j)*ps/(deltsita*R);
tmc(i,j)=uc*FN(i,j);
    end
end
Fp=sum(th(:))*deltsita*R*deltL;
Fmc=sum(tmc(:))*deltsita*R*deltL;
W=sqrt((Qx+M*g)^2+Qy^2);
fp=Fp/W;
fmc=Fmc/W;
end