function PK2=P2calcu(H,Hz,deltsita,deltL,P,n,m,VV,c,e,w,R,csj,fai,ps)
dbstop if error
aa1=zeros(n,m);bb1=zeros(n,m);cc1=zeros(n,m);dd1=zeros(n,m);ee1=zeros(n,m);ff1=zeros(n,m);kk=zeros(n,m);theita=zeros(n,m);VX=zeros(n,m);VY=zeros(n,m);
PK2=P;
GFG=(P*ps)*0.8*10^(-9);
ND=exp(GFG);
MD=1+(0.6*10^(-9)*P*ps)./(1+1.7*10^(-9)*P*ps);
%% 迭代系数 参见文献公式
for i=2:n
    for j=2:m
        theita(i,j)=csj+(i-1)*deltsita-fai;
        VX(i,j)=e*sin(theita(i,j))/(c*w);
        VY(i,j)=e*cos(theita(i,j))/(c*w);
        %         VX(i,j)=0;
        %         VY(i,j)=0;
        aa1(i,j)=(ND(i,j)/MD(i,j))*((H(i+1,j)+H(i,j))/2)^3;
        bb1(i,j)=(ND(i,j)/MD(i,j))*((H(i-1,j)+H(i,j))/2)^3;
        cc1(i,j)=(ND(i,j)/MD(i,j))*(R^2)*(deltsita^2)*(((H(i,j)+H(i,j+1))/2)^3)/(deltL^2);
        dd1(i,j)=(ND(i,j)/MD(i,j))*(R^2)*(deltsita^2)*(((H(i,j)+H(i,j-1))/2)^3)/(deltL^2);
        ee1(i,j)=(VV/ps)*((H(i+1,j)*MD(i+1,j)-H(i-1,j)*MD(i-1,j))*deltsita/2+2*(VX(i,j)*cos(theita(i,j))+VY(i,j)*sin(theita(i,j)))*(deltsita^2));
        ff1(i,j)=(aa1(i,j)+bb1(i,j)+cc1(i,j)+dd1(i,j));
    end
end
%% 迭代过程计算
for i=2:n
    for j=2:m
        % PK(i,j)=(1-namda)*P(i,j)+namda*(ff(i,j)-(aa(i,j)*P(i-1,j)+bb(i,j)*P(i+1,j)+dd(i,j)*P(i,j-1)+ee(i,j)*P(i,j+1)))/cc(i,j); %加速收敛
        kk(i,j)=(aa1(i,j)*P(i+1,j)+bb1(i,j)*P(i-1,j)+cc1(i,j)*P(i,j+1)+dd1(i,j)*P(i,j-1));
        PK2(i,j)=(kk(i,j)-ee1(i,j))/ff1(i,j);
        if ff1(i,j)==0
            PK2(i,j)=0;
        end
        if PK2(i,j)<0
            PK2(i,j)=0; %油膜破裂 置零
            break;
        end
        if Hz(i,j)==1
            PK2(i,j)=0;
        end
        if isnan(PK2(i,j))
            PK2(i,j)=P(i,j);
        end
    end
end
Pavera=mean(PK2(:));
if (isnan(Pavera))
    disp 压力为NaN，出错
end
if (isinf(Pavera))
    disp 压力为inf，出错
end
end