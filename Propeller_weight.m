function [Fm,Tm,Fn]=Propeller_weight(md,N,bmax,t02,t06,d,D,PD,W,rm,px,py)
%___________input variable__________________
% md 桨叶密度
% N 桨叶数量
% bmax 螺旋桨最大宽度
% t02 0.2R处螺旋桨的最大厚度
% t06 0.6R处螺旋桨的最大厚度
% d 桨毂直径
% D 螺旋桨直径
% PD 主机最大持续功率情况下螺旋桨收到的马力,单位kw
% W 螺旋桨在PD下的转速,单位r/mi
% rm 螺旋桨叶片的质心半径
% px 主轴偏心x距离
% py 主轴偏心y距离

% Fn 轴承承受的径向载荷
% Db 水润滑轴承内径
% Lb 水润滑轴承的轴向长度（宽）
% au bu cu du 摩擦系数公式计算参数
% vr 转子和水润滑轴承之间的相对滑动速度
%-----------output variable-----------------
% Fm 由于螺旋桨重力偏心引起的径向载荷
% Tm 由于螺旋桨重力偏心引起的扭矩
% Fnm 不平衡载荷对于轴承正压力的影响
%______________testing data____________________
md=7600;N=3;bmax=0.132955656;
t02=0.155114932;t06=0.070909683;
d=0.886371041;D=4.3432181;
PD=850; W=100; rm=t02*0.25+t06*0.75; %可用solidworks测算出螺旋桨质心的位置
px=0.0001;py=0.0002;
au=0.0029;bu=4.621;cu=1.015;du=-0.3596;%见论文水润滑橡胶轴承摩擦特性及其诱导的螺旋桨轴系振动研究
%________________Weight caculation________________
ta=0.5*t02+t06;
da=(1-d/D)*D;
Gb=0.169*md*N*bmax*da*ta; %桨叶的重量,单位kg
Lk=d+0.1; %毂长，单位m
K=1/13; %轮毂配合的锥度
d0=0.045+0.108*(PD/W)^(1/3)-K*Lk/2; %桨毂长度中央处轴径,m
Gm=(0.88-0.6*d0/d)*Lk*md*d^2; %桨毂重量,kg
Imp=9.48*md*N*bmax*ta*D^3; %必须在d/D小于0.18的情况，单位kg`m`s2
%_______________weight center caculation____________________
len=360; %划分点
a0=0:360/N:360*(N-1)/N; %各个叶片初始角度
x=zeros(len,1);y=zeros(len,1);rd=zeros(len,1); a=zeros(len,N);ar=zeros(len,1);
xa=zeros(len,N);ya=zeros(len,N);
for i=1:len
    for j=1:N
        a(i,j)=mod((a0(j)+(i-1)*360/len),360); %叶片质心角度
        xa(i,j)=rm*cosd(a(i,j));
        ya(i,j)=rm*sind(a(i,j));
    end
    x(i)=sum(xa(i,:))/N; %螺旋桨质心x坐标
    y(i)=sum(ya(i,:))/N; %螺旋桨质心y坐标
    if abs(x(i))<1e-16
        x(i)=0;
    end
    if abs(y(i))<1e-16
        y(i)=0;
    end
    rd(i)=sqrt((x(i)-px)^2+(y(i)-py)^2); %螺旋桨质心偏移旋转中心距离
    if y(i)>=0
        ar(i)=atan(y(i)/x(i));
        if isnan(ar(i))
            ar(i)=pi/2;
        end
    else
        ar(i)=atan(y(i)/x(i))+pi;
        if isnan(ar(i))
            ar(i)=pi*3/2;
        end
    end
end
% scatter (x,y)
%______________________Fm calculation________________________
Fm=Gb+Gm;
Tm=Gb*rd+Gm*sqrt(px^2+py^2);
Fc=Gb*rd*(W*pi/30)^2;
alfaf=ones(len,1);
% alfaf，转子配重安装处至水润滑橡胶艉轴承支点的动态载荷传递效率
jd=(0:360/len:360*(len-1)/len)';
Fnm=alfaf.*Fc.*sind(jd);%不平衡载荷对于轴承正压力的影响
plot(jd,Fn)
%不平衡载荷引起的摩擦系数变化
P=Fn/(Db*Lb); %原载荷引起的面压力
PM=Fnm/(Db*Lb); %不平衡载荷引起的面压力
uf=au*exp(bu*(cu+abs(vr))^(-1))*(P+PM)^du;%更新后的摩擦系数

end
