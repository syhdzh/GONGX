function Dh=DhCalcu(R,csj,DZ,JHcs,ks)
%计算不同形状的沟槽的膜厚
%参考文献 《径向推力一体式水润滑轴承多参数流固耦合润滑数值模型》
%%____________testing data________________
% R=0.15;DZ=6;csj=0:pi/720:2*pi-pi/720;JHcs=[0.012];ks=3;
%%____________formal fucntion________________
da=2*pi/DZ;m=length(csj);Dh=zeros(m,1);
%%长方形结构
if ks==1
    Ry=JHcs(1);depth=JHcs(2);
    Ja=Ry/R;
    for i=1:m
        fs=mod(csj(i),da);
        if (0<=fs)&&(fs<=Ja)
            Dh(i)=depth;
        elseif (Ja<fs)&&(fs<da-Ja)
            Dh(i)=0;
        elseif (da-Ja<=fs)&&(fs<=da)
            Dh(i)=depth;
        end
    end
    
%%三角形结构
elseif ks==2
    Ry=JHcs(1);depth=JHcs(2);
    Ja=Ry/R;
    for i=1:m
        fs=mod(csj(i),da);
        if (0<=fs)&&(fs<Ja)
            St1=R*cos(Ja)/cos(fs);
            Dh(i)=St1+depth*(1-tan(fs)/tan(Ja))-R;
        elseif (Ja<=fs)&&(fs<(da-Ja))
            Dh(i)=0;
        elseif ((da-Ja)<=fs)&&(fs<da)
            St1=R*cos(Ja)/cos(da-fs);
            Dh(i)=St1+depth*(1-tan(da-fs)/tan(Ja))-R;
        end
    end
    
%%半椭圆形结构
elseif ks==3
    Ry=JHcs;
    Ja=Ry/R;
    for i=1:m
        fs=mod(csj(i),da);
        if fs==0
            Dh(i)=Ry;
        elseif (0<fs)&&(fs<=Ja)
            St1=asin(R*sin(fs)/Ry);
            Dh(i)=Ry*sin(pi-fs-St1)/sin(fs)-R;
        elseif (Ja<fs)&&(fs<=da-Ja)
            Dh(i)=0;
        elseif (da-Ja<fs)&&(fs<=da)
            St1=asin(R*sin(da-fs)/Ry);
            Dh(i)=Ry*sin(pi-(da-fs)-St1)/sin(da-fs)-R;
        end
    end
    
%%圆形带倒角
elseif ks==4
    R3=JHcs(1);R4=JHcs(2);
    Ja1=0;
    Jc1=acos((R^2+(R+R4)^2-(R3+R4)^2)/(2*R*(R+R4)));
    Jo234=asin((R4+R)*sin(Jc1)/(R3+R4));
    O2B=sqrt(R^2+R3^2-2*R*R3*cos(Jo234));
    Jb1=asin(R3*sin(Jo234)/O2B);
    Jd1=da-Jc1;Je1=da-Jb1;Jf1=da;
    Dh=zeros(m,1);
    for i=1:m
        fs=mod(csj(i),da);
        if fs==0
            Dh(i)=R3;
        elseif (Ja1<fs)&&(fs<=Jb1)
            St1=asin(R*sin(fs-Ja1)/R3);
            Dh(i)=R3*sin(pi-(fs-Ja1)-St1)/sin(fs-Ja1)-R;
        elseif (Jb1<fs)&&(fs<=Jc1)
            St2=pi-asin((R+R4)*sin(Jc1-fs)/R4);
            Dh(i)=R4*sin(pi-(Jc1-fs)-St2)/sin(Jc1-fs)-R;
        elseif (Jc1<fs)&&(fs<=Jd1)
            Dh(i)=0;
        elseif (Jd1<fs)&&(fs<=Je1)
            St4=pi-asin((R+R4)*sin(fs-Jd1)/R4);
            Dh(i)=R4*sin(pi-(fs-Jd1)-St4)/sin(fs-Jd1)-R;
        elseif (Je1<fs)&&(fs<=Jf1)
            St5=asin(R*sin(Jf1-fs)/R3);
            Dh(i)=R3*sin(pi-(Jf1-fs)-St5)/sin(Jf1-fs)-R;
        end
    end
end
end


