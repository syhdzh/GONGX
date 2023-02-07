function E=Rubber(YBL,lh)
%计算橡胶的应力-应变关系，利用应力应变模型反推出弹性模量
%参考文献：基于 Mooney-Rivlin 本构模型的捣固凸轮橡胶弹簧变形分析 
%         橡胶压缩Mooney-Ri...lin本构模型参数拟合分析
%测试参数
% YBL=YB/c;
[n,m]=size(YBL);
for i=1:n
    for j=1:m
        if (YBL(i,j)/lh)<0.3 %Mooney-Rivlin模型
            C10=1.9519;C01=-0.3067;
            % r1=(1.+YBL(i,j)).^2;YL(i,j)=2*C10.*r1-2*C01/r1;
            YL(i,j)=C10.*(2+2.*YBL(i,j)-2./(1.+YBL(i,j)).^2)+C01.*(2-2./(1.+YBL(i,j)).^3);
            E(i,j)=10^6*YL(i,j)./YBL(i,j); %应力除以应变直接求解弹性模量
        else %yeoh 模型
            C10=1.8488;C20=-0.937;
            %C30=4.61*10^(-5);C01=0.00524;r1=(1.+YBL(i,j)).^2;YL(i,j)=2.*(C10+2*C20.*(r1-3)+3*C30.*(r1-3).^2).*r1-2*C01./r1;
            r1=1+YBL(i,j);
            YL(i,j)=2*(r1^2-1/r1)*(C10+2*(r1^2+2/r1-3)*C20);
            E(i,j)=10^6*YL(i,j)./YBL(i,j);
        end
    end
end
% plot(E)
% plot(YBL,YL)
end

% %测试参数
% P=ones(129,129);
% xs=1;
% [n,m]=size(P);
% syms YB
% if xs==1 %Mooney-Rivlin模型
%     C10=0.579;C01=4.826*10^(-2);
%     K1=C10*(I1-3)+C01*(I2-3);
%     K2=(J-1)^2/D1;
%     U=K1+K2;
%     YL=C10*(2+2*YB-2/(1+YB)^2)+C01*(2-2/(1+YB)^3);
% elseif xs==2 %yeoh 模型
%     C10=0.696;C20=-0.09147;C30=0.01724;
%     K1=C10*(I1-3)+C20*(I2-3)^2+C30*(I3-3)^3;
%     K2=(J-1)^2/D1+(J-1)^4/D2+(J-1)^6/D3;
%     U=K1+K2;
%     YL=2*(C10+2*C20+3*C30)*(1+YB)^2;
% elseif xs==3 %ogden 3参数模型
%     mu1=1.257;mu2=0.0004081;mu3=5.611*10^(-5);
%     alfa1=1.0985;alfa2=11.983;alfa3=-10.055;
%     Renda1=renda1*J^(-1/3);Renda2=renda2*J^(-1/3);Renda3=renda3*J^(-1/3);
%     K11=2*mu1*(Renda1^alfa1+Renda2^alfa1+Renda3^alfa1-3)/(alfa1^2);
%     K12=2*mu2*(Renda1^alfa2+Renda2^alfa2+Renda3^alfa2-3)/(alfa2^2);
%     K13=2*mu3*(Renda1^alfa3+Renda2^alfa3+Renda3^alfa3-3)/(alfa3^2);
%     K1=K11+K12+K13;
%     K2=(J-1)^2/D1+(J-1)^4/D2+(J-1)^6/D3;
%     U=K1+K2;
% elseif xs==4 %ogden 4参数模型
%     mu1=47.646;mu2=0.0004883;mu3=-82.364;mu4=36.078;
%     alfa1=-1.55;alfa2=9.306;alfa3=-1.825;alfa4=-2.082;
%     Renda1=renda1*J^(-1/3);Renda2=renda2*J^(-1/3);
%     Renda3=renda3*J^(-1/3);Renda4=renda4*J^(-1/3);
%     K11=2*mu1*(Renda1^alfa1+Renda2^alfa1+Renda3^alfa1-3)/(alfa1^2);
%     K12=2*mu2*(Renda1^alfa2+Renda2^alfa2+Renda3^alfa2-3)/(alfa2^2);
%     K13=2*mu3*(Renda1^alfa3+Renda2^alfa3+Renda3^alfa3-3)/(alfa3^2);
%     K14=2*mu4*(Renda1^alfa4+Renda2^alfa4+Renda3^alfa4-3)/(alfa3^2);
%     K1=K11+K12+K13+K14;
%     K2=(J-1)^2/D1+(J-1)^4/D2+(J-1)^6/D3+(J-1)^8/D4;
%     U=K1+K2;
% end
% end

