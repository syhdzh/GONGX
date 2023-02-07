function PK=Pcalcu(H,deltsita,deltL,P,Sx,Sy,Ss,Sc,n,m,VY,VX,csj,fai,R,VV,ps)
dbstop if error
aa=zeros(n,m);bb=zeros(n,m);cc=zeros(n,m);dd=zeros(n,m);kk=zeros(n,m);ee=zeros(n,m);
ff=zeros(n,m);theita=zeros(n,m);
PK=P;
%% ����ϵ�� �μ����׹�ʽ
for i=2:n
    for j=2:m
        theita(i,j)=csj+(i-1)*deltsita-fai;
%         VX(i,j)=e*sin(theita(i,j))/(c*w);
%         VY(i,j)=e*cos(theita(i,j))/(c*w);
        aa(i,j)=Sx*(H(i,j)^2*(H(i,j)-3/4*(H(i+1,j)-H(i-1,j)))/(deltsita^2));
        bb(i,j)=Sx*(H(i,j)^2*(H(i,j)+3/4*(H(i+1,j)-H(i-1,j)))/(deltsita^2));
        cc(i,j)=-2*H(i,j)^3*(1/(deltsita^2)+R^2/(deltL^2));
        dd(i,j)=Sy*R^2*(H(i,j)^2*(H(i,j)-3/4*(H(i,j+1)-H(i,j-1)))/(deltL^2));
        ee(i,j)=Sy*R^2*(H(i,j)^2*(H(i,j)+3/4*(H(i,j+1)-H(i,j-1)))/(deltL^2));
        ff(i,j)=(VV/ps)*(Sc*(H(i+1,j)-H(i-1,j))/(2*deltsita)+2*Ss*(VX*cos(theita(i,j))+VY*sin(theita(i,j))));
        %         ff(i,j)=Sc*(H(i+1,j)-H(i-1,j))/(2*deltsita); %�����������ж� ֻ�����˲̬���ȶ�ѹ���ֲ�
    end
end
%% �������̼���
for i=2:n
    for j=2:m
        % PK(i,j)=(1-namda)*P(i,j)+namda*(ff(i,j)-(aa(i,j)*P(i-1,j)+bb(i,j)*P(i+1,j)+dd(i,j)*P(i,j-1)+ee(i,j)*P(i,j+1)))/cc(i,j); %��������
        kk(i,j)=(aa(i,j)*P(i-1,j)+bb(i,j)*P(i+1,j)+dd(i,j)*P(i,j-1)+ee(i,j)*P(i,j+1));
        PK(i,j)=(ff(i,j)-kk(i,j))/cc(i,j);
        if cc(i,j)==0
            PK(i,j)=0;
        end
        if PK(i,j)<0
            PK(i,j)=-PK(i,j); %��Ĥ���� ����
            break;
        end
        if isnan(PK(i,j))
            PK(i,j)=P(i,j);
        end
    end
end
Pavera=mean(PK(:));
if (isnan(Pavera))
    disp ѹ��ΪNaN������
end
if (isinf(Pavera))
    disp ѹ��Ϊinf������
end
end