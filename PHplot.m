function PHplot(m,L,sitas,PT,HT)
% ѹ����Ĥ���ͼ
Ls=(0:m-1)*(L/m); %����
[SITA,LL]=meshgrid(sitas,Ls);
subplot(1,2,1);
mesh(SITA,LL,PT')
set(gca,'xtick',0:60:360);xlim([0 360]);xlabel('����Ƕ�(��)'); ylabel('���򳤶�(m)'); zlabel('ˮĤѹ��(Pa)');
subplot(1,2,2);
mesh(SITA,LL,HT')
set(gca,'xtick',0:60:360);xlim([0 360]);zlim([-inf inf]); xlabel('����Ƕ�(��)'); ylabel('���򳤶�(m)'); zlabel('ˮĤ���(m)');
saveas(gcf, 'Ĥ����ѹ���ֲ�.jpg')
end

