function PHplot(m,L,sitas,PT,HT)
% 压力和膜厚绘图
Ls=(0:m-1)*(L/m); %轴向
[SITA,LL]=meshgrid(sitas,Ls);
subplot(1,2,1);
mesh(SITA,LL,PT')
set(gca,'xtick',0:60:360);xlim([0 360]);xlabel('周向角度(°)'); ylabel('轴向长度(m)'); zlabel('水膜压力(Pa)');
subplot(1,2,2);
mesh(SITA,LL,HT')
set(gca,'xtick',0:60:360);xlim([0 360]);zlim([-inf inf]); xlabel('周向角度(°)'); ylabel('轴向长度(m)'); zlabel('水膜厚度(m)');
saveas(gcf, '膜厚与压力分布.jpg')
end

