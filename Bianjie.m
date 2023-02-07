function PK=Bianjie(Dc,ps,H,P,n,m,E,Hz)
PK=zeros(n+1,m+1);
DC=zeros(n+1,m+1);
for i=1:n+1
    for j=1:m+1
        if Hz(i,j)==1
            KK=Dc{i,j};
            DC(i,j)=KK(i,j);
            PK(i,j)=(-H(i,j)*pi*E/(2*ps)+DC(i,j)*P(i,j))/DC(i,j);
        end
    end
end
end

