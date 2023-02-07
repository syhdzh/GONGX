function Wy=WyCalcu(P,m,n,E0,ps,Dc,WY0c,lh)
P=P*ps;
Wy=zeros(n+1,m+1);
E2=Rubber(WY0c,lh);
E=E0*ones(n+1,m+1);
for i=1:n+1
    for j=1:m+1
        dc=Dc{i,j};
        if isnan(E2(i,j))
        else
            E(i,j)=E2(i,j);
        end
        Wy(i,j)=sum(sum(dc.*P))*2/(pi*E(i,j));
    end
end
end

