function m=mapa2(q,i,list,qt)
ind=2;
    t=strcat(q,num2str(i));
     if strcmp(q,'Q01')||strcmp(q,'Q02')||strcmp(q,'Q03')||strcmp(q,'Q04')||strcmp(q,'Q05')||strcmp(q,'Q06')||strcmp(q,'Q07')||strcmp(q,'Q08')
        m=strcat('Q',num2str(i));
     
     elseif sum(ismember(qt, t))==1
         m=strcat('Q0',num2str(i));
     else
         
    while(ind~=1)
    if sum(ismember(list, t))==1
     m=t;ind=1;
     else t=eraseBetween(t,2,2);
    end
    end
     end
   
end
