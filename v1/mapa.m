% This function acts like a transition function, it takes 4 arguments, a
% state q, an index i, a list of states and a list of fake states defined
% before.
function m=mapa(q,i,list,qt)
ind=2; % indicator function to enter the while loop
    t=strcat(q,num2str(i));
     if strcmp(q,'Q0')
        m=strcat('Q',num2str(i));
     elseif sum(ismember(qt, t))==1
% if (i_1,...,i_n) is an observable sequence, then create the
% transition: \delta(Q_{1...n-1},i_n)=Q_0 using the fake states.
         m='Q0';
     else
         
    while(ind~=1)
% apply step 24 in the algorithm
    if sum(ismember(list, t))==1
     m=t;ind=1;
     else t=eraseBetween(t,2,2);
    end
    end
     end
   
end
