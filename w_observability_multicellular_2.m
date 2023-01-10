%% Fix k and define system
clear all;
k=3;
rhoo=1;lambda=0.1;
T=0.0003;
A1=expm(T*[0,0,0;0,0,0;0,0,-2e+4]);
A2=expm(T*[0,0,-25e+3;0,0,0;2e+3,0,-2e+4]);
A3=expm(T*[0,0,25e+3;0,0,-25e+3;-2e+3,2e+3,-2e+4]);
A4=expm(T*[0,0,0;0,0,0-25e+3;0,2e+3,-2e+4]);
A5=expm(T*[0,0,0;0,0,25e+3;0,-2e+3,-2e+4]);
A6=expm(T*[0,0,-25e+3;0,0,25e+3;2e+3,-2e+3,-2e+4]);
A7=expm(T*[0,0,25e+3;0,0,0;-2e+3,0,-2e+4]);
A8=expm(T*[0,0,0;0,0,0;0,0,-2e+4]);
A={A1,A2,A3,A4,A5,A6,A7,A8};
C1=[0,0,1];
C={C1,C1,C1,C1,C1,C1,C1,C1};
%% Finding the observable sequences
n=size(A1,2);
if (k<n) error('No observable sequences for the specified k')
end
m=length(A);
I=cell(1,k-n+1); % Matrix cell containing observable sequences of length <=k, inc. constraints
for t=n:1:k % Compute observable sequences
l=[];
x={};
for i=1:t
    x=[x,1:m];

v=allcomb(x{:});%% all possible combinations
end
for i=1:m^t
    if rank(ma(v(i,:),A,C))==n
        I{t-n+1}=[I{t-n+1};v(i,:)];
    end
end


for i=1:size(I{t-n+1},1)
    for j=1:size(I{t-n+1},2)-1
        if(ismember(I{t-n+1}(i,j+1),adm(I{t-n+1}(i,j)))==0)
            l=[l;i];
        end
    end
end
l=unique(l,'rows'); 
if (isempty(I{t-n+1})==0)
I{t-n+1}(l,:)=[]; % remove non admissible observable sequences
end
end

%% remove excess observable sequences
for i=k-n+1:-1:2
    for j=i-1:-1:1
   I{i}=remov(I{j},I{i}); 
    end
end
s=0;
for i=1:length(I)
s=s+isempty(I{i});
end
if(s==length(I)) error('No observable sequences for the specified k')
end

%% create states
for i=1:m
    my_field = strcat('Q0',num2str(i));
    variableQa.(my_field) = sdpvar(n,n);
end
Q={[1:m]'};
for j=2:k-n+2
    Q{j}=unique(I{j-1}(:,1:j),'rows');
end
for j=1:k-n+1
    for i=j+1:length(I)+1
    Q{j}=unique([Q{j};Q{i}(:,1:size(Q{j},2))],'rows');
    end  
end
for i=1:k-1
Qp{i}=str2num(char(Q{i} + '0'));
end
for i=1:k-1
for j=1:size(Qp{i},1)
    my_field = strcat('Q',num2str(Qp{i}(j)));
    variableQ.(my_field) = sdpvar(n,n);
end
end
q=fieldnames(variableQ);
%% fake states that corresponds to the accepted state
qt={};
for i=1:k-n+1
    for j=1:size(I{i},1)
    qt=[qt;strcat('Q',regexprep(num2str(I{i}(j,:)),' ',''))];
    end
end

%% define transitions
for i=1:m
    for b=adm(i)
    my_field = strcat('Y0',num2str(i),num2str(b));
    variableY.(my_field) = sdpvar(n,1);
    end
end
for i=1:length(Qp)
for j=1:size(Qp{i},1)
    y=num2str(Qp{i}(j));
    for b=adm(str2num(y(end)))
    my_field = strcat('Y',num2str(Qp{i}(j)),num2str(b));
    variableY.(my_field) = sdpvar(n,1);
     end
end
end

%% zero transitions LMIs
F=[];
for i=1:m
F=[F,getfield(variableQa,strcat('Q0',num2str(i)))>=eye(n)];
end
for i=1:length(Qp)
    for j=1:size(Qp{i},1)
F=[F, getfield(variableQ,strcat('Q',num2str(Qp{i}(j))))>=eye(n)];
    end
end
%getfield(variableY,strcat('Y0',num2str(i)));
%Y=sdpvar(n,1);
%% initial and intermediate and final transitions LMIs
for i=1:m
    for b=adm(i)
F=[F,[getfield(variableQ,mapa2(strcat('Q0',num2str(i)),b,q,qt)),getfield(variableQ,mapa2(strcat('Q0',num2str(i)),b,q,qt))*A{b}-getfield(variableY,strcat('Y0',num2str(i),num2str(b)))*C{b};(getfield(variableQ,mapa2(strcat('Q0',num2str(i)),b,q,qt))*A{b}-getfield(variableY,strcat('Y0',num2str(i),num2str(b)))*C{b})',rhoo^2*getfield(variableQa,strcat('Q0',num2str(i)))]];
    end 
end
for i=1:length(Qp)
    for j=1:size(Qp{i},1)
        y=num2str(Qp{i}(j));
    for b=adm(str2num(y(end)))
        if (strcmp('Q01',mapa2(strcat('Q',y),b,q,qt))~=1)&&(strcmp('Q02',mapa2(strcat('Q',y),b,q,qt))~=1)&&(strcmp('Q03',mapa2(strcat('Q',y),b,q,qt))~=1)&&(strcmp('Q04',mapa2(strcat('Q',y),b,q,qt))~=1)&&(strcmp('Q05',mapa2(strcat('Q',y),b,q,qt))~=1)&&(strcmp('Q06',mapa2(strcat('Q',y),b,q,qt))~=1)&&(strcmp('Q07',mapa2(strcat('Q',y),b,q,qt))~=1)&&(strcmp('Q08',mapa2(strcat('Q',y),b,q,qt))~=1)
        F=[F,[getfield(variableQ,mapa2(strcat('Q',y),b,q,qt)),getfield(variableQ,mapa2(strcat('Q',y),b,q,qt))*A{b}-getfield(variableY,strcat('Y',y,num2str(b)))*C{b};(getfield(variableQ,mapa2(strcat('Q',y),b,q,qt))*A{b}-getfield(variableY,strcat('Y',y,num2str(b)))*C{b})',rhoo^2*getfield(variableQ,strcat('Q',y))]];
        else
        F=[F,[getfield(variableQa,strcat('Q0',num2str(b))),getfield(variableQa,strcat('Q0',num2str(b)))*A{b}-getfield(variableY,strcat('Y',y,num2str(b)))*C{b};(getfield(variableQa,strcat('Q0',num2str(b)))*A{b}-getfield(variableY,strcat('Y',y,num2str(b)))*C{b})',rhoo^2*lambda^2*getfield(variableQ,strcat('Q',y))]];
        end
    end
end
end
%% Solving the LMIs
su=0;
for i=1:m
    su=su+trace(getfield(variableQa,strcat('Q0',num2str(i))));
end
for i=1:length(Qp)
    for j=1:size(Qp{i},1)
su=su+trace(getfield(variableQ,strcat('Q',num2str(Qp{i}(j)))));
    end
end
optimize(F,su,sdpsettings('solver','sedumi'));
[primalfeas,dualfeas] = check(F);
check(F);

%% Computing the trajectories
for i=1:m
variableQ.(strcat('Q0',num2str(i))) = getfield(variableQ,strcat('Q',num2str(i)));
end
% Number of steps
N=20;
% The trajectory, an array containing the values of the state initialized by x_0
traj1=[-1.5;-0.5;2];
% V, an array containing the Lyapunov function values initialized by V(q_0,x_0)
X={'Q01','Q02','Q03','Q04','Q05','Q06','Q07','Q08'};
qu=string(X(randi(numel(X))));
qu=convertStringsToChars(qu);

V1=[sqrt(traj1'*value(getfield(variableQ,qu))*traj1)];
% The states q_t where the initial state is q_0
% Mode 1 is activated in the first 50 time steps
tmp=adm(str2num(qu(end)));
theta=[tmp(randi(numel(tmp)))];
for i=1:N
    temp =randi([1 3]);
    tt=adm(theta(end));
    theta=[theta tt(temp)];
end
% Update the switching signal value according to the automaton (Figure 2) at each time step
st={qu};yy=[];
% Update the trajectories array and V at each time step
test=[];
for i=1:N
        traj1=[traj1 (A{theta(i)}-(inv(value(getfield(variableQ,mapa2(qu,theta(i),q,qt))))*value(getfield(variableY,strcat('Y',qu(2:end),num2str(theta(i))))))*C{theta(i)})*traj1(:,end)]; 
    V1=[V1 sqrt(traj1(:,end)'*value(getfield(variableQ,mapa2(qu,theta(i),q,qt)))*traj1(:,end))];
    yy=[yy strcat('Y',qu(2:end),num2str(theta(i)))];
    qu=mapa2(qu,theta(i),q,qt);
    %test=[test;[theta(i),q]];
    st=[st qu];
    
end

%% temp
to={};
for i=1:size(st,2)
    tmp=cell2mat(st(i));
    if(strcmp(tmp,'Q01'))
    tmp=strcat('Q',char(949),'1');
    elseif(strcmp(tmp,'Q02'))
    tmp=strcat('Q',char(949),'2');
    elseif(strcmp(tmp,'Q03'))
    tmp=strcat('Q',char(949),'3');
    elseif(strcmp(tmp,'Q04'))
    tmp=strcat('Q',char(949),'4');
    elseif(strcmp(tmp,'Q05'))
    tmp=strcat('Q',char(949),'5');
    elseif(strcmp(tmp,'Q06'))
    tmp=strcat('Q',char(949),'6');
    elseif(strcmp(tmp,'Q07'))
    tmp=strcat('Q',char(949),'7');
    elseif(strcmp(tmp,'Q08'))
    tmp=strcat('Q',char(949),'8');
     end
to=[to tmp(2:end)];
end
st=to;
%% 2nd method
Q=sdpvar(n,n);
F=[Q>=eye(n)];
for i=1:m
    F=[F,A{i}'*Q*A{i}<=rhoo^2*Q];
end
optimize(F,trace(Q),sdpsettings('solver','sedumi'));
[primalfeas,dualfeas] = check(F);
check(F);
Q=value(Q);
%%
mx=[];
Po={};
for i=1:length(I)
    for(j=1:size(I{i},1))
    Po=[Po,mb(I{i}(j,:),A,C,rhoo)];
    end
end


%% dichotomy
y=sdpvar(1,1);
for i=1:length(Po)
F=[F,Po{i}>=y*Q];
end
optimize(F,-y,sdpsettings('solver','sedumi'));
[primalfeas,dualfeas] = check(F);
check(F);
y=value(y);
%% create states
for i=1:m
my_field = strcat('Q0',num2str(i));
  variableQQ.(my_field) = Q;
end
for t=1:k-1
x={};
for i=1:t
    x=[x,1:m];
end
  v=allcomb(x{:});
    for i=1:size(v,1)
    my_field = strcat('Q',regexprep(num2str(v(i,:)),' ',''));
  variableQQ.(my_field) = eye(n);
    end
end
%% fake states that corresponds to the accepted state
qt={};
for i=1:k-n+1
    for j=1:size(I{i},1)
    qt=[qt;strcat('Q',regexprep(num2str(I{i}(j,:)),' ',''))];
    end
end
%% matrices P computation
for i=1:8
variableQQ.(strcat('Q',num2str(i)))=rhoo^2*inv(A{i}')*variableQQ.('Q01')*inv(A{i})+((rhoo^2)/(y*lambda^2))*inv(A{i}')*C{i}'*C{i}*inv(A{i});
end
E=fieldnames(variableQQ);
for i=2*m+1:length(E)
s=E{i};s1=s(1:end-1);s2=s(end);   
variableQQ.(E{i})=rhoo^2*inv(A{str2num(s2)}')*variableQQ.(s1)*inv(A{str2num(s2)})+((rhoo^2)/(y*lambda^2))*inv(A{str2num(s2)}')*C{str2num(s2)}'*C{str2num(s2)}*inv(A{str2num(s2)});
end
%% define gains L
for j=1:m
for i=1:m
my_field = strcat('L0',num2str(i),num2str(j));
  variableLL.(my_field) = ones(n,size(C{1},1));
end
end
for t=2:k
x={};
for i=1:t
    x=[x,1:m];
end
  v=allcomb(x{:});
    for i=1:size(v,1)
    my_field = strcat('L',regexprep(num2str(v(i,:)),' ',''));
  variableLL.(my_field) = ones(n,size(C{1},1));
    end
end
%% gain computing
for j=1:m
for i=1:m
variableLL.(strcat('L0',num2str(j),num2str(i)))=A{i}*inv(variableQQ.('Q01'))*C{i}'*inv(y*lambda^2*eye(size(C{1},1))+C{i}*inv(variableQQ.('Q01'))*C{i}');
end
end
E=fieldnames(variableLL);
for i=m^2+1:length(E)
s=E{i};s1=s(2:end-1);s2=str2num(s(end)) ;   
variableLL.(E{i})=A{s2}*inv(variableQQ.(strcat('Q',s1)))*C{s2}'*inv(y*lambda^2*eye(size(C{1},1))+C{s2}*inv(variableQQ.(strcat('Q',s1)))*C{s2}');
end
%% Simulation
%variableQ.('Q0') = Q;
% Number of steps
% The trajectory, an array containing the values of the state initialized by x_0
traj2=[2;1.5;-0.5];
N=20;

% V, an array containing the Lyapunov function values initialized by V(q_0,x_0)
X={'Q01','Q02','Q03','Q04','Q05','Q06','Q07','Q08'};
qu=string(X(randi(numel(X))));
qu=convertStringsToChars(qu);

V2=[sqrt(traj2'*value(getfield(variableQQ,qu))*traj2)];
% The states q_t where the initial state is q_0
% Mode 1 is activated in the first 50 time steps
%{
tmp=adm(str2num(qu(end)));
theta=[tmp(randi(numel(tmp)))];
for i=1:N
    temp =randi([1 3]);
    tt=adm(theta(end));
    theta=[theta tt(temp)];
end
%}
% Update the switching signal value according to the automaton (Figure 2) at each time step
st={qu};ll=[];
% Update the trajectories array and V at each time step
test=[];
for i=1:N
          traj2=[traj2 (A{theta(i)}-variableLL.(strcat('L',qu(2:end),num2str(theta(i))))*C{theta(i)})*traj2(:,end)]; 
    V2=[V2 sqrt(traj2(:,end)'*variableQQ.(mapa2(qu,theta(i),q,qt))*traj2(:,end))];
    ll=[ll strcat('L',qu(2:end),num2str(theta(i)))];
    qu=mapa2(qu,theta(i),q,qt);
    %test=[test;[theta(i),q]];
    st=[st qu];
    
end
to={};
for i=1:size(st,2)
    tmp=cell2mat(st(i));
    if(strcmp(tmp,'Q01'))
    tmp=strcat('Q',char(949),'1');
    elseif(strcmp(tmp,'Q02'))
    tmp=strcat('Q',char(949),'2');
    elseif(strcmp(tmp,'Q03'))
    tmp=strcat('Q',char(949),'3');
    elseif(strcmp(tmp,'Q04'))
    tmp=strcat('Q',char(949),'4');
    elseif(strcmp(tmp,'Q05'))
    tmp=strcat('Q',char(949),'5');
    elseif(strcmp(tmp,'Q06'))
    tmp=strcat('Q',char(949),'6');
    elseif(strcmp(tmp,'Q07'))
    tmp=strcat('Q',char(949),'7');
    elseif(strcmp(tmp,'Q08'))
    tmp=strcat('Q',char(949),'8');
     end
to=[to tmp(2:end)];
end
st=to;
%% Plotting stuff

figure;
% Plot 2*2 subfigures on the same figure
subplot(5,1,1);% Plot the switching signal \theta_t
stairs([0:N],theta)
axis([0 N 0.9 m+0.1])
xlabel('t')
ylabel('\theta(t)')

subplot(5,1,2);
stairs(0:N,categorical(st));
xlabel('t')
ylabel('q_{t}')

subplot(5,1,3); % Plot the Lyapunov function V(q_t,x_t)
stairs(0:N,log(V1(1:N+1))/log(10))
hold on;
stairs(0:N,log(V2(1:N+1))/log(10))
xlabel('t')
ylabel('log(V(q_{t},x(t)))')
legend("LMI method","explicit construction")

subplot(5,1,4);

stairs([0:N],log(abs((traj1(1,1:N+1))))/log(10),'-')
hold on
stairs([0:N],log(abs(traj1(2,1:N+1)))/log(10),'--')
hold on
stairs([0:N],log(abs(traj1(3,1:N+1)))/log(10),'-.')
legend('log(|x_{1}(t)|)','log(|x_{2}(t)|)','log(|x_{3}(t)|)')
ylabel('log(|x_{i}(t)|)')
xlabel('t')
title('Trajectories using LMI method')

subplot(5,1,5);
stairs([0:N],log(abs((traj2(1,1:N+1))))/log(10),'-')
hold on
stairs([0:N],log(abs(traj2(2,1:N+1)))/log(10),'--')
hold on
stairs([0:N],log(abs(traj2(3,1:N+1)))/log(10),'-.')
legend('log(|x_{1}(t)|)','log(|x_{2}(t)|)','log(|x_{3}(t)|)')
ylabel('log(|x_{i}(t)|)')
xlabel('t')
title('Trajectories using explicit construction')

%------------------------------------------------------------------------------

