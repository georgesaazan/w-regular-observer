%% This produces the results of the Numerical Example in the paper
% This ouputs a figure corresponding to the fast switching signal, in order
% to get the second figure corresponding to the slow swtiching signal,
% check line 125.
%% define system and fix parameters
clear all;clc;
rhoo=1;lambda=0.1; %fix \rho and \lambda
k=3; %fik k the maximal length of observable sequences
T=0.0003;%define the sampling time and the system
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
n=size(A{1},2); %dimension of the system
m=length(A); %number of modes
p=size(C{1},1); %dimension of outputs
%% 1st method LMI solving
%% Compute the obserable and the minimal observable sequences
if (k<n) error('No observable sequences for the specified k')
end
O={}; %initialization of a cell that contains the observable sequences of length n to k separately
for t=n:1:k
x={};
for i=1:t
    x=[x,1:m];%intermediate variable for the next step
end
v=allcomb(x{:});%create all possible combinations of switching sequences of length up to t, n\le t \le k
for i=1:length(v)
    if rank(ma(v(i,:),A,C,1))==n %check if a combination is an observable sequence
        O=[O;v(i,:)];%add the observable sequence to I
    end
end
end
for i=1:length(O)
    if size(O{i},2)~=1
        for j=1:size(O{i},2)-1
            if ~ismember(O{i}(:,j+1),adm((O{i}(:,j))))
                O{i}=[]; break;
    end
        end
    end
end
O=O(~cellfun('isempty',O));
for i=1:length(O)
    S={};P={};
    if(size(O{i},2)~=1)
    for j=1:size(O{i},2)-1
    P=[P;O{i}(:,1:j)];
    end
    for j=size(O{i},2):-1:2
    S=[S;O{i}(:,j:end)];
    end
    for z=1:length(P)
        for zz=1:length(S)
    if(sum(find(cellfun(@(O) isequal(O, P{z}), O)))~=0||sum(find(cellfun(@(O) isequal(O, S{zz}), O)))~=0)
        O{i}=[]; 
    end
        end
    end
    end
end
O=O(~cellfun('isempty',O));
%% Steps 1 till 9 of the algorithm
for i=1:m
variableQ.(strcat('Q0',num2str(i))) = sdpvar(n,n); %initial state
end
delta=containers.Map;   
for i=1:m
    for j=adm(i)
    variableY.(strcat('Y0',num2str(i),num2str(j))) = sdpvar(n,p);
if(find(cellfun(@(O) isequal(O, j), O)))
     delta=[delta;containers.Map([strcat('Q0',num2str(i)),j],'Q0','UniformValues',false)];
else
     variableQ.(strcat('Q',num2str(j))) = sdpvar(n,n);
     delta=[delta;containers.Map([strcat('Q0',num2str(i)),j],strcat('Q',num2str(j)),'UniformValues',false)];
end
    end
end
%% Steps 10 till 20 of the algorithm
for i=1:length(O)
    l=size(O{i},2);
   if l>=2
       if l==2
         delta=[delta;containers.Map([strcat('Q',O{i}(:,1)),O{i}(:,2)],strcat('Q0',num2str(O{i}(:,2))),'UniformValues',false)];
        variableY.(strcat('Y',regexprep(num2str(O{i}(:,1:l)),' ',''))) = sdpvar(n,p);
       else
           for j=2:l-1
           variableQ.(strcat('Q',regexprep(num2str(O{i}(:,1:j)),' ',''))) = sdpvar(n,n);
           delta=[delta;containers.Map([strcat('Q',regexprep(num2str(O{i}(:,1:j-1)),' ','')),O{i}(:,j)],strcat('Q',regexprep(num2str(O{i}(:,1:j)),' ','')),'UniformValues',false)];
           variableY.(strcat('Y',strcat(regexprep(num2str(O{i}(:,1:j-1)),' ',''),num2str(O{i}(:,j))))) = sdpvar(n,p);
           end
           delta=[delta;containers.Map([strcat('Q',regexprep(num2str(O{i}(:,1:l-1)),' ','')),O{i}(:,l)],strcat('Q0',num2str(O{i}(:,l))),'UniformValues',false)];
        variableY.(strcat('Y',strcat(regexprep(num2str(O{i}(:,1:l-1)),' ',''),num2str(O{i}(:,l))))) = sdpvar(n,p);
       end
   end
end
q=fieldnames(variableQ); %extract the names of all elements in the structure variable Q
%% Step 21 of the algorithm
L={};
r=0;
for i=1:length(q)
   for sigma=adm(str2num(q{i}(end)))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
      if ~isKey(delta,[q{i},sigma])
       L=[L;{q{i},sigma}];
      end
    end
end
%% Steps 22 till the of the algorithm
if ~isempty(L)
 for i=1:size(L,1) 
  S={};
  x=[num2str(str2num(L{i,1}(2:end)))-'0',L{i,2}];
  for j=1:length(x)
   if(find(cellfun(@(O) isequal(O, x(:,j:length(x))), O)))
    delta=[delta;containers.Map([L{i,1},L{i,2}],'Q0','UniformValues',false)];
    variableY.(strcat('Y',regexprep(num2str(x),' ',''))) = sdpvar(n,p);
    break;
   end
  end
   for j=1:length(x)
    S=[S;strcat('Q',regexprep(num2str(x(:,j:length(x))),' ',''))];
    int=intersect(S,q);
    if ~isempty(int)
     val=cellfun(@(x) numel(x),int);
     out=int(val==max(val));%argmax
     delta=[delta;containers.Map([L{i,1},L{i,2}],out{:},'UniformValues',false)];
     variableY.(strcat('Y',regexprep(num2str(x),' ',''))) = sdpvar(n,p);
     break;
    end
   end 
 end       
end
%% LMIs corresponding to the first Lyapunov equation
F=[];
for i=1:length(q)
F=[F, variableQ.(q{i})>=eye(n)];
end
%% initial and intermediate and final transitions LMIs
%LMIs corresponding to the second and the third Lyapunov equation
%for each state q and i\in [0,...,m] find \delta(q,i) using the function
%mapa.
for i=1:length(q)
    for j=adm(str2num(q{i}(end)))
    x=q{i}(2:end); %For each state retain only the number part (i.e exclude 'Q')
    jp=num2str(j);tmp=delta([strcat('Q',x),j]);
        if ~strcmp('Q0',tmp(1:2)) %LMI corresponding to the second Lyapunov equation
        F=[F,[variableQ.(tmp),variableQ.(tmp)*A{j}-variableY.(strcat('Y',x,jp))*C{j};(variableQ.(tmp)*A{j}-variableY.(strcat('Y',x,jp))*C{j})',rhoo^2*variableQ.(strcat('Q',x))]];
        else %LMI corresponding to the third Lyapunov equation
        F=[F,[variableQ.(tmp),variableQ.(tmp)*A{j}-variableY.(strcat('Y',x,jp))*C{j};(variableQ.(tmp)*A{j}-variableY.(strcat('Y',x,jp))*C{j})',lambda^2*rhoo^2*variableQ.(strcat('Q',x))]];
        end
    end
end
%% solve the LMIs
su=0;
for i=1:length(q)
su=su+trace(variableQ.(q{i}));
end
optimize(F,su,sdpsettings('solver','sedumi'));
[primalfeas,dualfeas] = check(F);
check(F);
%% Computing the trajectories
% Number of steps
N=20;
traj1_norm=[];
% Number of visits to initial state
kappa=[0];
% Initialization of the trajectory
traj1=[2;1.5;-0.5];
traj1_norm=[traj1_norm,sqrt(traj1(:,end)'*traj1(:,end))];
% V, an array containing the Lyapunov function values initialized by V(Q_0,x_0)
qu=strcat('Q0',num2str(randi([1,m])));
V1=[sqrt(traj1'*value(variableQ.(qu))*traj1)];
% Generate a random switching signal
tmp=adm(str2num(qu(end)));
theta=[tmp(randi([1 4]))];
for i=1:N
    tt=adm(theta(end));
    theta=[theta tt(randi([1 4]))];
end
st={qu}; % the states at each instant
% Update the trajectories array and V at each time step
for i=1:N
    traj1=[traj1 (A{theta(i)}-inv(value(variableQ.(delta([qu,theta(i)]))))*value(variableY.(strcat('Y',qu(2:end),num2str(theta(i)))))*C{theta(i)})*traj1(:,end)]; 
    traj1_norm=[traj1_norm,sqrt(traj1(:,end)'*traj1(:,end))];
    V1=[V1 sqrt(traj1(:,end)'*value(variableQ.(delta([qu,theta(i)])))*traj1(:,end))];
    qu=delta([qu,theta(i)]);
    st=[st qu]; % the states at each instant
    kappa=[kappa kappa(end)+strcmp(qu(1:2),'Q0')];% \kappa(t) is the return index used to compute the accepting rate
end
%% 2nd method: Explicit computation of gains
%% Ensures that Assumption 2 is verified
Q=sdpvar(n,n);
F=[Q>=eye(n)];
for i=1:m
    F=[F,A{i}'*Q*A{i}<=rhoo^2*Q];
end
optimize(F,trace(Q),sdpsettings('solver','sedumi'));
[primalfeas,dualfeas] = check(F);
check(F);
Q=value(Q);
%% gamma computing as in Lemma 3
mx=[];
Po={};
for i=1:length(O)
    Po=[Po,ma(O{i},A,C,rhoo)];
end
y=sdpvar(1,1);
for i=1:length(Po)
F=[F,Po{i}'*Po{i}>=y*Q];
end
optimize(F,-y,sdpsettings('solver','sedumi')); % maximal \gamma
%satisfying the equation in Lemma 3.
[primalfeas,dualfeas] = check(F);
check(F);
y=value(y); %\gamma
%% create states and initialize them
for i=1:m % multiple initial states
variableQQ.(strcat('Q0',num2str(i))) = Q; %initialize Q0 as in eq. (6)
end
for t=1:k-1
x={};
for i=1:t
    x=[x,1:m];
end
% create the states indexed by the powerset of [1,...,m]
  v=allcomb(x{:});
    for i=1:size(v,1)
  variableQQ.(strcat('Q',regexprep(num2str(v(i,:)),' ',''))) = eye(n);
    end
end
%% matrices P computation
% as in eq. 7
E=fieldnames(variableQQ);
for i=1:length(E)
s=E{i};s1=s(1:end-1);s2=str2num(s(end));   
if(strcmp(s1,'Q0')==1) % the m first states are already defined and equal to Q
    continue;
end
if(strcmp(s1,'Q')==1)
s1='Q01';
end
variableQQ.(E{i})=rhoo^2*inv(A{s2}')*variableQQ.(s1)*inv(A{s2})+((rhoo^2)/(y*lambda^2))*inv(A{s2}')*C{s2}'*C{s2}*inv(A{s2}); % eq. 7
end
%% define gains L
% Define the gains using the powerset of [1,...,m].
for t=1:k
x={};
for i=1:t
    x=[x,1:m];
end
  v=allcomb(x{:});
    for i=1:size(v,1)
        if t==1 % replace L with L0 for referring to transitions starting from the initial state
            for j=1:m % can be replaced by j=adm(i)
            variableLL.(strcat('L0',num2str(v(i,end)),num2str(j))) = ones(n,p);
            end
            else
            variableLL.(strcat('L',regexprep(num2str(v(i,:)),' ',''))) = ones(n,p);
        end
    end
end
%% gain computing
% as in eq. 8
E=fieldnames(variableLL);
for i=1:length(E)
    s=E{i};s1=s(2:end-1);s2=str2num(s(end)) ;
    if (strcmp(s(1:2),'L0')==1) % always exceptions for initial state
        s1='01';
    end
variableLL.(E{i})=A{s2}*inv(variableQQ.(strcat('Q',s1)))*C{s2}'*inv(y*lambda^2*eye(size(C{1},1))+C{s2}*inv(variableQQ.(strcat('Q',s1)))*C{s2}');
end
%% Simulation
traj2_norm=[];
traj2=[2;1.5;-0.5]; %trajectories--- 2 for refers to second method
traj2_norm=[traj2_norm,sqrt(traj2(:,end)'*traj2(:,end))];
V2=[sqrt(traj2'*variableQ.(st{1})*traj2)];
for i=1:N
    qu=st{i};
    traj2=[traj2 (A{theta(i)}-variableLL.(strcat('L',qu(2:end),num2str(theta(i))))*C{theta(i)})*traj2(:,end)]; 
    traj2_norm=[traj2_norm,sqrt(traj2(:,end)'*traj2(:,end))];
    V2=[V2 sqrt(traj2(:,end)'*variableQQ.(delta([qu,theta(i)]))*traj2(:,end))];
end
%% temp for plotting
%exchange Q0 with Q_{\epsilon}
to={};
for i=1:size(st,2)
    tmp=cell2mat(st(i));
    if(strcmp(tmp(1:2),'Q0'))
    tmp=strcat('Q',char(949),tmp(end));
     end
to=[to tmp(2:end)];
end
st=to;
%% Plotting stuff
figure;
accept=[];
for i=0:N
    tmp=cell2mat(st(i+1));
    if strcmp(tmp(1),char(949))
        accept=[accept i];
    end
end

subplot1=subplot(2,1,1);% Plot the switching signal \theta_t
stairs([0:N],theta,'LineWidth',2)
hold on
plot(accept,theta(1+accept),'o','MarkerSize',8,'LineWidth',2);
axis([0 N 0.8 m+0.2])
set(subplot1,'FontSize',12);
xlabel('t','FontSize',14)
ylabel('\theta(t)','FontSize',14)

    
% subplot(3,1,2);% Plot the automaton states q_t
% stairs(0:N,categorical(st));
% xlabel('t')
% ylabel('q_{t}')


subplot3=subplot(2,1,2); % Plot the Lyapunov function V(q_t,x_t)
stairs(0:N,traj1_norm(1:N+1),'LineWidth',2)
set(gca,'YScale','log');
hold on;
stairs(0:N,traj2_norm(1:N+1),'LineWidth',2)
yticks([10^(-10) 10^(-5) 1])
set(gca,'YScale','log');
set(subplot3,'FontSize',12);
xlabel('t','FontSize',14)
ylabel('||e(t)||','FontSize',14)
legend("LMI","explicit",'Location','SouthWest')
%{
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
%}
%------------------------------------------------------------------------------

