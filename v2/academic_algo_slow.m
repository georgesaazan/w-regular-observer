%% This produces the results of the Numerical Example in the paper
% This ouputs a figure corresponding to the fast switching signal, in order
% to get the second figure corresponding to the slow swtiching signal,
% check line 125.
%% define system and fix parameters
clear all;clc;
rhoo=1.5;lambda=0.1; %fix \rho and \lambda
k=3; %fik k the maximal length of observable sequences
A1=eye(3);A2=1.5*[0,0,1;0,1,0;1,0,0]; %define the system
C1=[1,0,0];C2=[0,1,1];
A={A1,A2};
C={C1,C2};
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
    if rank(ma(v(i,:),A,C))==n %check if a combination is an observable sequence
        O=[O;v(i,:)];%add the observable sequence to I
    end
end
end
OO=O;
len_O=length(O);
for i=1:len_O
    S={};P={};
    if(size(O{i},2)~=1)
    for j=1:size(O{i},2)-1
    P=[P;O{i}(:,1:j)];
    end
    for j=size(O{i},2):-1:2
    S=[S;O{i}(:,j:end)];
    end
    if(sum(find(cellfun(@(O) isequal(O, P{:}), O)))~=0||sum(find(cellfun(@(O) isequal(O, S{:}), O)))~=0)
        O{i}={}; 
    end
    end
end
%% Steps 1 till 9 of the algorithm
I=1:m;
variableQ.('Q0') = sdpvar(n,n); %initial state, step 1
delta=containers.Map;   
for sigma=I %step 2 till step 9 
    variableY.(strcat('Y0',num2str(sigma))) = sdpvar(n,p);
if(find(cellfun(@(O) isequal(O, [sigma]), O)))
     delta=[delta;containers.Map(['Q0',sigma],'Q0','UniformValues',false)];
else
     variableQ.(strcat('Q',num2str(sigma))) = sdpvar(n,n);
     delta=[delta;containers.Map(['Q0',sigma],strcat('Q',num2str(sigma)),'UniformValues',false)];
end
end
%% Steps 10 till 20 of the algorithm
for i=1:length(O)
    l=size(O{i},2);
   if l>=2
       if l==2
         delta=[delta;containers.Map([strcat('Q',O{i}(:,1)),O{i}(:,2)],'Q0','UniformValues',false)];
        variableY.(strcat('Y',regexprep(num2str(O{i}(:,1:l)),' ',''))) = sdpvar(n,p);
       else
           for j=2:l-1
           variableQ.(strcat('Q',regexprep(num2str(O{i}(:,1:j)),' ',''))) = sdpvar(n,n);
           delta=[delta;containers.Map([strcat('Q',regexprep(num2str(O{i}(:,1:j-1)),' ','')),O{i}(:,j)],strcat('Q',regexprep(num2str(O{i}(:,1:j)),' ','')),'UniformValues',false)];
           variableY.(strcat('Y',strcat(regexprep(num2str(O{i}(:,1:j-1)),' ',''),num2str(O{i}(:,j))))) = sdpvar(n,p);
           end
           delta=[delta;containers.Map([strcat('Q',regexprep(num2str(O{i}(:,1:l-1)),' ','')),O{i}(:,l)],'Q0','UniformValues',false)];
        variableY.(strcat('Y',strcat(regexprep(num2str(O{i}(:,1:l-1)),' ',''),num2str(O{i}(:,l))))) = sdpvar(n,p);
       end
   end
end
q=fieldnames(variableQ); %extract the names of all elements in the structure variable Q
%% Step 21 of the algorithm
L={};
r=0;
for i=1:length(q)
   for sigma=1:m                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
      try
       delta([q{i},sigma]);
      catch err
       r=r+1;
       L{r}={q{i},sigma};
      end
    end
end
%% Steps 22 till the of the algorithm
if ~isempty(L)
 for i=1:length(L)  
  S={};
  x=[num2str(str2num(L{i}{1}(2:end)))-'0',L{i}{2}];
  for j=1:length(x)
   if(find(cellfun(@(O) isequal(O, x(:,j:length(x))), O)))
    delta=[delta;containers.Map([L{i}{1},L{i}{2}],'Q0','UniformValues',false)];
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
     delta=[delta;containers.Map([L{i}{1},L{i}{2}],out{:},'UniformValues',false)];
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
    for j=1:m
    x=q{i}(2:end); %For each state retain only the number part (i.e exclude 'Q')
    jp=num2str(j);
        if (strcmp('Q0',delta([strcat('Q',x),j]))~=1) %LMI corresponding to the second Lyapunov equation
        F=[F,[variableQ.(delta([strcat('Q',x),j])),variableQ.(delta([strcat('Q',x),j]))*A{j}-variableY.(strcat('Y',x,jp))*C{j};(variableQ.(delta([strcat('Q',x),j]))*A{j}-variableY.(strcat('Y',x,jp))*C{j})',rhoo^2*variableQ.(strcat('Q',x))]];
        else %LMI corresponding to the third Lyapunov equation
        F=[F,[variableQ.(delta([strcat('Q',x),j])),variableQ.(delta([strcat('Q',x),j]))*A{j}-variableY.(strcat('Y',x,jp))*C{j};(variableQ.(delta([strcat('Q',x),j]))*A{j}-variableY.(strcat('Y',x,jp))*C{j})',lambda^2*rhoo^2*variableQ.(strcat('Q',x))]];
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
N=40;
% Number of visits to initial state
kappa=[0];
% Initialization of the trajectory
traj1=[2;1.5;-0.5];
% V, an array containing the Lyapunov function values initialized by V(Q_0,x_0)
qu='Q0';
V1=[sqrt(traj1'*value(variableQ.(qu))*traj1)];
% Generate a random switching signal
theta=[randi([1 2])];
% for the slow switching signal comment out line 161 and uncomment out line
% 162.
%for i=1:N/4 theta=[theta O{randi([1,length(O)])}]; end
N=150;theta(1:30)=2;theta(31:40)=1;theta=[theta,theta,theta];theta(end:N+1)=2;
st={'Q0'};
% Update the trajectories array and V at each time step
for i=1:N
    traj1=[traj1 (A{theta(i)}-inv(value(variableQ.(delta([qu,theta(i)]))))*value(variableY.(strcat('Y',qu(2:end),num2str(theta(i)))))*C{theta(i)})*traj1(:,end)]; 
    V1=[V1 sqrt(traj1(:,end)'*value(variableQ.(delta([qu,theta(i)])))*traj1(:,end))];
    qu=delta([qu,theta(i)]);
    st=[st qu]; % the states at each instant
    kappa=[kappa kappa(end)+strcmp(qu,'Q0')];% \kappa(t) is the return index used to compute the accepting rate
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
    Po=[Po,mb(O{i},A,C,rhoo)];
end
y=sdpvar(1,1);
for i=1:length(Po)
F=[F,Po{i}>=y*Q];
end
optimize(F,-y,sdpsettings('solver','sedumi')); % maximal \gamma
%satisfying the equation in Lemma 3.
[primalfeas,dualfeas] = check(F);
check(F);
y=value(y); %\gamma
%% create states and initialize them
variableQQ.('Q0') = Q; %initialize Q0 as in eq. (6)
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
for i=2:length(E)
s=E{i};s1=s(1:end-1);s2=str2num(s(end)) ;   
if(strcmp(s1,'Q')==1)
s1='Q0';
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
            variableLL.(strcat('L0',regexprep(num2str(v(i,:)),' ',''))) = ones(n,p);
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
variableLL.(E{i})=A{s2}*inv(variableQQ.(strcat('Q',s1)))*C{s2}'*inv(y*lambda^2*eye(size(C{1},1))+C{s2}*inv(variableQQ.(strcat('Q',s1)))*C{s2}');
end
%% Simulation
traj2=[2;1.5;-0.5]; %trajectories--- 2 for refers to second method
V2=[sqrt(traj2'*variableQ.(st{1})*traj2)];
for i=1:N
    qu=st{i};
    traj2=[traj2 (A{theta(i)}-variableLL.(strcat('L',qu(2:end),num2str(theta(i))))*C{theta(i)})*traj2(:,end)]; 
    V2=[V2 sqrt(traj2(:,end)'*variableQQ.(delta([qu,theta(i)]))*traj2(:,end))];
end
%% temp for plotting
%exchange Q0 with Q_{\epsilon}
to={};
for i=1:size(st,2)
    tmp=cell2mat(st(i));
    if(strcmp(tmp,'Q0'))
    tmp=strcat('Q',char(949)); % 949 for \epsilon
     end
to=[to tmp(2:end)];
end
st=to;
%% Plotting stuff
figure;
% Plot 4*1 subfigures on the same figure
subplot1=subplot(3,1,1);% Plot the switching signal \theta_t
stairs([0:N],theta,'LineWidth',2)
axis([0 N 0.9 m+0.1])
set(subplot1,'FontSize',12);
xlabel('t','FontSize',14)
ylabel('\theta(t)','FontSize',14)


subplot2=subplot(3,1,2);% Plot the automaton states q_t
stairs(0:N,categorical(st),'LineWidth',2);
accept=(categorical(st)==categorical(st(1)));
time=[0:N];
hold on
plot(time(accept),categorical(st(accept)),'o','MarkerSize',8,'LineWidth',2);
set(subplot2,'FontSize',12);
xlabel('t','FontSize',14)
ylabel('q(t)','FontSize',14)

% subplot(4,1,3);% Plot the switching signal \theta_t
% stairs([0:N],kappa./[0:N])
% xlabel('t')
% ylabel('k(t)/t')

subplot3=subplot(3,1,3); % Plot the Lyapunov function V(q_t,x_t)
stairs(0:N,V1(1:N+1),'LineWidth',2)
set(gca,'YScale','log');
hold on;
stairs(0:N,V2(1:N+1),'LineWidth',2)
yticks([10^(-5) 1 10^(5) 10^10])
set(gca,'YScale','log');
set(subplot3,'FontSize',12);
xlabel('t','FontSize',14)
ylabel('||e(t)||','FontSize',14)
legend("LMI","explicit",'Location','NorthWest')
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

