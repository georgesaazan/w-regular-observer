%% This produces the results of the Numerical Example in the paper
% It is useful to edit theta on line 122 to produce the 2 figures
%% define system and fix parameters
clear all;clc;
rhoo=1.5;lambda=0.2; %fix \rho and \lambda
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
I=cell(1,k-n+1); %initialization of a cell that contains the observable sequences of length n to k separately
for t=n:1:k
x={};
for i=1:t
    x=[x,1:m];%intermediate variable for the next step
end
v=allcomb(x{:});%create all possible combinations of switching sequences of length up to t, n\le t \le k
for i=1:length(v)
    if rank(ma(v(i,:),A,C))==n %check if a combination is an observable sequence
        I{t-n+1}=[I{t-n+1};v(i,:)];%add the observable sequence to I
    end
end
end
if length(I)>=2 %Compute the minimal observable sequences when k>n
for i=length(I):-1:2
    for j=i-1:-1:1
   I{i}=remov(I{j},I{i});%remove the non-minimal observable sequences from the observable sequences of length up to k
    end
end
end
%check if there isn't observable sequences up to the specified k to end the program
s=0; 
for i=1:length(I)
s=s+isempty(I{i});
end
if(s==length(I)) error('No observable sequences for the specified k')
end
%% create states based on minimal observable sequences to build the LMIs
% To build our LMIs we must define a matrix Q per state and Y per transition
variableQ.('Q0') = sdpvar(n,n); %initial state
for i=1:m
   variableQ.(strcat('Q',num2str(i))) = sdpvar(n,n);  %states Q_1 ... Q_m
   variableY.(strcat('Y0',num2str(i))) = sdpvar(n,p); %transitions from Q_0 to Q_1, ... ,Q_m
end
for i=1:length(I)
for j=1:size(I{i},1)
for l=2:size(I{i},2)
    % for an observable sequence (i_1 ... i_n), define the state Q_{12},
    % the transition from Q_{1} to Q_{12} ... the state Q_{12...n-1} and the
    % transition from Q_{12...n-2} to Q_{12...n-1}
    variableQ.(strcat('Q',regexprep(num2str(I{i}(j,1:l-1)),' ',''))) = sdpvar(n,n);
    q=fieldnames(variableQ); %extract the names of all elements in the structure variable Q
    for b=1:m
      variableY.(strcat('Y',regexprep(num2str(I{i}(j,1:l-1)),' ',''),num2str(b))) = sdpvar(n,p);
    end
end
end
end
%% fake states corresponding to initial state
% for each observable sequence (i_1,...,i_n) create a fake state
% Q_{1...n} which will be used later.
qt={};
for i=1:length(I)
    for j=1:size(I{i},1)
    qt=[qt;strcat('Q',regexprep(num2str(I{i}(j,:)),' ',''))];
    end
end
%% define transitions
% Ensure that all transitions are created
% Useful when there is an index that is not at the beginning of any observable sequence
for i=2:m+1
for j=1:m
variableY.(strcat('Y',q{i}(2:end),num2str(j))) = sdpvar(n,p);
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
        if (strcmp('Q0',mapa(strcat('Q',x),j,q,qt))~=1) %LMI corresponding to the second Lyapunov equation
        F=[F,[variableQ.(mapa(strcat('Q',x),j,q,qt)),variableQ.(mapa(strcat('Q',x),j,q,qt))*A{j}-variableY.(strcat('Y',x,jp))*C{j};(variableQ.(mapa(strcat('Q',x),j,q,qt))*A{j}-variableY.(strcat('Y',x,jp))*C{j})',rhoo^2*variableQ.(strcat('Q',x))]];
        else %LMI corresponding to the third Lyapunov equation
      F=[F,[variableQ.(mapa(strcat('Q',x),j,q,qt)),variableQ.(mapa(strcat('Q',x),j,q,qt))*A{j}-variableY.(strcat('Y',x,jp))*C{j};(variableQ.(mapa(strcat('Q',x),j,q,qt))*A{j}-variableY.(strcat('Y',x,jp))*C{j})',lambda^2*rhoo^2*variableQ.(strcat('Q',x))]];
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
N=200;
% Initialization of the trajectory
traj1=[2;1.5;-0.5];
% V, an array containing the Lyapunov function values initialized by V(Q_0,x_0)
V1=[sqrt(traj1'*value(variableQ.('Q0'))*traj1)];
qu='Q0';
% Generate a random switching signal
theta=[];
% for the fast switching signal comment out the second line and vise versa
% for the slow switching signal.
theta=[randi([1 m])]; for i=1:N temp =randi([1 m]); theta=[theta ,temp]; end
%theta(1:40)=2;theta(41:60)=1;theta=[theta,theta,theta];theta(end:201)=2;
st={'Q0'};
% Update the trajectories array and V at each time step
for i=1:N
    traj1=[traj1 (A{theta(i)}-(inv(value(variableQ.(mapa(qu,theta(i),q,qt))))*value(variableY.(strcat('Y',qu(2:end),num2str(theta(i))))))*C{theta(i)})*traj1(:,end)]; 
    V1=[V1 sqrt(traj1(:,end)'*value(variableQ.(mapa(qu,theta(i),q,qt)))*traj1(:,end))];
    qu=mapa(qu,theta(i),q,qt);
    st=[st qu]; % the states at each instant
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
for i=1:length(I)
    for j=1:size(I{i},1)
    Po=[Po,mb(I{i}(j,:),A,C,rhoo)];
    end
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
    V2=[V2 sqrt(traj2(:,end)'*variableQQ.(mapa(qu,theta(i),q,qt))*traj2(:,end))];
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
% Plot 3*1 subfigures on the same figure
subplot(3,1,1);% Plot the switching signal \theta_t
stairs([0:N],theta,"color",'k')
axis([0 N 0.9 m+0.1])
xlabel('t')
ylabel('\theta(t)')

subplot(3,1,2);% Plot the automaton states q_t
stairs(0:N,categorical(st),"color",'k');
xlabel('t')
ylabel('q_{t}')

subplot(3,1,3); % Plot the Lyapunov function V(q_t,x_t)
stairs(0:N,log(V1(1:N+1))/log(10),'-.',"color",'k')
hold on;
stairs(0:N,log(V2(1:N+1))/log(10),"color",'k')
xlabel('t')
ylabel('log(V(q_{t},x(t)))')
legend("LMI method","explicit construction")
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
