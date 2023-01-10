%% This produces the results of Numerical Example 1 in the paper
%% Fix k and define system
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
%% Finding the observable sequences
if (k<n) error('No observable sequences for the specified k')
end
I=cell(1,k-n+1); %initialization of a cell that contains the observable sequences of length n to k separately
for t=n:1:k % Compute observable sequences
l=[];
x={};
for i=1:t
    x=[x,1:m];%intermediate variable for the next step
end
v=allcomb(x{:});%create all possible combinations of switching sequences of length up to t, n\le t\le k
for i=1:length(v)
    if rank(ma(v(i,:),A,C))==n %check if a combination is an observable sequence
        I{t-n+1}=[I{t-n+1};v(i,:)]; %add the observable sequence to I
    end
end
% The switching is constrained, therefore we should remove all inadmissible
% sequences
for i=1:size(I{t-n+1},1)
    for j=1:size(I{t-n+1},2)-1
        if(ismember(I{t-n+1}(i,j+1),adm(I{t-n+1}(i,j)))==0) %condition to check if the sequence is inadmissible
            l=[l;i]; %get the index of the inadmissible sequence
        end
    end
end
l=unique(l,'rows'); %remove the recurrent indexes
if (isempty(I{t-n+1})==0)
I{t-n+1}(l,:)=[]; % remove inadmissible observable sequences
end
end
%minimal observable sequences
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
%here we have m initial states Q01,...,Q08
for i=1:m
    variableQ.(strcat('Q0',num2str(i))) = sdpvar(n,n);
end
for i=1:m
   variableQ.(strcat('Q',num2str(i))) = sdpvar(n,n);  %states Q_1 ... Q_m
   for j=adm(i)
   variableY.(strcat('Y0',num2str(i),num2str(j))) = sdpvar(n,p); %transitions from Q01,...,Q08
   end
end
for i=1:length(I)
for j=1:size(I{i},1)
for l=2:size(I{i},2)
    % for an observable sequence (i_1 ... i_n), define the state Q_{12},
    % the transition from Q_{1} to Q_{12} ... the state Q_{12...n-1} and the
    % transition from Q_{12...n-2} to Q_{12...n-1}
    variableQ.(strcat('Q',regexprep(num2str(I{i}(j,1:l-1)),' ',''))) = sdpvar(n,n);
    q=fieldnames(variableQ); %extract the names of all elements in the structure variable Q
    for b=adm(I{i}(j,l-1)) % consider constraints for transitions
    variableY.(strcat('Y',regexprep(num2str(I{i}(j,1:l-1)),' ',''),num2str(b))) = sdpvar(n,p);
    end
end
end
end
%% fake states corresponding to initial state
% for each observable sequence (i_1,...,i_n) create a fake state
% Q_{1...n} which will be used later.
qt={};
for i=1:k-n+1
    for j=1:size(I{i},1)
    qt=[qt;strcat('Q',regexprep(num2str(I{i}(j,:)),' ',''))];
    end
end
%% define transitions
% Ensure that all transitions are created
% Useful when there is an index that is not at the beginning of any observable sequence
for i=m+1:2*m
for j=adm(str2num(q{i}(2:end)))%consider constraints
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
%mapa2.
for i=1:length(q)
    for j=adm(str2num(q{i}(end)))
    x=q{i}(2:end); %For each state retain only the number part (i.e exclude 'Q')
    jp=num2str(j);
        if (strcmp('Q0',mapa2(q{i},j,q,qt))~=1) %LMI corresponding to the second Lyapunov equation
        F=[F,[variableQ.(mapa2(q{i},j,q,qt)),variableQ.(mapa2(q{i},j,q,qt))*A{j}-variableY.(strcat('Y',x,jp))*C{j};(variableQ.(mapa2(q{i},j,q,qt))*A{j}-variableY.(strcat('Y',x,jp))*C{j})',rhoo^2*variableQ.(q{i})]];
        else %LMI corresponding to the third Lyapunov equation
      F=[F,[variableQ.(mapa2(q{i},j,q,qt)),variableQ.(mapa2(q{i},j,q,qt))*A{j}-variableY.(strcat('Y',x,jp))*C{j};(variableQ.(mapa2(q{i},j,q,qt))*A{j}-variableY.(strcat('Y',x,jp))*C{j})',lambda^2*rhoo^2*variableQ.(q{i})]];
        end
    end
end
%% Solve the LMIs
su=0;
for i=1:length(q)
su=su+trace(variableQ.(q{i}));
end
optimize(F,su,sdpsettings('solver','sedumi'));
[primalfeas,dualfeas] = check(F);
check(F);
%% Computing the trajectories
N=20;
% Number of steps
% Initialization of the trajectory
traj1=[-1.5;-0.5;2];
% V, an array containing the Lyapunov function values initialized by V(q_0,x_0)
%start at random initial state
X={'Q01','Q02','Q03','Q04','Q05','Q06','Q07','Q08'};
qu=string(X(randi(numel(X))));
qu=convertStringsToChars(qu);
V1=[sqrt(traj1'*value(variableQ.(qu))*traj1)];
% generate a random switching signal while respecting the constraints
tmp=adm(str2num(qu(end)));
theta=[tmp(randi([1 3]))];
for i=1:N
    tt=adm(theta(end));
    theta=[theta tt(randi([1 3]))];
end
st={qu}; % the states at each instant
% Update the trajectories array and V at each time step
for i=1:N
 traj1=[traj1 (A{theta(i)}-(inv(value(variableQ.(mapa2(qu,theta(i),q,qt))))*value(variableY.(strcat('Y',qu(2:end),num2str(theta(i))))))*C{theta(i)})*traj1(:,end)]; 
 V1=[V1 sqrt(traj1(:,end)'*value(variableQ.(mapa2(qu,theta(i),q,qt)))*traj1(:,end))];
 qu=mapa2(qu,theta(i),q,qt);
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
    for(j=1:size(I{i},1))
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
            for j=1:m % can be replaced by adm
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
traj2=[2;1.5;-0.5]; %trajectories--- 2 for refers to second method
V2=[sqrt(traj2'*variableQ.(st{1})*traj2)];
for i=1:N
    qu=st{i};
    traj2=[traj2 (A{theta(i)}-variableLL.(strcat('L',qu(2:end),num2str(theta(i))))*C{theta(i)})*traj2(:,end)]; 
    V2=[V2 sqrt(traj2(:,end)'*variableQQ.(mapa2(qu,theta(i),q,qt))*traj2(:,end))];
end
%% temp for plotting
%exchange Q0 with Q_{\epsilon}
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
% Plot 3*1 subfigures on the same figure
subplot(3,1,1);% Plot the switching signal \theta_t
stairs([0:N],theta)
axis([0 N 0.9 m+0.1])
xlabel('t')
ylabel('\theta(t)')

subplot(3,1,2); % Plot the automaton states q_t
stairs(0:N,categorical(st));
xlabel('t')
ylabel('q_{t}')

subplot(3,1,3); % Plot the Lyapunov function V(q_t,x_t)
stairs(0:N,log(V1(1:N+1))/log(10))
hold on;
stairs(0:N,log(V2(1:N+1))/log(10))
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

