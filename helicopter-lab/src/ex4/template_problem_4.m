% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Flåten

%% Initialization and model definition
clc;
%clear;

init; % Change this to the init file corresponding to your helicopter

% Discrete time system model. x = [lambda r p p_dot]'
delta_t	= 0.25; % sampling time
A1 = [1   delta_t       0                  0              0     0;
      0     1    -delta_t*K_2              0              0     0;
      0     0           1               delta_t           0     0;
      0     0  -delta_t*K_1*K_pp   1 - delta_t*K_1*K_pd   0     0;
      0     0           0                  0              1   delta_t;
      0     0           0                  0           -delta_t*K_3*K_ep   1-delta_t*K_3*K_ed];
B1 = [0     0           0          delta_t*K_1*K_pp       0     0;
      0     0           0                  0              0     delta_t*K_3*K_ep]';
global N mx
% Number of states and inputs
mx = size(A1,2); % Number of states (number of columns in A)
mu = size(B1,2); % Number of inputs(number of columns in B)

% Initial values
x1_0 = pi;                               % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               % e
x6_0 = 0;                               % e_dot
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';           % Initial values

% Time horizon and initialization
N  = 40;                                  % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization
z0(1:mx) = x0;

% Bounds
ul 	    = [-(30*pi/180) -inf]';                  % Lower bound on control
uu 	    = [(30*pi/180) inf]';                   % Upper bound on control

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul(1);                           % Lower bound on state x3
xu(3)   = uu(1);                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu); % hint: gen_constraints
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
P1 = [2 0;
      0 1];                                % Weight on input
Q = gen_q(Q1,P1,N,M);                  % Generate Q, hint: gen_q
c = zeros(N*mx+M*mu,1)';                        % Generate c, this is the linear constant term in the QP

%% Generate system matrixes for linear model
Aeq = gen_aeq(A1,B1,N,mx,mu);             % Generate A, hint: gen_aeq
beq = zeros(N*mx,1);                     % Generate b
beq(1:6) = A1*x0;             
%% Solve QP problem with linear model

fun = @(X) (1/2)*X'*Q*X;

options = optimoptions('fmincon','Algorithm','sqp');



tic
z = fmincon(fun,z0,[],[],Aeq,beq,vlb,vub,@gen_nonlcon, options);
%[z,lambda] = quadprog(Q,c,[],[],Aeq,beq,vlb,vub,z0); % hint: quadprog. Type 'doc quadprog' for more info 
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
  phi1=phi1+Q(i,i)*z(i)*z(i);
  PhiOut(i) = phi1;
end

%% Extract control inputs and states
u_p  = [z(N*mx+1:mu:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution
u_e  = [z(N*mx+2:mu:N*mx+M*mu);z(N*mx+M*mu)];

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution

num_variables = 5/delta_t;
zero_padding = zeros(num_variables*2,1);
unit_padding  = ones(num_variables*2,1);


u_p   = [zero_padding; u_p; zero_padding];
u_e   = [zero_padding; u_e; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];

[obstacle, temp] = gen_nonlcon(z);
obstacle = [0; zero_padding; obstacle; zero_padding] + x5;


t = 0:delta_t:delta_t*(length(u_p)-1);

%Generating optimal input for the simulation
x_star = [x1,x2,x3,x4,x5,x6];
u = [u_p u_e];
optimal_input = timeseries(u, t);
optimal_trajectory = timeseries(x_star, t);

%% Feedback Optimization
Q2 = zeros(mx,mx);
Q2(1,1) = 1;                            % Weight on state x1
Q2(2,2) = 1;                            % Weight on state x2
Q2(3,3) = 1;                            % Weight on state x3
Q2(4,4) = 1;                            % Weight on state x4
Q2(5,5) = 20;                            % Weight on state x5
Q2(6,6) = 2;                            % Weight on state x6
R = [1 0;                             % Weight on input
     0 0.5];

[K,S,e] = dlqr(A1,B1,Q2,R);

t = 0:delta_t:delta_t*(length(u_e)-1);

figure(2)
subplot(511)
stairs(t,u_e),grid
ylabel('u_e')
subplot(512)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
% subplot(513)
% plot(t,x2,'m',t,x2','mo'),grid
% ylabel('r')
subplot(513)
plot(t,obstacle,'m',t,x5','mo'),grid
%plot(t,obstacle,'m',t,obstacle','mo'),grid
xlabel('tid (s)'),ylabel('e')
subplot(514)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(515)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')


