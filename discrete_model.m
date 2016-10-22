clear all; close all;

%% Parameters

% Currently a problem with using even number. Only use odd sized grid.
n = 201;
m = 201;
N = n*m;

x = zeros(n*m,1);

Bi = zeros(n*m,1);
Bot = ones(n*m,1);
Bob = ones(n*m,1);

% Properties for Copper
cp = 390; % J/(kg-K)
p = 8920; % kg/m^3
tmax = 1000; % run time: seconds
ks = 40; % W/(m-K)
L = .05; % length: meter
th = 0.8e-3; % thickness: meter
Qs = 7; % Value of step power input
Qot = -0.1; % Value of step power loss Top
Qob = -0.01; % Value of step power loss Bottom

dt = tmax/1000;
t = 0:dt:tmax;

uQs = Qs * zeros(length(t),1); %step function; for all time, input = Qs
uQot = Qot * ones(length(t),1);
uQob = Qob * ones(length(t),1);

uQs(1:length(t),1) = Qs;

u = [uQs uQot uQob];

% u = Qs * zeros(length(t),1); %step function; for all time, input = Qs
% Qsu = Qs * ones(10,1);
% u(1:10) = Qsu;

dx = L/n; % length between each node
a = dx*th; % cross sectional area of node

C = cp*dx*a*p; % thermal capacitance
R = dx/(ks*a); % thermal resistance

%% Setup A Matrix
tic

Aa = (-(4*R^3)/(C*R^4))*ones(n*m,1);
Aa(1:n) = -(3*R^2)/(C*R^3);
Aa(n*m-(n-1):n*m) = -(3*R^2)/(C*R^3);
for i = 1:n*m
    if mod(i,n) == 0
        Aa(i) = -(3*R^2)/(C*R^3);
    elseif mod(i,n) == 1
        Aa(i) = -(3*R^2)/(C*R^3);
    else
    end
end

Ab = (1/(C*R))*ones(n*m-n,1);

Ac = (1/(C*R))*ones(n*m-1,1);
for i = 1:n*m-1
    if mod(i,n) == 0
        Ac(i) = 0;
    else
    end
end

Ad = (1/(C*R))*ones(n*m-1,1);
for i = 2:n*m-1
    if mod(i,n) == 1
        Ad(i-1) = 0;
    else
    end
end

v = zeros(n,1);
Aba = [v;Ab];
Abb = [Ab;v];
v = [0];
Ac = [Ac;v];
Ad = [v;Ad];
% Aa(1) = -1/(C*R);
% Aa(n*m) = -1/(C*R);

A = [Abb Ac Aa Ad Aba];
v = [-n -1 0 1 n]';

A = spdiags(A,v,N,N);

timeDiag = toc;

% Bi(ceil(m*n/2)) = 1/(C*R); % Center Element
% Bi(1:n) = 1/(C*R); % All Left Edge
% Bi(n) = 1/(C*R);
% Bi(m*n-(n-1):m*n) = 1/(C*R); % All Right Edge
% Bi(m*n) = 1/(C*R);
% Bi(m*n-(floor(n/2))) = 1/(C*R);
% Bi(n+2) = 1/(C*R);
% Bi(n) = 1/(C*R);
Bi(1654:2278) = 1/(C*R);

Bot = (1/(C*R))*Bot;

Bob = (1/(C*R))*Bob;

B = [Bi Bot Bob];

% B = B * 1/(C*R);

Cc = speye(m*n);

D = 0;

%% Simulation
tic
sys = ss(A,B,Cc,D);
y = lsim(sys,u,t);
timeSim = toc;

%% Reorder results into element grid and plot contours

for i = 0:n-1
    for j = 2:length(t)-1
        M(:,i+1,j) = y(j,i*n+1:i*n+n);
    end
end

M1 = M(:,:,2);

figure

g = 4;
for j = 1:g
    subplot(g/2,2,j)
    contourf(M(:,:,j*floor(length(t)/g)),10,'ShowText','on');
    xlabel('x position (node)'); ylabel('y position (node)')
    lbl = num2str(t(j*floor(length(t)/g)));
    lbl = [lbl,' seconds'];
    ttl = ['Time = ',lbl]; title(ttl);
end

% print('/Users/jonathon/Documents/MME Controls Theory/latex/math_model/figures/contour_trans','-depsc')

figure
contourf(M(:,:,length(t)-1),10,'ShowText','on');
grid on
xlabel('x position (node)'); ylabel('y position (node)')
lbl = num2str(t(end));
lbl = [lbl,' seconds'];
ttl = ['Time = ',lbl]; title(ttl);

figure
plot(t,y(:,1),'LineWidth',2)
hold on
plot(t,y(:,n),'-.','LineWidth',2)
plot(t,y(:,m*n),'LineWidth',2)
plot(t,y(:,m*n-(n-1)),'--','LineWidth',2)
plot(t,y(:,ceil(m*n/2)),'LineWidth',2)
xlabel('Time (seconds)'); ylabel('Temperature (Degrees Celcius');
legend('Bottom Left','Top Left','Top Right','Bottom Right','Center Element')
% print('/Users/jonathon/Documents/MME Controls Theory/latex/math_model/figures/time_temp','-depsc')

timeSim/60;

if timeSim <= 60
    disp(timeSim)
else
    disp(timeSim/60)
end
beep