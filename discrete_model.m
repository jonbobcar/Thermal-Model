clear all; close all;

% 2D Thermal Model for a thin plate with zoned internal heat generation

%% Printing Options
simTitle = 'test';
% Prompts whether to print a new plot when the script is run
% Will plot and save file with some descriptive filename
% Fix .eps printing to tikz printing someday
prompt = 'Should the Temperature Response and Transient Gradient plots be printed? [Y/n]';
prnt = input(prompt,'s');

%% Node Control from .csv
% heatZone = csvread('SnakeInput.csv');
% heatZoneSize = size(heatZone);
% heatZoneV = reshape(heatZone,numel(heatZone),1);
% inputNodes = sum(heatZoneV);
% 
% n = heatZoneSize(1);
% m = heatZoneSize(2);
% N = n*m;

%% Manual Node Control
n = 25;
m = 25;
N = n*m;

%% Parameters
tempInit = 0;

x0 = tempInit * ones(n*m,1);

% Properties for Copper
cp = 840; % J/(kg-K)
p = 1850; % kg/m^3
tmax = 10000; % run time: seconds
ks = 100; % W/(m-K)
L = .16; % length: meter
th = .18; % thickness: meter
Qs = 5; % Value of step power input
area = L^2;
Qos = 5/n; % Value of step power loss Side
h_t = 1; % Convection coefficient from top
h_b = 1; % Convection coefficient from bottom
Qot = 1; % Value of step power loss Top
Qob = 1; % Value of step power loss Bottom


% Time Step and Vector
dt = tmax/100;
t = 0:dt:tmax;

dx = L/n; % length between each node
a = dx*th; % cross sectional area of node

C = cp*dx*a*p; % thermal capacitance
R = dx/(ks*a); % thermal resistance

% Inputs vectors (power input and power loss)
uQs = zeros(length(t),1);
uQot = zeros(length(t),1);
uQob = zeros(length(t),1);
uQos = zeros(length(t),1);

% Power Input Vector (used to set time power is on or off)
uQs(1:length(t),1) = Qs;
uQos(1:length(t),1) = Qos;
uQot(1:length(t),1) = Qot;
uQob(1:length(t),1) = Qob;

% Combine inputs for lsim
u = [uQs];

% u = Qs * zeros(length(t),1); %step function; for all time, input = Qs
% Qsu = Qs * ones(10,1);
% u(1:10) = Qsu;

%% Setup A Matrix
tic

%% Main diagonal construction: coefficient of i,j element
% Coefficient value for internal elements are inserted at every location, 
% then replaced if the element is on an edge.
Aa = (-(4*R^3)/(C*R^4))*ones(n*m,1)+Qob+Qot; % Value of C_i,j for internal elements
Aa(1:n) = -(3*R^2)/(C*R^3); % Value of C_i,j for left edge replacement
Aa(n*m-(n-1):n*m) = -(3*R^2)/(C*R^3); % Value of C_i,j for right edge replacement

% This loop modifies the C_i,j value for edge elements along the top and
% bottom of the grid by checking if they are an even multiple of n or are 
% one larger than an even multiple of n.
for i = 1:n*m
    if mod(i,n) == 0
        Aa(i) = -(3*R^2)/(C*R^3); % Bottom row
    elseif mod(i,n) == 1
        Aa(i) = -(3*R^2)/(C*R^3); % Top row
    else
    end
end

% This block adjusts the corners to the correct coefficient values.
Aa(1) = -(2*R)/(C*R^2);
Aa(n) = Aa(1);
Aa(n*m) = Aa(1);
Aa(n*m-(n-1)) = Aa(1);

%% Auxiliary diagonal construction: coefficients of adjacent elements

% Ab will map to two diagonals directly next to the main diagonal and
% represents resistance between C_i,j and C_i,j(+/-)1 (elements east and
% west). The A Matrix is constructed by first going down the left column
% and finishing by going down the right column, so offsetting the diagonal
% by n rows takes care of removing resistances where there is no element to
% the east or west.
Ab = (1/(C*R))*ones(n*m-n,1);

% Ac will map to a diagonal which represents resistance between C_i,j and
% C_i-1,j (element south)
Ac = (1/(C*R))*ones(n*m-1,1);

% This loop removes elements from Ac at the bottom edge of the grid array,
% where there is no thermal element to the south.
for i = 1:n*m-1
    if mod(i,n) == 0
        Ac(i) = 0;
    else
    end
end

% Ad will map to a diagonal which represents resistance between C_i,j and
% C_i+1,j (element north)
Ad = (1/(C*R))*ones(n*m-1,1);

% This loop removes elements from Ad at the top edge of the grid array,
% where there is no thermal element to the north.
for i = 2:n*m-1
    if mod(i,n) == 1
        Ad(i-1) = 0;
    else
    end
end

%% Combine state variable coefficient values into sparse A matrix.
v = zeros(n,1); % zeroes required to offset east and west edges
Aba = [v;Ab]; % elements to the east
Abb = [Ab;v]; % elements to the west
v = [0]; % zero required to offset north and south edges
Ac = [Ac;v]; % elements to the south
Ad = [v;Ad]; % elements to the north
% Aa(1) = -1/(C*R);
% Aa(n*m) = -1/(C*R);

% Vectors required to create sparse A matrix
A = [Abb Ac Aa Ad Aba];
v = [-n -1 0 1 n]';

A = spdiags(A,v,N,N);

A_full = full(A);

timeDiag = toc;

% Input Vectors (Bi = Heat in; Bot,Bob,Bos = Heat out (top,bottom,side))
Bi = zeros(n*m,1);
% Bot = zeros(n*m,1);
% Bob = ones(n*m,1);
% Bos = zeros(n*m,1);

%% Control for where heat is input by changing B to weight specific nodes.

% Bi = 1/(C*R)*heatZoneV;

% Bi = ones(n*m,1); Bi = (1/(C*R))*Bi;
Bi(ceil(m*n/2)) = 1/(C*R); % Center Element
% Bi(1) = 1/(C*R); % First Element
% Bi(m*n-(n-1):m*n) = 1/(C*R); % All Right Edge
% Bi(n) = 1/(C*R);
% Bi(m*n-(floor(n/2))) = 1/(C*R);
% Bi(n+2) = 1/(C*R);
% Bi(n) = 1/(C*R);
% Bi(1654:2278) = 1/(C*R);

% Bos(1:n) = 1/(C*R); % All Left Edge
% Bos(m*n-(n-1):m*n) = 1/(R*C); % All Right Edge
% Bos(n*m) = 1/(C*R); % Last Element
% Bos(1) = 1/(C*R); % First Element

%% Input vectors for SS

% Bot(m*n-(n-1):m*n) = 1/(C*R);

% Bot = (1/(C*R))*Bot;
% Bob = (1/(C*R))*Bob;
% Bos = (1/(C*R))*Bos;
B = [Bi];

% B = B * 1/(C*R);

%% Output vectors for SS

Cc = speye(m*n);

D = 0;

%% Simulation
tic
sys = ss(A,B,Cc,D);
y = lsim(sys,u,t,x0);
timeSim = toc;

%% Reorder results into element grid and plot/print contours

for i = 0:n-1
    for j = 2:length(t)-1
        M(:,i+1,j) = y(j,i*n+1:i*n+n);
    end
end

tempFinal = y(end,:);
maxTemp = max(tempFinal);
minTemp = min(tempFinal);
tenCont = floor((maxTemp-minTemp)/10);
tenCont = 10;

M1 = M(:,:,2);

figure

g = 4;

for j = 1:g
    subplot(g/2,2,j)
    contourf(M(:,:,j*floor(length(t)/g)),'ShowText','on','LevelStepMode','manual','LevelStep',tenCont);
    xlabel('x position (node)'); ylabel('y position (node)')
    lbl = num2str(t(j*floor(length(t)/g)));
    lbl = [lbl,' seconds'];
    ttl = ['Time = ',lbl]; title(ttl);
    set(gca,'Ydir','reverse');
    axis equal
end
pathStr = '/Users/jonathon/Documents/Thesis/thermalmodel/latex/figures/CtrTnst';
tempInitStr = [',it' num2str(tempInit)];
powerInputStr = [',pi' num2str(Qs)];
pathStr = strcat(pathStr,simTitle);
if strcmp(prnt,'Y') == 1
    print(pathStr,'-depsc')
else
end

figure
contourf(M(:,:,length(t)-1),10,'ShowText','on');
grid on
xlabel('x position (node)'); ylabel('y position (node)')
lbl = num2str(t(end));
lbl = [lbl,' seconds'];
ttl = ['Time = ',lbl]; title(ttl);
set(gca,'Ydir','reverse');
axis equal

figure
plot(t,y(:,1),'LineWidth',2)
hold on
plot(t,y(:,n),'-.','LineWidth',2)
plot(t,y(:,m*n),'LineWidth',2)
plot(t,y(:,m*n-(n-1)),'--','LineWidth',2)
plot(t,y(:,ceil(m*n/2)),'LineWidth',2)
xlabel('Time (seconds)'); ylabel('Temperature (Degrees Celcius)');
legend('Bottom Left','Top Left','Top Right','Bottom Right','Center Element')
pathStr = '/Users/jonathon/Documents/Thesis/thermalmodel/latex/figures/TmpRspn';
pathStr = strcat(pathStr,simTitle);
if strcmp(prnt,'Y') == 1
    print(pathStr,'-depsc')
else
end

timeSim/60;
timeSimSec = num2str(timeSim,3);
timeSimMin = num2str(timeSim/60,3);
if timeSim <= 60
    disp(['Simulation Run time = ',timeSimSec,' seconds'])
else
    disp(['Simulation Run time = ',timeSimMin,' minutes'])
end
beep