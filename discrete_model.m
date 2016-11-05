clear; close all;

% 2D Thermal Model for a thin plate with zoned internal heat generation

%% Printing Options
simTitle = 'Nov2';
% Prompts whether to print a new plot when the script is run
% Will plot and save file with some descriptive filename
% Fix .eps printing to tikz printing someday
prompt = 'Should the Temperature Response and Transient Gradient plots be printed? [Y/n]';
prnt = input(prompt,'s');

%% Node Control from .csv
heatZone = csvread('SnakeInput.csv');
heatZoneSize = size(heatZone);
heatZoneV = reshape(heatZone,numel(heatZone),1);
inputNodes = sum(heatZoneV);

n = heatZoneSize(1);
m = heatZoneSize(2);
N = n*m;

%% Manual Node Control
% n = 49;
% m = 49;
% N = n*m;

%% Parameters
tempInit = 0;

x0 = tempInit * ones(n*m,1);

% Properties for somematerial (R)
cp = 840; % J/(kg-K)
p = 1850; % kg/m^3
tmax = 6000; % run time: seconds
ks = 1; % W/(m-K)
L = .16; % length: meter
th = .18; % thickness: meter
Qs = 150; % Value of step power input
area = L^2;
%%%Qos = -15; % Value of step power loss Side
h_t = 50; % Convection coefficient at top
h_b = 10; % Convection coefficient at bottom
Qot = -h_t*(area/N); % Value of step power loss Top
Qob = -h_b*(area/N); % Value of step power loss Bottom
%%%Qot = 1;
%%%Qob = 1;


% Time Step and Vector
dt = tmax/100;
t = 0:dt:tmax;

dx = L/n; % length between each node
a = dx*th; % cross sectional area of node

C = cp*dx*a*p; % thermal capacitance
R = dx/(ks*a); % thermal resistance

%% Setup A Matrix
tic

%% Main diagonal construction: coefficient of i,j element
% Coefficient value for internal elements are inserted at every location, 
% then replaced if the element is on an edge.
Aa = (-(4*R^3)/(C*R^4)+Qot+Qob)*ones(n*m,1); % Value of C_i,j for internal elements
Aa(1:n) = -(3*R^2)/(C*R^3)+Qot+Qob; % Value of C_i,j for left edge replacement
Aa(n*m-(n-1):n*m) = -(3*R^2)/(C*R^3)+Qot+Qob; % Value of C_i,j for right edge 
%    replacement

% This loop modifies the C_i,j value for edge elements along the top and
% bottom of the grid by checking if they are an even multiple of n or are 
% one larger than an even multiple of n.
for i = 1:n*m
    if mod(i,n) == 0
        Aa(i) = -(3*R^2)/(C*R^3)+Qot+Qob; % Bottom row
    elseif mod(i,n) == 1
        Aa(i) = -(3*R^2)/(C*R^3)+Qot+Qob; % Top row
    else
    end
end

% This block adjusts the corners to the correct coefficient values.
Aa(1) = -(2*R)/(C*R^2)+Qot+Qob;
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
%%%Bot = zeros(n*m,1);
%%%Bob = ones(n*m,1);
%%%Bos = zeros(n*m,1);

%% Control for where heat is input by changing B to weight specific nodes.

Bi = 1/(C*R)*heatZoneV;

% Bi = ones(n*m,1); Bi = (1/(C*R))*Bi;
% Bi(ceil(m*n/2)) = 1/(C*R); % Center Element
% Bi(1) = 1/(C*R); % First Element
% Bi(m*n-(n-1):m*n) = 1/(C*R); % All Right Edge
% Bi(1:n) = 1/(C*R); % All Left Edge
% Bi(m*n-(floor(n/2))) = 1/(C*R);
% Bi(n+2) = 1/(C*R);
% Bi(N) = 1/(C*R); % Last Element
% Bi(29:70) = 1/(C*R); % Some Random Elements
% Bi(587:615) = 1/(C*R); % Some Random Elements

% Bos(1:n) = 1/(C*R); % All Left Edge
% Bos(m*n-(n-1):m*n) = 1/(R*C); % All Right Edge
% Bos(n*m) = 1/(C*R); % Last Element
% Bos(1) = 1/(C*R); % First Element

%% Input vectors for SS

% Bot(m*n-(n-1):m*n) = 1/(C*R);

%%%Bot = (1/(C*R))*Bot;
%%%Bob = (1/(C*R))*Bob;
%%%Bos = (1/(C*R))*Bos;
%%%B = [Bi Bot Bob Bos];

B = [Bi];

% B = B * 1/(C*R);

%% Power Input
% Inputs vectors (power input and power loss)
uQs = zeros(length(t),1);
%%%uQot = zeros(length(t),1);
%%%uQob = zeros(length(t),1);
%%%uQos = zeros(length(t),1);

% Power Input Vector (used to set time power is on or off)
uQs(1:length(t),1) = Qs/(nnz(Bi));
%%%uQos(1:length(t)/3,1) = Qos/(nnz(Bi));

% uQs((length(t)-length(t)/3):length(t),1) = Qs/(nnz(Bi));
%uQos((length(t)-length(t)/3):length(t),1) = Qos/(nnz(Bi));
%%%uQot(1:length(t),1) = Qot/N;
%%%uQob(1:length(t),1) = Qob/N;

% Combine inputs for lsim
%%%u = [uQs uQos uQot uQob];
u = [uQs];

% u = Qs * zeros(length(t),1); %step function; for all time, input = Qs
% Qsu = Qs * ones(10,1);
% u(1:10) = Qsu;

%% Output vectors for SS

Cc = speye(m*n);

D = 0;

%% Simulation
tic
sys = ss(A,B,Cc,D);
y = lsim(sys,u,t,x0);
timeSim = toc;

%% Reorder results into element grid and plot/print contours
timecode = num2str(now);

for i = 0:n-1
    for j = 2:length(t)-1
        M(:,i+1,j) = y(j,i*n+1:i*n+n);
    end
end

tempFinal = y(end,:);
maxTemp = max(tempFinal);
minTemp = min(tempFinal);
tenCont = floor((maxTemp-minTemp)/10);

M1 = M(:,:,2);

plots_script

beep