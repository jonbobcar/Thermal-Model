function [A,B,C,D,u] = buildSS(R,C,Qot,Qob,Qs,t)

heatZone = csvread('SnakeInput.csv');
heatZoneSize = size(heatZone);
heatZoneV = reshape(heatZone,numel(heatZone),1);
inputNodes = sum(heatZoneV);

n = heatZoneSize(1);
m = heatZoneSize(2);
N = n*m;

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
u = uQs;

% u = Qs * zeros(length(t),1); %step function; for all time, input = Qs
% Qsu = Qs * ones(10,1);
% u(1:10) = Qsu;

%% Output vectors for SS

C = speye(m*n);

D = 0;