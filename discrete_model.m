clear; close all;

% 2D Thermal Model for a thin plate with zoned internal heat generation

%% Printing Options
simTitle = 'Mar 9';
% Prompts whether to print a new plot when the script is run
% Will plot and save file with some descriptive filename
% Fix .eps printing to tikz printing someday
prompt = 'Output plots to file? ["Yes" / Literally any other input]';
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
tempInitial = 0;

initialStates = tempInitial * ones(n*m,1);

% Properties for the material
materialHeatCapacity = 840; % J/(kg-K)
materialDensity = 1850; % kg/m^3
symRunTime = 6000; % run time: seconds
materialConductivity = 10; % W/(m-K)
plateLength = .16; % length: meter
plateThickness = .18; % thickness: meter
powerInput = 190; % Value of step power input (Watts)
plateArea = plateLength^2;
%%%powerLossSides = -15; % Value of step power loss Side
convectionCoeffTop = 15; % Convection coefficient at top
convectionCoeffBottom = 2; % Convection coefficient at bottom
powerLossTob = -convectionCoeffTop*(plateArea/N); % Value of step power loss Top
powerLossBottom = -convectionCoeffBottom*(plateArea/N); % Value of step power loss Bottom
%%%powerLossTob = 1; % Step value used for debug
%%%powerLossBottom = 1; % Step value used for debug


% Time Step and Vector
dt = symRunTime/100;
timeVector = 0:dt:symRunTime;

dx = plateLength/n; % length between each node
nodeCrossSection = dx*plateThickness; % cross sectional area of node

materialCapacitance = materialHeatCapacity * dx * nodeCrossSection * materialDensity; % thermal capacitance
materialResistance = dx / (materialConductivity * nodeCrossSection); % thermal resistance

%% Setup A Matrix
tic

%% Main diagonal construction: coefficient of i,j element
% Coefficient value for internal elements are inserted at every location, 
% then replaced if the element is on an edge.
tempA1 = (-(4*materialResistance^3)/(materialCapacitance*materialResistance^4)+powerLossTob+powerLossBottom)*ones(n*m,1); % Value of C_i,j for internal elements
tempA1(1:n) = -(3*materialResistance^2)/(materialCapacitance*materialResistance^3)+powerLossTob+powerLossBottom; % Value of C_i,j for left edge replacement
tempA1(n*m-(n-1):n*m) = -(3*materialResistance^2)/(materialCapacitance*materialResistance^3)+powerLossTob+powerLossBottom; % Value of C_i,j for right edge 
%    replacement

% This loop modifies the C_i,j value for edge elements along the top and
% bottom of the grid by checking if they are an even multiple of n or are 
% one larger than an even multiple of n.
for i = 1:n*m
    if mod(i,n) == 0
        tempA1(i) = -(3*materialResistance^2)/(materialCapacitance*materialResistance^3)+powerLossTob+powerLossBottom; % Bottom row
    elseif mod(i,n) == 1
        tempA1(i) = -(3*materialResistance^2)/(materialCapacitance*materialResistance^3)+powerLossTob+powerLossBottom; % Top row
    else
    end
end

% This block adjusts the corners to the correct coefficient values.
tempA1(1) = -(2*materialResistance)/(materialCapacitance*materialResistance^2)+powerLossTob+powerLossBottom;
tempA1(n) = tempA1(1);
tempA1(n*m) = tempA1(1);
tempA1(n*m-(n-1)) = tempA1(1);

%% Auxiliary diagonal construction: coefficients of adjacent elements

% tempA2 will map to two diagonals directly next to the main diagonal and
% represents resistance between C_i,j and C_i,j(+/-)1 (elements east and
% west). The A Matrix is constructed by first going down the left column
% and finishing by going down the right column, so offsetting the diagonal
% by n rows takes care of removing resistances where there is no element to
% the east or west.
tempA2 = (1/(materialCapacitance*materialResistance))*ones(n*m-n,1);

% tempA3 will map to a diagonal which represents resistance between C_i,j and
% C_i-1,j (element south)
tempA3 = (1/(materialCapacitance*materialResistance))*ones(n*m-1,1);

% This loop removes elements from tempA3 at the bottom edge of the grid array,
% where there is no thermal element to the south.
for i = 1:n*m-1
    if mod(i,n) == 0
        tempA3(i) = 0;
    else
    end
end

% tempA4 will map to a diagonal which represents resistance between C_i,j and
% C_i+1,j (element north)
tempA4 = (1/(materialCapacitance*materialResistance))*ones(n*m-1,1);

% This loop removes elements from tempA4 at the top edge of the grid array,
% where there is no thermal element to the north.
for i = 2:n*m-1
    if mod(i,n) == 1
        tempA4(i-1) = 0;
    else
    end
end

%% Combine state variable coefficient values into sparse A matrix.
v = zeros(n,1); % zeroes required to offset east and west edges
tempAEast = [v;tempA2]; % elements to the east
tempAWest = [tempA2;v]; % elements to the west
v = [0]; % zero required to offset north and south edges
tempA3 = [tempA3;v]; % elements to the south
tempA4 = [v;tempA4]; % elements to the north
% tempA1(1) = -1/(materialCapacitance*materialResistance);
% tempA1(n*m) = -1/(materialCapacitance*materialResistance);

% Vectors required to create sparse A matrix
A = [tempAWest tempA3 tempA1 tempA4 tempAEast];
v = [-n -1 0 1 n]';

A = spdiags(A,v,N,N);

A_full = full(A); % For visual debugging of A matrix

timeDiag = toc;

% Input Vectors (sysBComponents = Heat in; sysBOutTop,sysBOutBottom,sysBOutSide = Heat out (top,bottom,side))
sysBComponents = zeros(n*m,1);
%%%sysBOutTop = zeros(n*m,1); Legacy code; Replaced with convection instead of step power loss.
%%%sysBOutBottom = ones(n*m,1); Legacy code; Replaced with convection instead of step power loss.
%%%sysBOutSide = zeros(n*m,1); Legacy code; Replaced with convection instead of step power loss.

%% Control for where heat is input by changing sysB to weight specific nodes. Comment/Uncomment to select. Currently very clumsy. Better to use .csv file.

sysBComponents = 1/(materialCapacitance*materialResistance)*heatZoneV; % Loads sysB with location information from .csv file.

% sysBComponents = ones(n*m,1); sysBComponents = (1/(materialCapacitance*materialResistance))*sysBComponents;
% sysBComponents(ceil(m*n/2)) = 1/(materialCapacitance*materialResistance); % Center Element
% sysBComponents(1) = 1/(materialCapacitance*materialResistance); % First Element
% sysBComponents(m*n-(n-1):m*n) = 1/(materialCapacitance*materialResistance); % All Right Edge
% sysBComponents(1:n) = 1/(materialCapacitance*materialResistance); % All Left Edge
% sysBComponents(m*n-(floor(n/2))) = 1/(materialCapacitance*materialResistance);
% sysBComponents(n+2) = 1/(materialCapacitance*materialResistance);
% sysBComponents(N) = 1/(materialCapacitance*materialResistance); % Last Element
% sysBComponents(29:70) = 1/(materialCapacitance*materialResistance); % Some Random Elements
% sysBComponents(587:615) = 1/(materialCapacitance*materialResistance); % Some Random Elements

% sysBOutSide(1:n) = 1/(materialCapacitance*materialResistance); % All Left Edge
% sysBOutSide(m*n-(n-1):m*n) = 1/(materialResistance*materialCapacitance); % All Right Edge
% sysBOutSide(n*m) = 1/(materialCapacitance*materialResistance); % Last Element
% sysBOutSide(1) = 1/(materialCapacitance*materialResistance); % First Element

%% Input vectors for SS

% sysBOutTop(m*n-(n-1):m*n) = 1/(materialCapacitance*materialResistance);

%%%sysBOutTop = (1/(materialCapacitance*materialResistance))*sysBOutTop;
%%%sysBOutBottom = (1/(materialCapacitance*materialResistance))*sysBOutBottom;
%%%sysBOutSide = (1/(materialCapacitance*materialResistance))*sysBOutSide;
%%%sysB = [sysBComponents sysBOutTop sysBOutBottom sysBOutSide];

sysB = [sysBComponents];

% sysB = sysB * 1/(materialCapacitance*materialResistance);

%% Power Input
% Inputs vectors (power input and power loss)
sysUPowerInput = zeros(length(timeVector),1);
%%%uQot = zeros(length(timeVector),1); Legacy code; Replaced with convection instead of step power loss.
%%%uQob = zeros(length(timeVector),1); Legacy code; Replaced with convection instead of step power loss.
%%%uQos = zeros(length(timeVector),1); Legacy code; Replaced with convection instead of step power loss.

% Power Input Vector (used to set time power is on or off)
sysUPowerInput(1:length(timeVector),1) = powerInput/(nnz(sysBComponents));
%%%uQos(1:length(timeVector)/3,1) = powerLossSides/(nnz(sysBComponents));

% sysUPowerInput((length(timeVector)-length(timeVector)/3):length(timeVector),1) = powerInput/(nnz(sysBComponents));
%uQos((length(timeVector)-length(timeVector)/3):length(timeVector),1) = powerLossSides/(nnz(sysBComponents));
%%%uQot(1:length(timeVector),1) = powerLossTob/N;
%%%uQob(1:length(timeVector),1) = powerLossBottom/N;

% Combine inputs for lsim
%%%sysU = [sysUPowerInput uQos uQot uQob];
sysU = [sysUPowerInput];

% sysU = powerInput * zeros(length(timeVector),1); %step function; for all time, input = powerInput
% Qsu = powerInput * ones(10,1);
% sysU(1:10) = Qsu;

%% Output vectors for SS

sysC = speye(m*n);

sysD = 0;

%% Simulation
tic
sys = ss(A,sysB,sysC,sysD);
sysOutput = lsim(sys,sysU,timeVector,initialStates);
timeSim = toc;

%% Reorder results into element grid and plot/print contours
timecode = num2str(now);

for i = 0:n-1
    for j = 2:length(timeVector)-1
        M(:,i+1,j) = sysOutput(j,i*n+1:i*n+n);
    end
end

tempFinal = sysOutput(end,:);
maxTemp = max(tempFinal);
minTemp = min(tempFinal);

tenCont = floor((maxTemp-minTemp)/10);

M1 = M(:,:,2);

plots_script

beep