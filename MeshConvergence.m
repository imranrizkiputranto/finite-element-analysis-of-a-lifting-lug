clc
clear
close all

% Plotting Mesh Convergence Graphs
%

%% Quadratic Tet Mesh, C3D10H
data = xlsread('ElementTypeVariation.xlsx');

SeedSizes(:,1) = data(117:135,8);
SeedSizesInverse = 1./SeedSizes;
StressOutside(:,1) = data(117:135,12);
Deflection(:,1) = data(117:135,13);
n = length(SeedSizesInverse);

% Calculating Residuals
% y1(:,1) = [155.776 157.567 154.556 152.091];
% difference(:,1) = (StressOutside - y1).^2;
% RMS = sqrt(sum(difference)/(n-2));

%% Quadratic Wedge Mesh, C3D15H

WedgeSeedSizes(:,1) = data(75:93,8);
WedgeSeedSizesInverse = 1./WedgeSeedSizes;
WedgeStressOutside(:,1) = data(75:93,12);
WedgeDeflection(:,1) = data(75:93,13);
n = length(WedgeSeedSizesInverse);

%% Quadratic Hex Mesh, C3D10I

HexSeedSizes(:,1) = data(44:62,8);
HexSeedSizesInverse = 1./HexSeedSizes;
HexStressOutside(:,1) = data(44:62,12);
HexDeflection(:,1) = data(44:62,13);
n = length(HexSeedSizesInverse);

%% Plotting Stress

% Plot Stress Against 1/Element
Tet = fit(SeedSizesInverse,StressOutside,'poly2');
plot(SeedSizesInverse,StressOutside,'k-o')
hold on
title('Maximum Stress Against 1/Element Size')
plot(WedgeSeedSizesInverse,WedgeStressOutside, 'r-o')
plot(HexSeedSizesInverse,HexStressOutside,'b-o')
xlabel('1/Element Size')
ylabel('Maximum Stress Outside Contact Region [MPa]')
%scatter(SeedSizesInverse,StressOutside,'MarkerFaceColor',[0 .75 .75])
%scatter(WedgeSeedSizesInverse,WedgeStressOutside,'MarkerFaceColor',[0 .25 .25])
%scatter(HexSeedSizesInverse,HexStressOutside,'MarkerFaceColor',[0 .5 .5])
legend('Tet Mesh','Wedge Mesh','Hex Mesh')% 'Tet Mesh Points','Wedge Mesh Points','Hex Mesh Points')

%% Plotting Deflection

figure 
plot(SeedSizesInverse,Deflection,'k-o')
hold on
title('Maximum Deflection Against 1/Element Size')
plot(WedgeSeedSizesInverse,WedgeDeflection, 'r-o')
plot(HexSeedSizesInverse,HexDeflection,'b-o')
xlabel('1/Element Size')
ylabel('Maximum Deflection [mm]')
legend('Tet Mesh','Wedge Mesh','Hex Mesh')

%% Plotting Best Fit

function [fitresult, gof] = createFit1(x, y)
%% Fit: 'untitled fit 1'.

 
x = SeedSizesInverse;

y = StressOutside(:,1);

 

[xData, yData] = prepareCurveData(x,y);


% Set up fittype and options.

ft = fittype( 'power2' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

opts.Display = 'Off';

opts.StartPoint = [284.601596584229 0.206279953048744 -0.806336149906312];

 

% Fit model to data.

[fitresult, gof] = fit( xData, yData, ft, opts );

 

% Plot fit with data.

figure( 'Name', 'untitled fit 1' );

h = plot( fitresult, xData, yData );

%legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );

% Label axes

xlabel( '1/Element Size', 'Interpreter', 'none' );

ylabel( 'Maximum Stress [MPa]', 'Interpreter', 'none' );

grid on

 
%% 
cftool
end
