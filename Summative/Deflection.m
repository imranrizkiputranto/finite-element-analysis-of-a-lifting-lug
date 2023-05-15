clear
clc
close all

% Calculating the analytical deflection for a cantilever beam and comparing
% mesh studies
% 
% Imran Rizki Putranto
% 30 November 2022

% Parameters

stress = 249.8; % Max stress at contact area

% Linear Hex Deflection
linearhexdata = xlsread('ElementTypeVariation.xlsx');
elements = linearhexdata(:,2);
linearmaxstress = linearhexdata(:,3);
quadmaxstress = linearhexdata(1:35,10);

hold on
plot(elements,linearmaxstress)
plot(elements,quadmaxstress)
yline(stress)
hold off
 
title('Mesh Convergence Study for Linear and Quadratic Element Orders')
xlabel('Number of Elements')
ylabel('Stress (MPa)')
legend('Linear Mesh', 'Quadratic Mesh', 'Analytical Stress')