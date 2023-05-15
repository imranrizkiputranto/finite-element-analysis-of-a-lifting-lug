% Calculating Hertzian Contact Stress

% Calculation parameters

F = 58860; % Force in Newtons
l = 50; % Length of pin in mm
mu = 0.3; % Poisson's Ratio
E = 209000; % Young's Modulus in MPa
d = 20; % Pin and Lug Diameter in mm

% Contact half width
b1 = ((2*F)/(pi*l));
b2 = (((1-mu^2)/E)+((1-mu^2)/E))/((1/d)+(1/d));
b = sqrt(b1*b2) 

% Max pressure
P = 2*F/(pi*b*l)

% Principal Stresses
% StressX = -2*mu*P*(sqrt