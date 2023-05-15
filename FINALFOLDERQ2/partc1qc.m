close all
clc
clear

%% test

E = 70*10^9;
v = 0.3;
eps_xx = fn_strain(2055792,0);
eps_yy = fn_strain(2055792,90);

a = [sqrt(3)/2 1/2; -1/2 sqrt(3)/2]; % Transformation matrix at 30 degrees
at = transpose(a); % Transpose of a

sigmax = (E/(1-v^2))*(eps_xx + 0.3*eps_yy);
sigmay = (E/(1-v^2))*(eps_yy + 0.3*eps_xx);

sigmatensor = [sigmax 0; 0 sigmay]; % Stress tensor 
sigmatransform = a*sigmatensor*at; % Transformed tensor at 30 degrees

% Finding Principal Stress
sigmaprinc1 = [];
sigmap1 = [-400*10^6:10^6:400*10^6]; % Range of sigma_p values to try 

% Finding principal stress 1 for -400MPa to 400MPa
for j = 1:length(sigmap1)
    sigmaprinc1(j) = (sigmatransform(1,1)+sigmap1(j)+sigmatransform(2,2))/2 + sqrt(((sigmatransform(1,1)+sigmap1(j)-sigmatransform(2,2))/2)^2 + (sigmatransform(1,2))^2);
end

% Making array: [ sigma_p principalstress_1]
sigma1 = [];
sigma1 = [sigmap1; sigmaprinc1];
sigma1 = transpose(sigma1);

% Finding principal stress 2 for -400MPa to 400MPa
sigmaprinc2 = [];
sigmap2 = [-400*10^6:10^6:400*10^6];

for j = 1:length(sigmap2)
    sigmaprinc2(j) = (sigmatransform(1,1)+sigmap2(j)+sigmatransform(2,2))/2 - sqrt(((sigmatransform(1,1)+sigmap2(j)-sigmatransform(2,2))/2)^2 + (sigmatransform(1,2))^2);
end

% Making array: [ sigma_p principalstress_2]
sigma2 = [];
sigma2 = [sigmap2; sigmaprinc2];
sigma2 = transpose(sigma2);

figure
plot(sigma1(:,1)./10^6,sigma1(:,2)./10^6)
hold on
plot(sigma2(:,1)./10^6,sigma2(:,2)./10^6)
title('Plot Showing How Principal Stresses 1 and 2 vary with Applied Stress P')
xlabel('Applied Stress P [MPa]')
ylabel('Principal Stress Values')
legend('Principal Stress 1', 'Principal Stress 2')

% For Tresca criterion
trescamax = [];
for i = 1:length(sigmap1)
    tresc1 = abs(sigma1(i,2) - sigma2(i,2));
    tresc2 = sigma2(i,2);
    tresc3 = sigma1(i,2);
    tresca = [tresc1 tresc2 tresc3];
    trescamax(i) = max(tresca);
end

trescamax = transpose(trescamax);
trescamax = [transpose(sigmap1) trescamax];

% Finding valid sigma_p for tresca criterion giving stress below yield
validsigmap = [];
for i = 1:length(sigmap1)
    if trescamax(i,2) <= 400*10^6
        validsigmap(i) = trescamax(i,1);
    end
end

validsigmap = transpose(validsigmap);
% validsigmap = validsigmap(161:758);
safetyfactor_tresc = (400*10^6)./trescamax(:,2);

% For Von Mises

vonmises = [];

for i = 1:length(sigmap1)
    diff1 = (sigmaprinc1(i) - sigmaprinc2(i))^2;
    diff2 = (sigmaprinc2(i))^2;
    diff3 = (-sigmaprinc1(i))^2;
    vonmises(i) = (1/sqrt(2))*sqrt(diff1 + diff2 + diff3);
end

vonmises = [transpose(sigmap1) transpose(vonmises)];

% Finding valid sigmap
validsigmap_v = [];
for i = 1:length(sigmap1)
    if vonmises(i,2) <= 400*10^6
        validsigmap_v(i) = vonmises(i,1);
    end
end

validsigmap_v = transpose(validsigmap_v);
% validsigmap_v = validsigmap_v(103:800);

safetyfactor_vonmises = (400*10^6)./vonmises(:,2);

figure
plot(sigmap1./10^6,safetyfactor_tresc)
hold on
plot(sigmap1./10^6,safetyfactor_vonmises)
yline(1)
title('Plot of Safety Factors for Both Yield Criteria Against Applied Stress')
xlabel('Applied Stress P [MPa]')
ylabel('Safety Factor')
legend('Tresca','Von Mises','Valid Safety Factor Limit')

figure
hold on
plot(trescamax(:,1)/10^6,trescamax(:,2)/10^6)
plot(vonmises(:,1)/10^6,vonmises(:,2)/10^6)
title('Plot of Von Mises and Tresca Stress Criterion for a Range of Sigma P')
xlabel('Applied Stress P [MPa]')
ylabel('Stress [MPa]')
yline(400)
legend('Tresca Criterion','Von Mises Criterion','Yield Strength of Material')

%% Brittle Fracture
% Assuming design stress = normal to crack, datum at 30 degrees rotation

transformedsigmayy = [];
cl = 0.0025;
i = 1;
a = {};
finalsigmatensor = {};

for sigmap = -400*10^6:10^6:400*10^6
    finalsigmatensor{i} = [sigmatransform(1,1)+sigmap, sigmatransform(1,2); sigmatransform(2,1) sigmatransform(2,2)];
    for thet = 1:1:361
        a{thet} = [cosd(thet-1) sind(thet-1); -sind(thet-1) cosd(thet-1)];
        rotatedstresstensor{i,thet} = a{thet}*finalsigmatensor{i}*transpose(a{thet});
    end
    i = i+1;
end

% LEFM
sigmay = [];
for i = 1:length([-400*10^6:10^6:400*10^6])
    for j = 1:length(0:pi/180:2*pi)
        Kc(j,i) = ((rotatedstresstensor{i,j}(2,2))*sqrt(pi*cl));
    end
end

sigma_p = [-400*10^6:10^6:400*10^6];
angle = [30:1:390];
[X,Y] = meshgrid(sigma_p,angle);

figure
s = contourf(X./10^6,Y,Kc/10^6,'LevelList',20); 
hold on
%colorbar
title('Plot of LEFM Stress Intensity Factor for Every Value of Applied Stress and Crack Orientation')
xlabel('Applied Stress P [MPa]')
ylabel('Angle [Degrees]')
% ylabel(colorbar, 'Stress Intensity Factor [MPa√m]', 'Fontsize',11);

% Contour plot for LEFM Validity
stress_intensity_lefm = [];
for i = 1:length([-400*10^6:10^6:400*10^6])
    for j = 1:length(0:pi/180:2*pi)
        LEFM = (4/pi)*(Kc(j,i)/(400*10^6))^2;
        stress_intensity_lefm(j,i) = LEFM;
    end
end

figure 
h = contourf(X,Y,stress_intensity_lefm,'LevelList',2.5*10^-3);
colorbar
title('Plot of LEFM Validity Equation Value for Every Value of Applied Stress and Crack Orientation')
xlabel('Applied Stress P [MPa]')
ylabel('Angle [Degrees]')
% ylabel(colorbar, 'LEFM Validity Equation Value', 'Fontsize',11);

% Plot for K against LEFM Validity
max_stress_intensity = [];
for i = 1:length(sigma_p)
    max_stress_intensity(i,1) = max(Kc(:,i));
    max_stress_intensity(i,2) = (4/pi)*(max_stress_intensity(i,1)/(400*10^6))^2;
end

figure
plot(max_stress_intensity(:,1),max_stress_intensity(:,2)*1000)
yline(2.5)
title('Maximum LEFM Validity Value from All Possible Crack Orientations for Each Value of Applied Stress')
xlabel('Stress Intensity Factor, K')
ylabel('LEFM Validity Equation')
legend('LEFM Validity Equation', 'Half crack length, a = 2.5mm')

%% EPFM
sigmay = [];
for i = 1:length([-400*10^6:10^6:400*10^6])
    for j = 1:length(0:pi/180:2*pi)
%         if (4/pi)*(Kc(j,i)/(400*10^6))^2 > 0.0025
            Kc(j,i) = ((rotatedstresstensor{i,j}(2,2))*sqrt(pi*cl))/sqrt(1-((1/2)*((rotatedstresstensor{i,j}(2,2))/(400*10^6))^2));
%         else
%             continue
%         end
    end
end

sigma_p = [-400*10^6:10^6:400*10^6];
angle = [30:1:390];
[X,Y] = meshgrid(sigma_p,angle);

% Plotting EPFM Contour
figure

% Tresca
pgon = polyshape([-297 -297 -240 -240],[390 30 30 390]);
plot(pgon)
hold on

pgon2 = polyshape([400 400 362 362],[390 30 30 390]);
plot(pgon2)

% Vonmises
pgon3 = polyshape([-400 -400 -297 -297],[390 30 30 390]);
plot(pgon3)

s = contourf(X./10^6,Y,Kc/10^6,'LevelList',20);
title('Plot of EPFM Stress Intensity Factor for Every Value of Applied Stress and Crack Orientation')
xlabel('Applied Stress P [MPa]')
ylabel('Angle [Degrees]')
% ylabel(colorbar, 'Stress Intensity Factor [MPa√m]', 'Fontsize',11);
legend('Tresca Lower Limit', 'Tresca Upper Limit', 'Von Mises Lower Limit')

%% For Curved Plot
finalsigmatensor = [];
i = 1;
for sigmap = -400*10^6:10^6:400*10^6
    for thet = 0:pi/180:2*pi
        finalsigmatensor = [sigmatransform(1,1)+sigmap, sigmatransform(1,2); sigmatransform(2,1) sigmatransform(2,2)];
        a = [cos(thet) sin(thet); -sin(thet) cos(thet)];
        transformedtensor = a*finalsigmatensor*transpose(a);
        transformedsigmayy(i,1) = sigmap;
        transformedsigmayy(i,2) = thet;
        transformedsigmayy(i,3) = transformedtensor(2,2);
        i = i+1;
    end
end

finalfracturematrix =[];
j=1;

for cl = 0.0025
    for i = 1:length(transformedsigmayy)
        if (4/pi)*((transformedsigmayy(i,3)*sqrt(pi*cl))/(400*10^6))^2 > 0.0025
            finalfracturematrix(j,1) = transformedsigmayy(i,1); % sigmap
            finalfracturematrix(j,2) = transformedsigmayy(i,2); % theta
            finalfracturematrix(j,3) = cl; % length a (half crack length)
            finalfracturematrix(j,4) = (transformedsigmayy(i,3)*sqrt(pi*cl))/sqrt(1-((1/2)*(transformedsigmayy(i,3)/(400*10^6))^2)); % K if < 20
            j = j+1;
        elseif (4/pi)*((transformedsigmayy(i,3)*sqrt(pi*cl))/(400*10^6))^2 <= 0.0025
            finalfracturematrix(j,1) = transformedsigmayy(i,1); % sigmap
            finalfracturematrix(j,2) = transformedsigmayy(i,2); % theta
            finalfracturematrix(j,3) = cl; % length a (half crack length)
            finalfracturematrix(j,4) = (transformedsigmayy(i,3)*sqrt(pi*cl)); % K if < 20
            j = j+1;
        end
    end
end

% && transformedsigmayy(i,3)*sqrt(pi*cl) <= 20*10^6

[M,I] = max(finalfracturematrix(:,4));

maxvaluescell = {};
plotmatrix = [];
sigmapp = -400*10^6;
sigmappp = [-400*10^6:10^6:400*10^6];
sigmapp = transpose(sigmapp);
j = 1;
k = 1;

while sigmapp <= 400*10^6
    for i = 1:length(finalfracturematrix)
        if finalfracturematrix(i,1) == sigmapp
           maxvaluescell{k}(j,1) = sigmapp;
           maxvaluescell{k}(j,2) = finalfracturematrix(i,2);
           maxvaluescell{k}(j,3) = finalfracturematrix(i,4);
        else
            continue
        end
        j = j+1;
    end
    j = 1;
    sigmapp = sigmapp + 10^6;
    k = k+1;
end

maxvaluescell = transpose(maxvaluescell);

plotmatrix(:,1) = sigmappp;
k = 1;
for i = 1:length(maxvaluescell)
    [plotmatrix(k,3), IDX] = max(maxvaluescell{i}(:,3)); % [max stress, index]
    plotmatrix(k,2) = maxvaluescell{i}(IDX,2);
    k = k+1;
end

%% Plotting for K Limits
figure
plot(plotmatrix(:,1)./10^6,plotmatrix(:,3)./10^6)
hold on
yline(20)
title('Plot of Maximum Stress Intensity Factors for Each Value of Applied Stress P')
xlabel('Applied Stress P [MPa]')
ylabel('Stress Intensity Factor K [MPa√m]')
% figure
% scatter3(plotmatrix(:,1),plotmatrix(:,2),plotmatrix(:,3),'MarkerFaceColor',[0 .75 .75])


%% Polar plot - Determine from matrices ^ run above section first
% cl = 0.0025m

% Upper limit = 106 MPa sigma p for fracture
trescaupper = [];
for i = 1:length(maxvaluescell)
    if maxvaluescell{i}(1,1) == 106*10^6
        trescaupper(:,1) = maxvaluescell{i}(:,2) + pi/6;
        trescaupper(:,2) = maxvaluescell{i}(:,3);
    end
end

figure
polarplot(trescaupper(:,1), trescaupper(:,2))
title('Stress Intensity Factor K for Different Crack Orientations for Upper Limit')    

% % Lower Limit = 26 MPa sigma p for fracture
% trescalower = [];
% for i = 1:length(maxvaluescell)
%     if maxvaluescell{i}(1,1) == 26*10^6
%         trescalower(:,1) = maxvaluescell{i}(:,2) + pi/6;
%         trescalower(:,2) = maxvaluescell{i}(:,3);
%     end
% end
% 
% figure
% polarplot(trescalower(:,1), trescalower(:,2))
% title('Stress Intensity Factor K for Different Crack Orientations for Lower Limit') 

