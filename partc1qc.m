close all
clc
clear

a = [sqrt(3)/2 1/2; -1/2 sqrt(3)/2];
at = transpose(a);
sigmax = -45.385*10^6;
sigmay = 175.385*10^6;
sigmatensor = [sigmax 0; 0 sigmay];
sigmatransform = a*sigmatensor*at;

% Finding Principal Stress
sigma1 = [];
sigmaprinc1 = [];
sigmap1 = [-400*10^6:10^6:400*10^6];
for j = 1:length(sigmap1)
    sigmaprinc1(j) = (sigmatransform(1,1)+sigmap1(j)+sigmatransform(2,2))/2 + sqrt(((sigmatransform(1,1)+sigmap1(j)-sigmatransform(2,2))/2)^2 + (sigmatransform(1,2))^2);
end


sigma1 = [sigmap1; sigmaprinc1];
sigma1 = transpose(sigma1);

sigma2 = [];
sigmaprinc2 = [];
sigmap2 = [-400*10^6:10^6:400*10^6];

for j = 1:length(sigmap2)
    sigmaprinc2(j) = (sigmatransform(1,1)+sigmap2(j)+sigmatransform(2,2))/2 - sqrt(((sigmatransform(1,1)+sigmap2(j)-sigmatransform(2,2))/2)^2 + (sigmatransform(1,2))^2);
end

sigma2 = [sigmap2; sigmaprinc2];
sigma2 = transpose(sigma2);

% For Tresca
tresca = [];
tresc1 = [];
tresc2 = [];
tresc3 = [];
trescamax = [];

for i = 1:length(sigmap1)
    tresc1(i) = abs(sigma1(i,2) - sigma2(i,2));
    tresc2(i) = sigma2(i,2);
    tresc3(i) = sigma1(i,2);
    tresca = [tresc1(i) tresc2(i) tresc3(i)];
    trescamax(i) = max(tresca);
end
trescamax = transpose(trescamax);
trescamax = [transpose(sigmap1) trescamax];

% Finding valid sigmap
validsigmap = [];

for i = 1:length(sigmap1)
    if trescamax(i,2) >= -400*10^6 && trescamax(i,2) <= 400*10^6
        validsigmap(i) = trescamax(i,1);
    end
end

validsigmap = transpose(validsigmap);
validsigmap = validsigmap(161:758);
safetyfactor_tresc = (400*10^6)./trescamax(161:758,2);

% For Von Mises

vonmises = [];

for i = 1:length(sigmap1)
    diff1 = (sigma1(i,2) - sigma2(i,2))^2;
    diff2 = (sigma2(i,2))^2;
    diff3 = (sigma1(i,2))^2;
    vonmises(i) = (1/sqrt(2))*sqrt(diff1 + diff2 + diff3);
end

vonmises = [transpose(sigmap1) transpose(vonmises)];

% Finding valid sigmap
validsigmap_v = [];
for i = 1:length(sigmap1)
    if vonmises(i,2) >= -400*10^6 && vonmises(i,2) <= 400*10^6
        validsigmap_v(i) = vonmises(i,1);
    end
end

validsigmap_v = transpose(validsigmap_v);
validsigmap_v = validsigmap_v(103:800);

safetyfactor_vonmises = (400*10^6)./vonmises(103:800,2);

% For fracture 

