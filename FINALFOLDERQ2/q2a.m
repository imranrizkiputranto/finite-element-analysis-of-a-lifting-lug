clc
close all
clear

%% Plotting a against N
%
%

a = [0.3:0.005:10];
a_i = a(1);
a = transpose(a);
t  = 10;
sigma_max = 200*10^6;

% Probability
pod = [];
for i = 1:length(a)
    pod(i,1) = a(i);
    pod(i,2) = fn_pod(2055792,a(i));
end

% plot of POD against crack length
plot(a*2,pod(:,2))
title('Plot of Probability of detection against crack length')
xlabel('Crack length, 2a [mm]')
ylabel('Probability of detection')

%% Finding number of cycles until failure

% Paris law integral

% N = 1/(C*(Y*sigma_max/10^6*sqrt(pi/a)))^m;
C = 10^-12;
a = 0.0003; % Initial half crack length
a_f = 0.01; % 1mm initial half crack length
t = 0.01; % Thickness of pipe
FatigueLife = [];
i = 1;
m = 4;

while a <= 0.01
    FatigueLife(i,1) = a*1000; % Initial half crack length in mm
    fun = @(a) 1./(C.*((0.728 + 0.373.*(a/t).^2 - 0.029.*(a/t).^4)*(sigma_max./10^6).*sqrt(pi.*a)).^m); % Number of cycle
    FatigueLife(i,2) = integral(fun,a,0.01); % Fatigue life at specific initial crack length
    a = a + 0.00001; % 0.01mm iteration
    i = i+1;
end

figure
plot(FatigueLife(:,1)*2,FatigueLife(:,2))
title('Plot of number of cycles until failure against initial crack length')
xlabel('Initial crack length, 2a [mm]')
ylabel('Number of cycles until failure')

%% Paris Law
parislaw = [];
sigma_max = 200*10^6;
c_length = 0.0003;
i = 1;
t = 0.01;
C = 10^-12;
m = 4;
while c_length <= 0.01
    Y = 0.728 + 0.373*(c_length/t)^2 - 0.029*(c_length/t)^4;
    K = Y*(sigma_max/10^6)*sqrt(pi*c_length);
    parislaw(i,1) = c_length; % initial crack length (mm)
    parislaw(i,2) = K; % stress intensity factor, MPa(m)^1/2
    parislaw(i,3) = C*(K)^m; % length increase per cycle (m)
    c_length = c_length + parislaw(i,3); % length after increase during cycle (mm)
    parislaw(i,4) = c_length*10^3;
    parislaw(i,5) = i;
    i = i+1;
end

%% Plot
figure
title('Plot of crack length after each cycle against number of cycles that have occured')
plot(parislaw(:,5),parislaw(:,4)*2)
xlabel('Number of cycles that have occured')
ylabel('Crack length after each cycle, 2a [mm]')

figure
title('Plot of Crack Growth Rate against Stress Intensity Range')
plot(parislaw(:,2),parislaw(:,3)*10^3)
xlabel('Stress Intensity Factor [MPa√m]')
ylabel('Crack Growth Rate in Log Scale [mm/cycle]')
set(gca, 'YScale', 'log')
%% Plotting prob of undetection against inspection interval

fail = [];
o = 1;
for cycles = 8000:1:12000
    insp_int = [];
    p_undet = [];
    for i = 1:length(parislaw)
        p_undet(i,1) = 1 - fn_pod(2055792,parislaw(i,1)*1000); % Probability of undetected crack
    end
    
    init_p_undet = [];
    for i = 1:length(parislaw)
        init_p_undet(i,1) = 1 - fn_pod(2055792,parislaw(i,1)*1000); % Probability of undetected crack
    end
    
    inspection_iterations = [];
    for i = 1:length(parislaw)
        insp_int(i,1) = parislaw(i,5)/cycles; % Inspections per 1000 cycles;
        if rem(insp_int(i,1),1) == 0
            inspection_iterations = insp_int(i,1);
        end
    end
    
    m = 0;
    
    for i = 1:length(p_undet)
        if (i/cycles) == floor(i/cycles) % if an inspection takes place
            m = m+1; % Number of inspections increase by 1
        end
    
        n = i/inspection_iterations;
        if n > 1 && m > 1
            if (i/cycles) == floor(i/cycles)
                p_undet(i,1) = p_undet(i,1) * p_undet(cycles*(m-1),1);
            elseif (i/cycles) ~= floor(i/cycles)
                p_undet(i,1) = p_undet(i,1) * p_undet(cycles*(m),1);
            end
        elseif n < 1 && m < 1
            p_undet(i,1) = p_undet(i,1);
        elseif n >= 1 && m == 1
            if (i/cycles) == floor(i/cycles)
                p_undet(i,1) = p_undet(i,1);
            elseif (i/cycles) ~= floor(i/cycles)
                p_undet(i,1) = p_undet(i,1) * p_undet(cycles*(m),1);
            end
        end
    end
    
    test = [];
    j = 1;
    for i = 1:length(p_undet)
        if i/cycles == floor(i/cycles)
            test(j,1) = i/cycles;
            test(j,2) = p_undet(i,1);
            j = j+1;
        else
            continue
        end
    end
    
    fail(o,2) = test(end,2); % Prob failure
    fail(o,1) = cycles;
    o = o+1;
end

%%
figure
scatter(fail(:,1),fail(:,2),'MarkerFaceColor',[0 .75 .75])
yline(0.01)
title('Probability of Failure Against Inspection Interval')
xlabel('Inspection Interval')
ylabel('Probability of Failure')

%% Matrix for scatter plot - Prob = 0.01
% scatter_undet = [];
% j = 1;
% for i = 1:length(p_undet)
%     if (i/cycles) == floor(i/cycles)
%         scatter_undet(i,1) = p_undet(i,1);
%     else
%         scatter_undet(i,1) = nan;
%     end
% end
%     
% fail = test(end,2);
% 
% % iterations = [0:1:inspection_iterations];
% fprintf('Probability of failure is %s.\n',fail);
% figure
% scatter(parislaw(:,5), scatter_undet, 'MarkerFaceColor',[0 .75 .75])
% %[scatterpoly,gof] = fit(insp_int,p_undet,'poly2');
% % plot(scatterpoly,insp_int,scatter_undet);
% hold on
% plot(parislaw(:,5),init_p_undet)
% legend('1','2');

%% Plotting Prob undetected after each inspection against number of
% inspections for normal cycle - no surge

cycles = 9682;
insp_int = [];
p_undet = [];
for i = 1:length(parislaw)
    p_undet(i,1) = 1 - fn_pod(2055792,parislaw(i,1)*1000); % Probability of undetected crack
end

init_p_undet = [];
for i = 1:length(parislaw)
    init_p_undet(i,1) = 1 - fn_pod(2055792,parislaw(i,1)*1000); % Probability of undetected crack
end

inspection_iterations = [];
for i = 1:length(parislaw)
    insp_int(i,1) = parislaw(i,5)/cycles; % Inspections per 1000 cycles;
    if rem(insp_int(i,1),1) == 0
        inspection_iterations = insp_int(i,1);
    end
end

m = 0;

for i = 1:length(p_undet)
    if (i/cycles) == floor(i/cycles) % if an inspection takes place
        m = m+1; % Number of inspections increase by 1
    end

    n = i/inspection_iterations;
    if n > 1 && m > 1
        if (i/cycles) == floor(i/cycles)
            p_undet(i,1) = p_undet(i,1) * p_undet(cycles*(m-1),1);
        elseif (i/cycles) ~= floor(i/cycles)
            p_undet(i,1) = p_undet(i,1) * p_undet(cycles*(m),1);
        end
    elseif n < 1 && m < 1
        p_undet(i,1) = p_undet(i,1);
    elseif n >= 1 && m == 1
        if (i/cycles) == floor(i/cycles)
            p_undet(i,1) = p_undet(i,1);
        elseif (i/cycles) ~= floor(i/cycles)
            p_undet(i,1) = p_undet(i,1) * p_undet(cycles*(m),1);
        end
    end
end

test = [];
j = 1;
for i = 1:length(p_undet)
    if i/cycles == floor(i/cycles)
        test(j,1) = i/cycles;
        test(j,2) = p_undet(i,1);
        j = j+1;
    else
        continue
    end
end


figure
scatter(test(:,1),test(:,2), 'MarkerFaceColor',[0 .75 .75])
title('Probability of Not Detecting a Crack At Each Inspection')
xlabel('Number of Inspections')
ylabel('Probability of Not Detecting a Crack At Each Inspection')

%% Plotting Prob Undetected After Each Inspection For Surge

% p_surge = 80; % MPa
% p_axialsurge = (80*0.1)/(2*0.01); % MPa
% 
% % Create Single Surge
% SurgeOccurence = 173990;
% 
% % Surge Paris Law
% ParisSurge= [];
% cLengthSurge = 0.001;
% i = 1;
% t = 0.01;
% C = 10^-12;
% m = 4;
% 
% while cLengthSurge <= 0.01
%     if i == SurgeOccurence
%         Y = 0.728 + 0.373*(cLengthSurge/t)^2 - 0.029*(cLengthSurge/t)^4;
%         K = Y*(p_axialsurge)*sqrt(pi*cLengthSurge);
%         ParisSurge(i,1) = cLengthSurge; % initial crack length (mm)
%         ParisSurge(i,2) = K; % stress intensity factor, MPa(m)^1/2
%         ParisSurge(i,3) = C*(K)^m; % length increase per cycle (m)
%         cLengthSurge = cLengthSurge + ParisSurge(i,3); % length after increase during cycle (mm)
%         ParisSurge(i,4) = cLengthSurge*10^3;
%         ParisSurge(i,5) = i;
%         i = i+1;
%     else
%         Y = 0.728 + 0.373*(cLengthSurge/t)^2 - 0.029*(cLengthSurge/t)^4;
%         K = Y*(sigma_max/10^6)*sqrt(pi*cLengthSurge);
%         ParisSurge(i,1) = cLengthSurge; % initial crack length (mm)
%         ParisSurge(i,2) = K; % stress intensity factor, MPa(m)^1/2
%         ParisSurge(i,3) = C*(K)^m; % length increase per cycle (m)
%         cLengthSurge = cLengthSurge + ParisSurge(i,3); % length after increase during cycle (mm)
%         ParisSurge(i,4) = cLengthSurge*10^3;
%         ParisSurge(i,5) = i;
%         i = i+1;
%     end
% end
% 
% SurgeUndet = [];
% for i = 1:length(ParisSurge)
%     SurgeUndet(i,1) = 1 - fn_pod(2055792,ParisSurge(i,1)*1000); % Probability of undetected crack
% end
% 
% SurgeInspInt = [];
% SurgeInspectionIterations = [];
% for i = 1:length(ParisSurge)
%     SurgeInspInt(i,1) = ParisSurge(i,5)/cycles; % Inspections per 1000 cycles;
%     if rem(SurgeInspInt(i,1),1) == 0
%         SurgeInspectionIterations = SurgeInspInt(i,1);
%     end
% end
% 
% m = 0;
% for i = 1:length(SurgeUndet)
%     if (i/cycles) == floor(i/cycles) % if an inspection takes place
%         m = m+1; % Number of inspections increase by 1
%     end
% 
%     n = i/SurgeInspectionIterations;
%     if n > 1 && m > 1
%         if (i/cycles) == floor(i/cycles)
%             SurgeUndet(i,1) = SurgeUndet(i,1) * SurgeUndet(cycles*(m-1),1);
%         elseif (i/cycles) ~= floor(i/cycles)
%             SurgeUndet(i,1) = SurgeUndet(i,1) * SurgeUndet(cycles*(m),1);
%         end
% 
%     elseif n < 1 && m < 1
%         SurgeUndet(i,1) = SurgeUndet(i,1);
%     elseif n >= 1 && m == 1
%         if (i/cycles) == floor(i/cycles)
%             SurgeUndet(i,1) = SurgeUndet(i,1);
%         elseif (i/cycles) ~= floor(i/cycles)
%             SurgeUndet(i,1) = SurgeUndet(i,1) * SurgeUndet(cycles*(m),1);
%         end
%     end
% end
% 
% SurgeTest = [];
% j = 1;
% for i = 1:length(SurgeUndet)
%     if i/cycles == floor(i/cycles)
%         SurgeTest(j,1) = i/cycles;
%         SurgeTest(j,2) = SurgeUndet(i,1);
%         j = j+1;
%     else
%         continue
%     end
% end
% 
% figure
% % scatter(parislaw(:,5), scatter_undet, 'MarkerFaceColor',[0 .75 .75])
% hold on
% plot(ParisSurge(:,5), ParisSurge(:,2))

%% Multiple Surges

p_surge = 80; % MPa
p_axialsurge = (80*0.1)/(2*0.01); % MPa

SurgeFrequency = 100; % 1 Surge per N cycles
SurgeOccurence = 1; % First surge during first cycle

% Surge Paris Law
MultiParisSurge = [];
cLengthMultiSurge = 0.0003;
SurgeCounter = 1;
i = 1;
t = 0.01;
C = 10^-12;
m = 4;

while cLengthMultiSurge <= 0.01
    if i == SurgeCounter
        Y = 0.728 + 0.373*(cLengthMultiSurge/t)^2 - 0.029*(cLengthMultiSurge/t)^4;
        K = Y*(p_axialsurge)*sqrt(pi*cLengthMultiSurge);
        MultiParisSurge(i,1) = cLengthMultiSurge; % initial crack length (mm)
        MultiParisSurge(i,2) = K; % stress intensity factor, MPa(m)^1/2
        MultiParisSurge(i,3) = C*(K)^m; % length increase per cycle (m)
        cLengthMultiSurge = cLengthMultiSurge + MultiParisSurge(i,3); % length after increase during cycle (mm)
        MultiParisSurge(i,4) = cLengthMultiSurge*10^3;
        MultiParisSurge(i,5) = i;
        i = i+1;
        SurgeCounter = SurgeCounter + SurgeFrequency;
    else
        Y = 0.728 + 0.373*(cLengthMultiSurge/t)^2 - 0.029*(cLengthMultiSurge/t)^4;
        K = Y*(sigma_max/10^6)*sqrt(pi*cLengthMultiSurge);
        MultiParisSurge(i,1) = cLengthMultiSurge; % initial crack length (mm)
        MultiParisSurge(i,2) = K; % stress intensity factor, MPa(m)^1/2
        MultiParisSurge(i,3) = C*(K)^m; % length increase per cycle (m)
        cLengthMultiSurge = cLengthMultiSurge + MultiParisSurge(i,3); % length after increase during cycle (mm)
        MultiParisSurge(i,4) = cLengthMultiSurge*10^3;
        MultiParisSurge(i,5) = i;
        i = i+1;
    end
end

%% Plotting Prob Undetection Against Ins Int

failmulti = [];
omult = 1;
for cyclesmulti = 8000:1:12000
    MultiSurgeUndet = [];
    MultiSurgeInspInt = [];
    for i = 1:length(MultiParisSurge)
        MultiSurgeUndet(i,1) = 1 - fn_pod(2055792,MultiParisSurge(i,1)*1000); % Probability of undetected crack
    end

    MultiSurgeInspectionIterations = [];
    for i = 1:length(MultiParisSurge)
        MultiSurgeInspInt(i,1) = MultiParisSurge(i,5)/cyclesmulti; % Inspections per N cycles;
        if rem(MultiSurgeInspInt(i,1),1) == 0
            MultiSurgeInspectionIterations = MultiSurgeInspInt(i,1);
        end
    end
    
    m = 0;
    for i = 1:length(MultiSurgeUndet)
        if (i/cyclesmulti) == floor(i/cyclesmulti) % if an inspection takes place
            m = m+1; % Number of inspections increase by 1
        end
    
        n = i/MultiSurgeInspectionIterations;
        if n > 1 && m > 1
            if (i/cyclesmulti) == floor(i/cyclesmulti)
                MultiSurgeUndet(i,1) = MultiSurgeUndet(i,1) * MultiSurgeUndet(cyclesmulti*(m-1),1);
            elseif (i/cyclesmulti) ~= floor(i/cyclesmulti)
                MultiSurgeUndet(i,1) = MultiSurgeUndet(i,1) * MultiSurgeUndet(cyclesmulti*(m),1);
            end
    
        elseif n < 1 && m < 1
            MultiSurgeUndet(i,1) = MultiSurgeUndet(i,1);
        elseif n >= 1 && m == 1
            if (i/cyclesmulti) == floor(i/cyclesmulti)
                MultiSurgeUndet(i,1) = MultiSurgeUndet(i,1);
            elseif (i/cyclesmulti) ~= floor(i/cyclesmulti)
                MultiSurgeUndet(i,1) = MultiSurgeUndet(i,1) * MultiSurgeUndet(cyclesmulti*m,1);
            end
        end
    end
    
    MultiSurgeTest = [];
    j = 1;
    for i = 1:length(MultiSurgeUndet)
        if i/cyclesmulti == floor(i/cyclesmulti)
            MultiSurgeTest(j,1) = i/cyclesmulti;
            MultiSurgeTest(j,2) = MultiSurgeUndet(i,1);
            j = j+1;
        else
            continue
        end
    end

    failmulti(omult,2) = MultiSurgeTest(end,2);
    failmulti(omult,1) = cyclesmulti;
    omult = omult + 1;
end

figure
scatter(MultiSurgeTest(:,1),MultiSurgeTest(:,2), 'MarkerFaceColor',[0 .75 .75])
hold on
yline(0.01)
title('Probability of Failure Against Inspection Interval')
xlabel('Inspection Interval')
ylabel('Probability of Failure')

%% Matrix for scatter plot
MultiSurgeScatter = [];
j = 1;
for i = 1:length(MultiSurgeUndet)
    if (i/cycles) == floor(i/cycles)
        MultiSurgeScatter(i,1) = MultiSurgeUndet(i,1);
    else
        MultiSurgeScatter(i,1) = nan;
    end
end

% Plotting for Multiple Surges Against Number of Cycles
figure
scatter(MultiParisSurge(:,5), MultiSurgeScatter,'MarkerFaceColor',[0 .75 .75])
%[scatterpoly,gof] = fit(insp_int,p_undet,'poly2');
% plot(scatterpoly,insp_int,scatter_undet);

figure
scatter(MultiSurgeTest(:,1), MultiSurgeTest(:,2), 'MarkerFaceColor',[0 .75 .75])

%% Creating Probability of Surges
% % Assuming probability of surge during every cycle is 0.5
% % Independent of outcome of previous cycles
% 
% cycles1 = 3000;
% 
% isSurge = 0;
% p_surge = 80; % MPa
% p_axialsurge = (80*0.1)/(2*0.01); % MPa
% 
% % Probability Surge Paris Law
% ProbParisSurge = [];
% cLengthProb = 0.0003;
% i = 1;
% t = 0.01;
% C = 10^-12;
% m = 4;
% 
% while cLengthProb <= 0.01
%     if rand < 0.5
%         isSurge = 0;
%     elseif rand >=0.5
%         isSurge = 1;
%     end
%     
%     if isSurge == 1
%         Y = 0.728 + 0.373*(cLengthProb/t)^2 - 0.029*(cLengthProb/t)^4;
%         K = Y*(p_axialsurge)*sqrt(pi*cLengthProb);
%         ProbParisSurge(i,1) = cLengthProb; % initial crack length (mm)
%         ProbParisSurge(i,2) = K; % stress intensity factor, MPa(m)^1/2
%         ProbParisSurge(i,3) = C*(K)^m; % length increase per cycle (m)
%         cLengthProb = cLengthProb + ProbParisSurge(i,3); % length after increase during cycle (mm)
%         ProbParisSurge(i,4) = cLengthProb*10^3;
%         ProbParisSurge(i,5) = i;
%         i = i+1;
%     elseif isSurge == 0
%         Y = 0.728 + 0.373*(cLengthProb/t)^2 - 0.029*(cLengthProb/t)^4;
%         K = Y*(sigma_max/10^6)*sqrt(pi*cLengthProb);
%         ProbParisSurge(i,1) = cLengthProb; % initial crack length (mm)
%         ProbParisSurge(i,2) = K; % stress intensity factor, MPa(m)^1/2
%         ProbParisSurge(i,3) = C*(K)^m; % length increase per cycle (m)
%         cLengthProb = cLengthProb + ProbParisSurge(i,3); % length after increase during cycle (mm)
%         ProbParisSurge(i,4) = cLengthProb*10^3;
%         ProbParisSurge(i,5) = i;
%         i = i+1;
%     end
% end
% 
% ProbSurgeUndet = [];
% for i = 1:length(ProbParisSurge)
%     ProbSurgeUndet(i,1) = 1 - fn_pod(2055792,ProbParisSurge(i,1)*1000); % Probability of undetected crack
% end
% 
% ProbSurgeInspInt = [];
% ProbSurgeInspectionIterations = [];
% for i = 1:length(ProbParisSurge)
%     ProbSurgeInspInt(i,1) = ProbParisSurge(i,5)/cycles1; % Inspections per 1000 cycles;
%     if rem(ProbSurgeInspInt(i,1),1) == 0
%         ProbSurgeInspectionIterations = ProbSurgeInspInt(i,1);
%     end
% end
% 
% m = 0;
% for i = 1:length(ProbSurgeUndet)
%     if (i/cycles1) == floor(i/cycles1) % if an inspection takes place
%         m = m+1; % Number of inspections increase by 1
%     end
% 
%     n = i/ProbSurgeInspectionIterations;
%     if n > 1 && m > 1
%         if (i/cycles1) == floor(i/cycles1)
%             ProbSurgeUndet(i,1) = ProbSurgeUndet(i,1) * ProbSurgeUndet(cycles1*(m-1),1);
%         elseif (i/cycles1) ~= floor(i/cycles1)
%             ProbSurgeUndet(i,1) = ProbSurgeUndet(i,1) * ProbSurgeUndet(cycles1*(m),1);
%         end
% 
%     elseif n < 1 && m < 1
%         ProbSurgeUndet(i,1) = ProbSurgeUndet(i,1);
%     elseif n >= 1 && m == 1
%         if (i/cycles1) == floor(i/cycles1)
%             ProbSurgeUndet(i,1) = ProbSurgeUndet(i,1);
%         elseif (i/cycles1) ~= floor(i/cycles1)
%             ProbSurgeUndet(i,1) = ProbSurgeUndet(i,1) * ProbSurgeUndet(cycles1*(m),1);
%         end
%     end
% end
% 
% ProbSurgeTest = [];
% j = 1;
% for i = 1:length(ProbSurgeUndet)
%     if i/cycles1 == floor(i/cycles1)
%         ProbSurgeTest(j,1) = i/cycles1;
%         ProbSurgeTest(j,2) = ProbSurgeUndet(i,1);
%         j = j+1;
%     else
%         continue
%     end
% end
% 
% % Matrix for scatter plot
% ProbSurgeScatter = [];
% j = 1;
% for i = 1:length(ProbSurgeUndet)
%     if (i/cycles1) == floor(i/cycles1)
%         ProbSurgeScatter(i,1) = ProbSurgeUndet(i,1);
%     else
%         ProbSurgeScatter(i,1) = nan;
%     end
% end
% 
% % Plotting for Probability Surges Against Number of Cycles
% figure
% scatter(ProbParisSurge(:,5), ProbSurgeScatter,'MarkerFaceColor',[0 .75 .75])
% %[scatterpoly,gof] = fit(insp_int,p_undet,'poly2');
% % plot(scatterpoly,insp_int,scatter_undet);
% 
% figure
% scatter(ProbSurgeTest(:,1), ProbSurgeTest(:,2), 'MarkerFaceColor',[0 .75 .75])

