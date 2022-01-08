function makeFigures(I,C,H,InfCurv,HD,l)
%Show results by figures

figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
plot(0:l,I,'LineWidth',1.5)
xlim([0 l+1])
xlabel("Days")
ylabel("Counts")
title("Number of infected individuals")

subplot(2,2,2)
hold on
plot(C(2:end),'LineWidth',1.5)
plot(H(2:end),'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250])
plot(C(2:end)+H(2:end),'LineWidth',1.5,'Color',[0.8500 0.3250 0.0980])
xlim([0 l+1])
legend("Community","Household","Total","Location","northwest")
xlabel("Days")
ylabel("Counts")
title("Number of transmissions")
hold off

subplot(2,2,3)
CL = [0,0,0;0,0.447000000000000,0.741000000000000;0.850000000000000,0.325000000000000,0.0980000000000000;0.929000000000000,0.694000000000000,0.125000000000000;0.494000000000000,0.184000000000000,0.556000000000000;0.466000000000000,0.674000000000000,0.188000000000000;0.301000000000000,0.745000000000000,0.933000000000000;0.635000000000000,0.0780000000000000,0.184000000000000];
hold on
for i = 1:8
    stairs(InfCurv(:,i),'LineWidth',1.5,'Color',CL(i,:))
end
hold off
xlim([0 l])
xlabel("Days")
ylabel("Probability of remaining uninfected")
title("Probability of remaining uninfected with respect to household sizes")
lg=legend("1","2","3","4","5","6","7","8","Location","southwest");
title(lg,'Size')

subplot(2,2,4)
bp = bar(HD,'hist');
xlim([0 9])
bp.FaceColor = [0.9290 0.6940 0.1250];
bp.FaceAlpha = 0.618;
xlabel("Household size")
ylabel("Percentage")
title("Household size distribution")

end

