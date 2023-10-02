function [] = Visualization(CSlip_Gps,CSlip_Bds,MP_GpsSequence,MP_BdsSequence,allGpsData,allBdsData)%designed by hzLiu,2023.7.14
clf;
%% Cycle-Slip : Single PRN TimeSequence
for i=1:length(CSlip_Gps)
    clf;
    figure(1);
    PRN=CSlip_Gps{i,1}(1,1);
    if PRN<10
        Sat=strcat('G0',int2str(PRN));
    else
        Sat=strcat('G',int2str(PRN));
    end
    name=strcat(Sat,'-','Wide-Lane AMB Time Sequence');
    set(gcf,'Position',[50 50 1000 560])
    %MW: Wide-Lane AMB of Healthy Epoch
    scatter(CSlip_Gps{i,1}(CSlip_Gps{i,1}(:,3)==0&(abs(CSlip_Gps{i,1}(:,4))<1000),2),CSlip_Gps{i,1}(CSlip_Gps{i,1}(:,3)==0&(abs(CSlip_Gps{i,1}(:,4))<1000),4),6,[0.28 0.57 0.54],"filled");hold on
    %MW: Wide-Lane AMB of Epoch with Cycle-Slips
    scatter(CSlip_Gps{i,1}(CSlip_Gps{i,1}(:,3)==1&(abs(CSlip_Gps{i,1}(:,4))<1000),2),CSlip_Gps{i,1}(CSlip_Gps{i,1}(:,3)==1&(abs(CSlip_Gps{i,1}(:,4))<1000),4),15,[0.73 0.47 0.58],"filled");hold on
    %MW: Wide-Lane AMB of Epoch with Outliers
    scatter(CSlip_Gps{i,1}(CSlip_Gps{i,1}(:,3)==2&(abs(CSlip_Gps{i,1}(:,4))<1000),2),CSlip_Gps{i,1}(CSlip_Gps{i,1}(:,3)==2&(abs(CSlip_Gps{i,1}(:,4))<1000),4),15,[0.62 0.49 0.31],"filled");hold on

    legend1=legend('$\bf{Healthy Epoch}$','$\bf{CycleSlips}$','$\bf{Outliers}$','interpreter','latex','FontSize',10.5);
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',10.5);
    xlabel('$\bf{Time(s)}$','interpreter','latex','FontSize', 17)
    ylabel('$\bf{Value}$','interpreter','latex','FontSize', 17)
    title({'$\bf{WideLane-AMB-TimeSequence(ValidEpoch)}$'}, 'interpreter','latex','FontSize', 19);
    set(gca,'linewidth',1.5);
    set(gca,'XGrid','on','XMinorGrid','off','YGrid','on','YMinorGrid','off');
    set(gca,'fontsize',13,'fontname','Times','FontWeight','bold')
    box on
    cd ..\img\Gps\
    saveas(gcf, name, 'png');
    cd ..\..\code\
    hold off
end
for i=1:length(CSlip_Bds)
    clf;
    figure(1);
    PRN=CSlip_Bds{i,1}(1,1);
    if PRN<10
        Sat=strcat('C0',int2str(PRN));
    else
        Sat=strcat('C',int2str(PRN));
    end
    name=strcat(Sat,'-','Wide-Lane AMB Time Sequence');
    set(gcf,'Position',[50 50 1000 560])
    %MW: Wide-Lane AMB of Healthy Epoch
    scatter(CSlip_Bds{i,1}(CSlip_Bds{i,1}(:,3)==0&(abs(CSlip_Bds{i,1}(:,4))<1000),2),CSlip_Bds{i,1}(CSlip_Bds{i,1}(:,3)==0&(abs(CSlip_Bds{i,1}(:,4))<1000),4),6,[0.28 0.57 0.54],"filled");hold on
    %MW: Wide-Lane AMB of Epoch with Cycle-Slips
    scatter(CSlip_Bds{i,1}(CSlip_Bds{i,1}(:,3)==1&(abs(CSlip_Bds{i,1}(:,4))<1000),2),CSlip_Bds{i,1}(CSlip_Bds{i,1}(:,3)==1&(abs(CSlip_Bds{i,1}(:,4))<1000),4),15,[0.73 0.47 0.58],"filled");hold on
    %MW: Wide-Lane AMB of Epoch with Outliers
    scatter(CSlip_Bds{i,1}(CSlip_Bds{i,1}(:,3)==2&(abs(CSlip_Bds{i,1}(:,4))<1000),2),CSlip_Bds{i,1}(CSlip_Bds{i,1}(:,3)==2&(abs(CSlip_Bds{i,1}(:,4))<1000),4),15,[0.62 0.49 0.31],"filled");hold on

    legend1=legend('$\bf{Healthy Epoch}$','$\bf{CycleSlips}$','$\bf{Outliers}$','interpreter','latex','FontSize',10.5);
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',10.5);
    xlabel('$\bf{Time(sec)}$','interpreter','latex','FontSize', 17)
    ylabel('$\bf{Value(cycle}$','interpreter','latex','FontSize', 17)
    title({'$\bf{WideLane-AMB-TimeSequence(ValidEpoch)}$'}, 'interpreter','latex','FontSize', 19);
    set(gca,'linewidth',1.5);
    set(gca,'XGrid','on','XMinorGrid','off','YGrid','on','YMinorGrid','off');
    set(gca,'fontsize',13,'fontname','Times','FontWeight','bold')

    box on
    cd ..\img\Bds\
    saveas(gcf, name, 'png');
    cd ..\..\code\
    hold off
end
%% SNR : Single PRN TimeSequence
for i=1:length(allGpsData)
    clf;
    figure(1);
    PRN=allGpsData{i,1}.SatelliteID(1,1);
    if PRN<10
        Sat=strcat('G0',int2str(PRN));
    else
        Sat=strcat('G',int2str(PRN));
    end
    name=strcat(Sat,'-','L1L2L5 SNR Time Sequence');
    set(gcf,'Position',[50 50 1000 560])   
    
    %SNR: L1
    plot(allGpsData{i,1}.Time,allGpsData{i,1}.S1C,'Color',[0.28 0.57 0.54],LineWidth=2);hold on
    %SNR: L2
    plot(allGpsData{i,1}.Time,allGpsData{i,1}.S2W,'Color',[0.73 0.47 0.58],LineWidth=2);hold on
    %SNR: L5
    plot(allGpsData{i,1}.Time,allGpsData{i,1}.S5X,'Color',[0.62 0.49 0.31],LineWidth=2);hold on

    legend1=legend('$\bf{L1:S1C}$','$\bf{L2:S2W}$','$\bf{L5:S5X}$','interpreter','latex','FontSize',10.5);
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',10.5);
    xlabel('$\bf{Time(sec)}$','interpreter','latex','FontSize', 17)
    ylabel('$\bf{SNR(dBHz)}$','interpreter','latex','FontSize', 17)
    title({'$\bf{L1/L2/L5-SNR-Time Sequence}$'}, 'interpreter','latex','FontSize', 19);
    set(gca,'linewidth',1.5);
    set(gca,'XGrid','on','XMinorGrid','off','YGrid','on','YMinorGrid','off');
    set(gca,'fontsize',13,'fontname','Times','FontWeight','bold')
    box on
    cd ..\img\Gps\
    saveas(gcf, name, 'png');
    cd ..\..\code\
    hold off
end
for i=1:length(allBdsData)
    clf;
    figure(1);
    PRN=allBdsData{i,1}.SatelliteID(1,1);
    if PRN<10
        Sat=strcat('C0',int2str(PRN));
    else
        Sat=strcat('C',int2str(PRN));
    end
    name=strcat(Sat,'-','B1B2B3 SNR Time Sequence');
    set(gcf,'Position',[50 50 1000 560])
    %SNR: B1
    plot(allBdsData{i,1}.Time,allBdsData{i,1}.S2I,'Color',[0.28 0.57 0.54],LineWidth=2);hold on
    %SNR: B2
    plot(allBdsData{i,1}.Time,allBdsData{i,1}.S7I,'Color',[0.62 0.49 0.31],LineWidth=2);hold on
    %SNR: B3
    plot(allBdsData{i,1}.Time,allBdsData{i,1}.S6I,'Color',[0.73 0.47 0.58],LineWidth=2);hold on

    legend1=legend('$\bf{B1:S2I}$','$\bf{B2:S7I}$','$\bf{B3:S6I}$','interpreter','latex','FontSize',10.5);
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',10.5);
    xlabel('$\bf{Time(sec)}$','interpreter','latex','FontSize', 17)
    ylabel('$\bf{SNR(dBHz)}$','interpreter','latex','FontSize', 17)
    title({'$\bf{B1/B2/B3-SNR-Time Sequence}$'}, 'interpreter','latex','FontSize', 19);
    set(gca,'linewidth',1.5);
    set(gca,'XGrid','on','XMinorGrid','off','YGrid','on','YMinorGrid','off');
    set(gca,'fontsize',13,'fontname','Times','FontWeight','bold')
    box on
    cd ..\img\Bds\
    saveas(gcf, name, 'png');
    cd ..\..\code\
    hold off
end

%% Multipath Effects : Single PRN TimeSequence
for i=1:length(MP_GpsSequence)
    clf;
    figure(1);
    PRN=MP_GpsSequence{i,1}(1,1);
    if PRN<10
        Sat=strcat('G0',int2str(PRN));
    else
        Sat=strcat('G',int2str(PRN));
    end
    name=strcat(Sat,'-','Multipath Effects Time Sequence');
    set(gcf,'Position',[50 50 1000 560])
    
    %Multipath: L1
    plot(MP_GpsSequence{i,1}(:,2),MP_GpsSequence{i,1}(:,3),'Color',[0.28 0.57 0.54],LineWidth=2);hold on
    %Multipath: L2
    plot(MP_GpsSequence{i,1}(:,2),MP_GpsSequence{i,1}(:,4),'Color',[0.73 0.47 0.58],LineWidth=2);hold on

    legend1=legend('$\bf{L1:Multipath}$','$\bf{L2:Multipath}$','interpreter','latex','FontSize',10.5);
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',10.5);
    xlabel('$\bf{Time(sec)}$','interpreter','latex','FontSize', 17)
    ylabel('$\bf{Value(m)}$','interpreter','latex','FontSize', 17)
    title({'$\bf{Multipath Effects-Time Sequence(HealthyEpoch)}$'}, 'interpreter','latex','FontSize', 19);
    set(gca,'linewidth',1.5);
    set(gca,'XGrid','on','XMinorGrid','off','YGrid','on','YMinorGrid','off');
    set(gca,'fontsize',13,'fontname','Times','FontWeight','bold')
    box on
    cd ..\img\Gps\
    saveas(gcf, name, 'png');
    cd ..\..\code\
    hold off
end
for i=1:length(MP_BdsSequence)
    clf;
    figure(1);
    PRN=MP_BdsSequence{i,1}(1,1);
    if PRN<10
        Sat=strcat('C0',int2str(PRN));
    else
        Sat=strcat('C',int2str(PRN));
    end
    name=strcat(Sat,'-','Multipath Effects Time Sequence');
    set(gcf,'Position',[50 50 1000 560])
    
    %Multipath: B1
    plot(MP_BdsSequence{i,1}(:,2),MP_BdsSequence{i,1}(:,3),'Color',[0.28 0.57 0.54],LineWidth=2);hold on
    %Multipath: B3
    plot(MP_BdsSequence{i,1}(:,2),MP_BdsSequence{i,1}(:,4),'Color',[0.73 0.47 0.58],LineWidth=2);hold on

    legend1=legend('$\bf{B1:Multipath}$','$\bf{B3:Multipath}$','interpreter','latex','FontSize',10.5);
    set(legend1,'LineWidth',1,'Interpreter','latex','FontSize',10.5);
    xlabel('$\bf{Time(sec)}$','interpreter','latex','FontSize', 17)
    ylabel('$\bf{Value(m)}$','interpreter','latex','FontSize', 17)
    title({'$\bf{Multipath Effects-Time Sequence(HealthyEpoch)}$'}, 'interpreter','latex','FontSize', 19);
    set(gca,'linewidth',1.5);
    set(gca,'XGrid','on','XMinorGrid','off','YGrid','on','YMinorGrid','off');
    set(gca,'fontsize',13,'fontname','Times','FontWeight','bold')
    box on
    cd ..\img\Bds\
    saveas(gcf, name, 'png');
    cd ..\..\code\
    hold off
end

end
