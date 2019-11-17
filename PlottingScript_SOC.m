clc
clear all
close all

%load SimData_SOC_Pres_2018.02.13_Olga.mat
%load SimData_SOC_GOR_2018.02.14_Olga.mat
load SimData_SOC2_Nominal_2019.03.01.mat


plotStates = 0;
plotInputs = 1;
plotOutputs = 1;
plotMonitors = 1;


col1 = [8, 118, 184]/255;
col2 = [36, 162, 14]/255;
col3 = [233, 11, 0]/255;
u = u_SOC;
nU = length(u(1,:));

NIT = length(x(:,1));
dt = 100;
time = simTime/3600;%(0:dt:(NIT-1)*dt)/3600;

if(plotStates)
    nX = length(x(1,:));
    for ix=1:nX
        figure('Name',states{ix})
        clf
        plot(time,x(:,ix),'LineWidth',1)
        legend(states{ix})
        hold on
        plot(time,xMin(ix)*ones(NIT,1),'-r','LineWidth',2)
        plot(time,xMax(ix)*ones(NIT,1),'-r','LineWidth',2)
        plot(time,xOptimal(ix)*ones(NIT,1),'--','Color',[0 .5 0],'LineWidth',2)
    end
end

figure('Name','SOC')
clf
subplot(2,1,1)
plot(time,y1,time,r1*ones(size(time)),'LineWidth',1)
legend('c1','r1')
subplot(2,1,2)
plot(time,y2,time,r2*ones(size(time)),'LineWidth',1)
legend('c2','r2')

figure('Name','J')
clf
subplot(2,1,1)
plot(time,J_SOC,'LineWidth',1)
legend('J')
subplot(2,1,2)
plot(time,Ju_SOC(:,1),'LineWidth',1)
hold on
plot(time,Ju_SOC(:,2),'LineWidth',1)
legend('Ju1','Ju2')

if(plotInputs)
    nu = length(u(1,:));
    for iu=1:nu
        figure('Name',controlTags{iu})
        clf
        plot(time,u(:,iu),'LineWidth',1)
        legend(controlTags{iu})
        hold on
        %plot(time,uMin(iu)*ones(NIT,1),'-r','LineWidth',2)
        %plot(time,uMax(iu)*ones(NIT,1),'-r','LineWidth',2)
        plot(time,uOptimal(iu)*ones(NIT,1),'--','Color',[0 .5 0],'LineWidth',1)
    end
end

if(plotOutputs)
    ny = length(yOlgaMeas(1,:));
    idx = [2 3 1]+nu-3;
    for iy=1:ny
        figure('Name',measurementTags{iy})
        clf
        plot(time,yOlgaMeas(:,iy),'LineWidth',1)
        hold on
        plot(time,yKfMeas(:,iy),'--','Color',col3,'LineWidth',1)
        legend(measurementTags{iy},measurementTagsKf{iy})
        if(iy>(ny-3))
            plot(time,u(:,idx(iy-(ny-3))),'--k','LineWidth',1)
        end
    end
end

if(plotMonitors)
    ny = length(yOlgaMoni(1,:));
    for iy=1:ny
        figure('Name',measurementTagsMonitor{iy})
        clf
        plot(time,yOlgaMoni(:,iy),'LineWidth',1)
        hold on
        plot(time,yKfMoni(:,iy),'--','Color',col2,'LineWidth',1)
        legend(measurementTagsMonitor{iy},measurementTagsMonitorKf{iy})
    end
end