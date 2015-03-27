function PlotEKFResults(EKFData,measurementData)

getData  = load(EKFData);
getData2 = load(measurementData);

KF_OUT   = [getData.EKF.data];
measData = [getData2.AO.data];

time    = KF_OUT(:,1);
predAO  = KF_OUT(:,2:7);
states  = KF_OUT(:,8:20);
covs    = sqrt(KF_OUT(:,27:39));
nu      = KF_OUT(:,40:45);
S = KF_OUT(:,46:end); 
 

name_string   = {'CX_dE';...
'CY_dR'; 'CY_dE';...
'CZ_dE'; 'CZ_dA';...
'Cl_dA'; 'Cl_dR'; 'Cl_dE';...
'Cm_dE'; 'Cm_dA';...
'Cn_dA'; 'Cn_dR'; 'Cn_dE'};

figure
subplot(3,1,1)
plot(measData(:,1),measData(:,2),'r')
hold on
plot(time(:,1),predAO(:,1),'b')
legend('Measured','Predicted')
title('ax')
hold on
subplot(3,1,2)
plot(measData(:,1),measData(:,3),'r')
hold on
plot(time(:,1),predAO(:,2),'b')
legend('Measured','Predicted')
title('ay')
subplot(3,1,3)
plot(measData(:,1),measData(:,4),'r')
hold on
plot(time(:,1),predAO(:,3),'b')
legend('Measured','Predicted')
title('az')

figure
subplot(3,1,1)
plot(measData(:,1),measData(:,5),'r')
hold on
plot(time(:,1),predAO(:,4),'b')
legend('Measured','Predicted')
title('P')
hold on
subplot(3,1,2)
plot(measData(:,1),measData(:,6),'r')
hold on
plot(time(:,1),predAO(:,5),'b')
legend('Measured','Predicted')
title('Q')
subplot(3,1,3)
plot(measData(:,1),measData(:,7),'r')
hold on
plot(time(:,1),predAO(:,6),'b')
legend('Measured','Predicted')
title('R')

figure
subplot(3,1,1)
plot(time,nu(:,1),'b')
hold on
plot(time,2*sqrt(S(:,1)),'-r')
hold on
plot(time,-2*sqrt(S(:,1)),'-r')
subplot(3,1,2)
plot(time,nu(:,2),'b')
hold on
plot(time,2*sqrt(S(:,2)),'-r')
hold on
plot(time,-2*sqrt(S(:,2)),'-r')
subplot(3,1,3)
plot(time,nu(:,3),'b')
hold on
plot(time,2*sqrt(S(:,3)),'-r')
hold on
plot(time,-2*sqrt(S(:,3)),'-r')

figure
subplot(3,1,1)
plot(time,nu(:,4),'b')
hold on
plot(time,2*sqrt(S(:,4)),'-r')
hold on
plot(time,-2*sqrt(S(:,4)),'-r')
subplot(3,1,2)
plot(time,nu(:,5),'b')
hold on
plot(time,2*sqrt(S(:,5)),'-r')
hold on
plot(time,-2*sqrt(S(:,5)),'-r')
subplot(3,1,3)
plot(time,nu(:,6),'b')
hold on
plot(time,2*sqrt(S(:,6)),'-r')
hold on
plot(time,-2*sqrt(S(:,6)),'-r')



for i = 1:13
    errors(:,i)= 1-states(:,i);% - aeroCoeffs(i);
    figure
    plot(time,errors(:,i),'b')
    hold on
    plot(time,2*covs(:,i),'--r')
    hold on
    plot(time,-2*covs(:,i),'--r')    
    legend(name_string{i})
end

figure
% plot(time,states)


