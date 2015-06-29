function PlotEKFResults(EKFData,measurementData)

getData     = load(EKFData);
getData2    = load(measurementData);

KF_OUT      = [getData.EKF.data];
measData    = [getData2.AO.data];

time        = KF_OUT(:,1);
predAO      = KF_OUT(:,2:7);
states      = KF_OUT(:,8:31);
covs        = sqrt(KF_OUT(:,38:61));
innovations = KF_OUT(:,62:67);
S           = sqrt(KF_OUT(:,68:73));


name_string   = {'CX_dE1'; 'CX_dE2';...
    'CY_dR1'; 'CY_dR2';'CY_dE1'; ...
    'CZ_dE1'; 'CZ_dE2'; 'CZ_dA1';...
    'Cl_dA1'; 'Cl_dA2'; 'Cl_dA3';...
    'Cl_dR1'; 'Cl_dR2'; 'Cl_dE1';...
    'Cm_dE1'; 'Cm_dE2'; 'Cm_dE3'; 'Cm_dA1';...
    'Cn_dA1'; 'Cn_dA2'; 'Cn_dR1';...
    'Cn_dR2'; 'Cn_dE1';'Cn_dE2'};

% figure
% subplot(3,1,1)
% plot(measData(:,1),measData(:,2),'r')
% hold on
% plot(time(:,1),predAO(:,1),'b')
% legend('Measured','Predicted')
% title('ax')
% hold on
% subplot(3,1,2)
% plot(measData(:,1),measData(:,3),'r')
% hold on
% plot(time(:,1),predAO(:,2),'b')
% legend('Measured','Predicted')
% title('ay')
% subplot(3,1,3)
% plot(measData(:,1),measData(:,4),'r')
% hold on
% plot(time(:,1),predAO(:,3),'b')
% legend('Measured','Predicted')
% title('az')
%
% figure
% subplot(3,1,1)
% plot(measData(:,1),measData(:,5),'r')
% hold on
% plot(time(:,1),predAO(:,4),'b')
% legend('Measured','Predicted')
% title('P')
% hold on
% subplot(3,1,2)
% plot(measData(:,1),measData(:,6),'r')
% hold on
% plot(time(:,1),predAO(:,5),'b')
% legend('Measured','Predicted')
% title('Q')
% subplot(3,1,3)
% plot(measData(:,1),measData(:,7),'r')
% hold on
% plot(time(:,1),predAO(:,6),'b')
% legend('Measured','Predicted')
% title('R')

%Plot innovations
figure
subplot(3,1,1)
hold on
plot(time,innovations(:,1),'b')
hold on
plot(time,2*S(:,1),'-r')
hold on
plot(time,-2*S(:,1),'-r')
title('Accelerations')
hold on
subplot(3,1,2)
hold on
plot(time,innovations(:,2),'b')
hold on
plot(time,2*S(:,2),'-r')
hold on
plot(time,-2*S(:,2),'-r')
hold on
subplot(3,1,3)
hold on
plot(time,innovations(:,3),'b')
hold on
plot(time,2*S(:,3),'-r')
hold on
plot(time,-2*S(:,3),'-r')


figure
subplot(3,1,1)
hold on
plot(time,innovations(:,4),'b')
hold on
plot(time,2*S(:,4),'-r')
hold on
plot(time,-2*S(:,4),'-r')
title('Angular Rates')
hold on
subplot(3,1,2)
hold on
plot(time,innovations(:,5),'b')
hold on
plot(time,2*S(:,5),'-r')
hold on
plot(time,-2*S(:,5),'-r')
hold on
subplot(3,1,3)
hold on
plot(time,innovations(:,6),'b')
hold on
plot(time,2*S(:,6),'-r')
hold on
plot(time,-2*S(:,6),'-r')

for i = 1:24
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
plot(time,states)


