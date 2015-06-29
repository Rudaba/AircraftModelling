function [states_OUT, cov_OUT] = simulateIMMEKF(accelerations,omega,ALP,BET,u,Vt,Height)

persistent IMM_states IMM_cov transMatrix modeProbs Q R

%Intialise filter
if isempty(IMM_states)
    
     [IMM_states, IMM_cov, transMatrix, modeProbs, Q, R] = initialiseIMMEKF(accelerations,omega)
    
end

%Determine number of models
m = size(IMM_states,2);

%Number of states
n = size(IMM_states,1);

%*****Interacting/Mixing*****
c_j = zeros(1,m);
for j = 1:m
    for i = 1:m
        c_j(j) = c_j(j) + transMatrix(i,j)*modeProbs(i);
    end
end

modeProbs_ij = zeros(m,m);
for i = 1:m
    for j = 1:m
        modeProbs_ij(i,j) = 1/c_j(j)*transMatrix(i,j)*modeProbs(i);
    end
end

%Calculate mixed state and covariance
X0j = zeros(n,m);
for j = 1:m
    for i = 1:m
        X0j(:,j) = X0j(:,j) + IMM_states(:,i)*modeProbs_ij(i,j); 
    end
end

P0j = zeros(n,n,m);
for j = 1:m
    for i = 1:m
        P0j(:,:,j) = P0j(:,:,j) + modeProbs_ij(i,j)*(IMM_cov(:,:,i) + (IMM_states(:,i)-X0j(:,j))*(IMM_states(:,i)-X0j(:,j))'); 
    end
end

%Perform Kalman filter updates on each filter and calculate likelihoods
for i = 1:m
    [IMM_states(:,i),IMM_cov(:,:,i), measurementPred, S, innovation]  = IMM_EKFPredictAndUpdate(X0j(:,i),P0j(:,:,i),i,accelerations,omega,Q(:,:,i),R,ALP,BET,u,Vt,Height);
                                                                                                
    lamda(i)                               = gaussPDF(measurement, measurementPred, S); 
end

%Mode probability update
c = 0;
for j = 1:m
    c = c + lamda(j) * c_j(j);
end

for j = 1:m
    modeProbs(j) = 1/c * lamda(j) * c_j(j); 
end

%Estimate and Covariance Combination
X_Filter = zeros(n,1);
for j = 1:m
    X_Filter = X_Filter + IMM_states(:,j) * modeProbs(j);
end

P_Filter = zeros(n,n);
for j = 1:m
    P_Filter = P_Filter +  modeProbs(j) * (IMM_cov(:,:,j) + (IMM_states(:,j)-X_Filter)*(IMM_states(:,j)-X_Filter)');
end