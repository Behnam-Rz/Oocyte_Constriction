%% Please Don't change the following three lines 
clear 
close all 
clc 
%% Parameters 
%%% In the following lines of codes the slopes of the critical pressure
%%% with respect to surface tension (m_tau) and elasticity (m_E) are
%%% studied for differnt values of d (constriction diameter) 

%%%%%%%%%%%%%%%%%% Study Parametrs (You can modify these) %%%%%%%%%%%%%%%%%%
d_min = 40;                     % Min constriction diameter          
d_max = 70;                     % Max constriction diameter 

%%%%%%%%%%%%%% Other Parameters (Based on fitted results) %%%%%%%%%%%%%%%%%%
th = 24.4;              % Effective thickness 
alpha = 9*pi/180;       % Entrance angle   
D0 = 74.4;              % Coell diameter 

E = 1.5;                % Elasticity (not applicable for calculating the slopes) 
gamma = 1;              % Surface Tension (not applicable for calculating the slopes) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Please Don't change the following Lines 
dvec = [d_min : (d_max-d_min)/100 : d_max]';
m_E = zeros(length(dvec),1);            % Initializing m_E vector 
m_tau = zeros(length(dvec),1);          % Initializing m_tau vector 
P_cr = zeros(length(dvec),1);           % Initializing Pcr vector

for ii = 1:length(dvec)
    d = dvec(ii); 
    DL_cr = cos(alpha)*(sin(alpha)*(4*D0^3 - (2-cot(alpha))*d^3)/(1+sin(alpha))^2)^(1/3);
    A_cr = pi*DL_cr*DL_cr/4*(1 + (1+sin(alpha))^2/cos(alpha)/cos(alpha) + 1/tan(alpha)) + pi*d^2/4*(2 - 1/tan(alpha));
    slope_tau0 = 10*2*(2/d - 2*cos(alpha)/DL_cr); 
    slope_E = slope_tau0*th*(A_cr/(pi*D0^2) - 1); 
    Pcr = E*slope_E + gamma*slope_tau0; 

    %%% Outputs 
    m_E(ii) = slope_E; 
    m_tau0(ii) = slope_tau0; 
    P_cr(ii) = Pcr; 
end

%% Plotting 
plot(dvec,m_E,'r','LineWidth',2)
hold on 
plot(dvec,m_tau0,'b','LineWidth',2)
legend('m_E','m_\tau')
