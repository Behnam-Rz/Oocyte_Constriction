clear 
close all 
clc 
%% Parameters 
%%% In the following lines of codes the slopes of the critical pressure
%%% with respect to surface tension (m_tau) and elasticity (m_E) are
%%% studied for differnt values of alpha (entrance angle) 

%%%%%%%%%%%%%%%%%% Study Parametrs (You can modify these) %%%%%%%%%%%%%%%%%%
alpha_min = 5;                % Min entrance angle 
alpha_max = 90;               % Max entrance angle 

%%%%%%%%%%%%%% Other Parameters (Based on fitted results) %%%%%%%%%%%%%%%%%%
th = 24.4;              % Effective thickness 
d = 50;                 % Constriction diameter   
D0 = 74.4;              % Cell diameter 

E = 1.5;                % Elasticity (not applicable for calculating the slopes) 
gamma = 1;              % Surface Tension (not applicable for calculating the slopes) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% m_E and m_tau0 vs alpha 
alphavec = [alpha_min : (alpha_max-alpha_min)/100 : alpha_max]'; 
m_E = zeros(length(alphavec),1);            % Initializing m_E vector 
m_tau = zeros(length(alphavec),1);          % Initializing m_tau vector 
P_cr = zeros(length(alphavec),1);           % Initializing Pcr vector

for ii = 1:length(alphavec)
    alpha = alphavec(ii)*pi/180; 
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
plot(alphavec,m_E,'r','LineWidth',2)
hold on 
plot(alphavec,m_tau0,'b','LineWidth',2)
