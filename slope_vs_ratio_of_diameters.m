close all 
clear 
clc
%% Constant Parameters during the whole simulations
th = 24.4;          % The effective thickness 

%% slopes of (E and \gamma_0) vs D0/d
%%%%%%%%%%%%%%% Don't change the following parameters %%%%%%%%%%%%%%%%%%
alpha = 9*pi/180;       %

D0_d_min = 1.2; 
D0_d_max = 2.0; 
D0_d_vec = D0_d_min : (D0_d_max - D0_d_min)/100 :  D0_d_max; 

D0vec = [65 75 85 105 125]; 
colorvec = {'k','b','r','g','m'}; 

figure (4); 
set(gcf,'Units','inches')
set(gcf,'Color','white')

for jj = 1 : length(D0vec)
    D0 = D0vec(jj); 
    for ii = 1:length(D0_d_vec)
        d = D0/D0_d_vec(ii); 
        DL_cr = cos(alpha)*(sin(alpha)*(4*D0^3 - (2-cot(alpha))*d^3)/(1+sin(alpha))^2)^(1/3);
        A_cr = pi*DL_cr*DL_cr/4*(1 + (1+sin(alpha))^2/cos(alpha)/cos(alpha) + 1/tan(alpha)) + pi*d^2/4*(2 - 1/tan(alpha));
        slope_gamma = 10*2*(2/d - 2*cos(alpha)/DL_cr); 
        slope_E = slope_gamma*th*(A_cr/(pi*D0^2) - 1); 
    
        %%% Outputs 
        m_E(ii) = slope_E; 
        m_tau(ii) = slope_gamma; 
    end
    
    %%%%%%%%%%%% m_tau vs D0 %%%%%%%%%%%%%
    subplot(1,2,1); 
    box off 
    hold on 
    set ( gca, 'xdir', 'reverse' )
    
    plot(D0_d_vec,m_tau,'color',colorvec{jj},'LineWidth',3)
    
    xlabel('D_0/d')
    ylabel('m_{\tau} (bar.m/N)')

    if jj == length(D0vec)
        legend(strcat('D_0 =',num2str(D0vec(1))),strcat('D_0 =',num2str(D0vec(2))),strcat('D_0 =',num2str(D0vec(3))),strcat('D_0 =',num2str(D0vec(4))),strcat('D_0 =',num2str(D0vec(5))),'Location','northeast','edgecolor','none')
    end
    
    %%%%%%%%%%%% m_E vs D0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,2); 
    box off 
    hold on 
    set ( gca, 'xdir', 'reverse' )
    
    plot(D0_d_vec,m_E,'color',colorvec{jj},'LineWidth',3)
    
    xlabel('D_0/d')
    ylabel('m_E (mbar/kPa)')

    if jj == length(D0vec)
        legend(strcat('D_0 =',num2str(D0vec(1))),strcat('D_0 =',num2str(D0vec(2))),strcat('D_0 =',num2str(D0vec(3))),strcat('D_0 =',num2str(D0vec(4))),strcat('D_0 =',num2str(D0vec(5))),'Location','northeast','edgecolor','none')
    end
end
