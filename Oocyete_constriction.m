close all 
clear 
clc 
%% Parameters 
%%%%%%%%%%%%%%%%%%%%% Model Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
D0 = 74.4;                  % Initial diameter of the cell (um)    


E = 1.5;                    % The Elasticity (kPa) 
tau0 = 1;                   % The surface tension (nN/um)
th = 24.4;                  % Mean shell thickness (um)


d = 50;                     % The diameter of the constriction (um)                                    
alpha = 9;                  % Inclination angle (degrees)   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model Ouputs 
% Pres_mbar: The pressure vector (mbar)
% Xf_shell: The vector of the coordinates of the front of the shell Xf
% X_rear: The vector of the coordinates of the rear of the shell Xr
% Width: The width vector (Width = Xf - Xr) 
% Height: The height vector (Height = 2*Rr)
% Aspect: The aspect ratio vector (AR = (Xf - Xr)/(2*Rr))

[Pres_mbar, Xf_shell, X_rear, Width, Height, Aspect] = model(D0,E,tau0,th,d,alpha);

%% Plotting Aspect ratio vs Pressure (Figure 3B) 
plot(Pres_mbar,Aspect,'k','LineWidth',2.5)


%% Functions

%%% Main Model function 
function [Pres_mbar, Xf_shell, X_rear, Width, Height, Aspect] = model(D0,E,tau0,th,d,alpha)
    
    alpha = alpha*pi/180;           % Conversion to radian
    A0 = pi*D0^2;                   % Initia area of the cell 
    V0 = pi*D0^3/6;                 % Initial volume of the cell 
    opts = optimset('TolX', 1e-20);

    %%%%%%%%%%%%%%%%%%%%%%% First Stage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    W = 300;                        % Length of the tapered section 
    X0 = 0;                         
    E_m = th * E;                   % Effective Elasticity (kPa.um) 
    Xfcap0 = X0 - W + D0;           % Initial Guess for Xfcap = Xf - hf 
    % We will finally delete the inadmissible values (negative pressure) in
    % the post processing stage 
    
    % First stage input (Xfcap = Xf - hf) 
    Xfcapvec = Xfcap0 : (X0-Xfcap0)/1000 : X0;      
    
    % Initializing the vectors of the model outputs 
    Dfvec1 = zeros(1,length(Xfcapvec));
    hfvec1 = zeros(1,length(Xfcapvec));
    Drvec1 = zeros(1,length(Xfcapvec));
    hrvec1 = zeros(1,length(Xfcapvec));
    Xrcapvec1 = zeros(1,length(Xfcapvec));
    Xfvec1 = zeros(1,length(Xfcapvec));
    Xrvec1 = zeros(1,length(Xfcapvec));
    Rfvec1 = zeros(1,length(Xfcapvec));
    Rrvec1 = zeros(1,length(Xfcapvec));
    widthvec1 = zeros(1,length(Xfcapvec));
    heightvec1 = zeros(1,length(Xfcapvec));
    Pvec1 = zeros(1,length(Xfcapvec));          % Pressure vector (kPa) 
    
    Xrcap0 = Xfcap0 - D0/2;                     % Initial Guess for the first use of fzero command 
    
    for ii = 1:length(Xfcapvec)
    
        Xfcap = Xfcapvec(ii);           % Input of stage 1 in each iteration 
        Df = -2*tan(alpha)*(Xfcap - X0) + d;
        hf = Df/2*(-tan(alpha) + (1 + tan(alpha)^2)^0.5);
        [Xrcap,fval] = fzero(@(Xrcap) firstphase(Xrcap,alpha,d,Df,hf,D0,X0),Xrcap0,opts);  % Incompressibilty eq
        Xrcap0 = Xrcap;                 % Updating the initial guess for the next iteration 
        Dr = -2*tan(alpha)*(Xrcap - X0) + d;
        hr = Dr/2*(tan(alpha) + (1 + tan(alpha)^2)^0.5);
        Xf = Xfcap + hf; 
        Xr = Xrcap - hr; 
    
        %%% The radii for using the formula 
        Rf = ((Df/2)^2+ hf^2)/(2*hf); 
        Rr = ((Dr/2)^2+ hr^2)/(2*hr); 
    
        %%% Areas 
        Ar = pi*(Dr^2/4 + hr^2); 
        Af = pi*(Df^2/4 + hf^2);
        Acone = pi/tan(alpha)*(Dr^2/4 - Df^2/4); 
        A = Ar + Af + Acone; 
        tau = tau0 + E_m*(A/A0 - 1); 
    
        %%% Outputs 
        Xrvec1(ii) = Xr; 
        Dfvec1(ii) = Df; 
        hfvec1(ii) = hf;
        Xrcapvec1(ii) = Xrcap; 
        Drvec1(ii) = Dr; 
        hrvec1(ii) = hr; 
        Rfvec1(ii) = Rf; 
        Rrvec1(ii) = Rr; 
        Xfvec1(ii) = Xf; 
        widthvec1(ii) = Xf - Xr;
        heightvec1(ii) = 2*Rr; 
    
        Pvec1(ii) = 2*tau*(1/Rf - 1/Rr); 
    end

    %%%%%%%%%%%%%% Second Stage (Prior to P_Xf=d/2) %%%%%%%%%%%%%%%%%%%%%%%
    
    Df = d;                                 % In this phase DR = d but hR is not related to DR through tangent relation 
    
    % Second stage input 
    hfvec2 = hf : (d/2-hf)/1000 : d/2;                    
    
    % Initializing the vectors of the model outputs 
    Drvec2 = zeros(1,length(hfvec2));
    hrvec2 = zeros(1,length(hfvec2));
    Xrcapvec2 = zeros(1,length(hfvec2));
    Xfvec2 = zeros(1,length(hfvec2));
    Xrvec2 = zeros(1,length(hfvec2));
    Rfvec2 = zeros(1,length(hfvec2));
    Rrvec2 = zeros(1,length(hfvec2));
    widthvec2 = zeros(1,length(hfvec2));
    heightvec2 = zeros(1,length(hfvec2));
    Pvec2 = zeros(1,length(hfvec2));              % Pressure vector (kPa) 
    
    Xrcap0 = Xrcap;                               % Initial Guess for fzero command
    
    for ii = 1:length(hfvec2)
        hf = hfvec2(ii); 
        Xrcap = fzero(@(Xrcap) secondphase(Xrcap,alpha,d,Df,hf,D0,X0),Xrcap0);      % Incompressibilty equation 
        Xrcap0 = Xrcap;                          % Updating the initial guess for the next iteration 
        Dr = -2*tan(alpha)*(Xrcap - X0) + d;
        hr = Dr/2*(tan(alpha) + (1 + tan(alpha)^2)^0.5);
        Xf = Xfcap + hf;
        Xr = Xrcap - hr; 
    
        %%% The radii for using the formula 
        Rf = ((Df/2)^2+ hf^2)/(2*hf); 
        Rr = ((Dr/2)^2+ hr^2)/(2*hr); 
    
        %%% Areas 
        Ar = pi*(Dr^2/4 + hr^2); 
        Af = pi*(Df^2/4 + hf^2);
        Acone = pi/tan(alpha)*(Dr^2/4 - Df^2/4); 
        A = Ar + Af + Acone; 
        tau = tau0 + E_m*(A/A0 - 1); 
    
    
        %%% Outpus 
        Xrvec2(ii) = Xr; 
        Xrcapvec2(ii) = Xrcap; 
        Drvec2(ii) = Dr; 
        hrvec2(ii) = hr; 
        Rfvec2(ii) = Rf; 
        Rrvec2(ii) = Rr; 
        Xfvec2(ii) = Xf;
        heightvec2(ii) = 2*Rr;
        widthvec2(ii) = Xf - Xr; 
    
        Pvec2(ii) = 2*tau*(1/Rf - 1/Rr);
    end

    %%%%%%%%%%%%%%%% Third Phase (Not formulated in the manuscript) %%%%%%%%% 
    % Since the maximum pressure in this stage corresponds to the model
    % maximum pressure plotted in Figure 3B (different from P_Xf=d/2) we also coded this stage  
    
    ymax = pi/6*(D0^3 - d^3)/(pi*d^2/4);            % Maximum height of the cylinder (when the partial conde disappears)          
   
    % Third stage input (The cylinder height starting from zero reaching to
    % its maximum value) 
    yvec = 0 : ymax/1000 : ymax; 
    
    Df = d; 
    hf = d/2;
    
    % Initializing the model output vector 
    Drvec3 = zeros(1,length(yvec));
    hrvec3 = zeros(1,length(yvec));
    Xrcapvec3 = zeros(1,length(yvec));
    Xfvec3 = zeros(1,length(yvec));
    Xrvec3 = zeros(1,length(yvec));
    Rfvec3 = zeros(1,length(yvec));
    Rrvec3 = zeros(1,length(yvec));
    heightvec3 = zeros(1,length(yvec));
    widthvec3 = zeros(1,length(yvec));
    
    Pvec3 = zeros(1,length(yvec));              % Pressure (kPa) 
    
    
    for ii = 1:length(yvec)
        y = yvec(ii); 
        Xrcap = fzero(@(Xrcap) thirdphase(Xrcap,alpha,d,Df,hf,D0,X0,y),Xrcap0);         % Incompressibilty Equation 
        Dr = -2*tan(alpha)*(Xrcap - X0) + d;
        hr = Dr/2*(tan(alpha) + (1 + tan(alpha)^2)^0.5);
        Xf = X0 + y + d/2; 
        Xr = Xrcap - hr; 
    
        %%% The radii for using the formula 
        Rf = d/2; 
        Rr = ((Dr/2)^2+ hr^2)/(2*hr); 
    
        Ar = pi*(Dr^2/4 + hr^2); 
        Af = pi*(Df^2/4 + hf^2);
        Acone = pi/tan(alpha)*(Dr^2/4 - Df^2/4); 
        Ay = pi*Df*y; 
        A = Ar + Af + Acone + Ay; 
        tau = tau0 + E_m*(A/A0 - 1); 
    
        Xrvec3(ii) = Xr; 
        Xfvec3(ii) = Xf; 
        Xrcapvec3(ii) = Xrcap; 
        Drvec3(ii) = Dr; 
        hrvec3(ii) = hr; 
        Rfvec3(ii) = Rf; 
        Rrvec3(ii) = Rr; 
        heightvec3(ii) = 2*Rr;
        widthvec3(ii) = Xf - Xr; 
    
        Pvec3(ii) = 2*tau*(1/Rf - 1/Rr);
    
    end

    %%%% Post Processing 
    Pressure_model = [Pvec1 Pvec2 Pvec3];           % Pressure vector (kPa)
    X_rear = [Xrvec1 Xrvec2 Xrvec3];                % Xr vector 
    Width = [widthvec1 widthvec2 widthvec3];        % Width vector 
    Height = 2*[Rrvec1 Rrvec2 Rrvec3];              % Height vector 
    Xf_shell = [Xfvec1 Xfvec2 Xfvec3];              % Xf vector 


    %%%% Finding the start of the experiments (Pressure >= 0)
    aval = 0; 
    [~,pzero] = find(Pressure_model > aval,1,'first'); 

    %%%% Finding the maximum pressure 
    [fmax,maxind] = max(Pressure_model); 

    %%%% Conversion from kPa to mbar and choosing the admissible values (0 < Pressure < P_max)  
    Pres_mbar = 10*Pressure_model(pzero : maxind);
    X_rear = X_rear(pzero : maxind); 
    Width = Width(pzero : maxind); 
    Height = Height(pzero : maxind); 
    Aspect = (Width)./(Height); 
    Xf_shell = Xf_shell(pzero : maxind);   
end

 
% First stage incompressibilty equation 
function F = firstphase(Xrcap,alpha,d,Df,hf,D0,X0)
  
    Vf = pi*hf*(3*Df^2/4 + hf^2)/6; 
    Dr = -2*tan(alpha)*(Xrcap - X0) + d;
    hr = Dr/2*(tan(alpha) + (1 + tan(alpha)^2)^0.5);
    Vr = pi*hr*(3*Dr^2/4 + hr^2)/6; 
    Vcone = (pi/24/tan(alpha))*(Dr^3 - Df^3);
    V0 = pi*D0^3/6; 
    F = Vf + Vr + Vcone - V0; 

end

% Second stage incompressibilty equation 
function F = secondphase(Xrcap,alpha,d,Df,hf,D0,X0)

    Vf = pi*hf/6*(3*Df^2/4 + hf^2); 
    Dr = -2*tan(alpha)*(Xrcap - X0) + d;
    hr = Dr/2*(tan(alpha) + (1 + tan(alpha)^2)^0.5);
    Vr = pi*hr/6*(3*Dr^2/4 + hr^2); 
    Vcone = (pi/24/tan(alpha))*(Dr^3 - Df^3);
    V0 = pi*D0^3/6; 
    F = Vf + Vr + Vcone - V0; 

end

% Third stage incompressibilty equation (Not formulated in the manuscript) 
function F = thirdphase(Xrcap,alpha,d,Df,hf,D0,X0,y)

    Vf = pi*hf/6*(3*Df^2/4 + hf^2); 
    Dr = -2*tan(alpha)*(Xrcap - X0) + d;
    hr = Dr/2*(tan(alpha) + (1 + tan(alpha)^2)^0.5);
    Vr = pi*hr/6*(3*Dr^2/4 + hr^2); 
    Vy = pi*d^2*y/4; 
    Vcone = (pi/24/tan(alpha))*(Dr^3 - Df^3);
    V0 = pi*D0^3/6; 
    F = Vf + Vr + Vcone + Vy - V0; 

end
