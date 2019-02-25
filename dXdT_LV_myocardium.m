function f = dXdT_LV_myocardium(~,X,P_in,P_LV,dPdT,fact1,fact2,fact3)

P_im1 = 1.2*0.167*P_LV;
P_im2 = 1.2*0.500*P_LV;
P_im3 = 1.2*0.833*P_LV;
dPim1_dt = 1.2*0.167*dPdT;
dPim2_dt = 1.2*0.500*dPdT;
dPim3_dt = 1.2*0.833*dPdT;
P_RA = 0; % right atrial pressure (mmHg)

% PARAMETERS
C_PA = 0.0013/3;  % mL / mmHg
L_PA = 2.0; % ?????
R_PA = 4; % mmHg / (mL / sec)
R_PV = 2; % mmHg / (mL / sec)
C_PV = 0.0254/3; % mL / mmHg
R0m = 44; % mmHg / (mL / sec)
R01 = 1.2*R0m;
R02 = 0.5*R0m;
V01 = 2.5/9; % mL
V02 = 8.0/9; % mL
C1 = 0.013/9; % mL / mmHg
C2 = 0.254/9; % mL / mmHg
gamma = 0.75; 
cf1 = 0.55; % epi/endo compliance factor
rf1 = 1.28; % epi/endo resistance factor
cf2 = 0.68; % epi/mid compliance factor
rf2 = 1.12; % epi/mid resistance factor

% STATE VARIABLES
P_PA = X(1); % penetrating artery pressure
Q_PA = X(2); % inlet flow penetrating artery
P11  = X(3); 
P21  = X(4);
P12  = X(5); 
P22  = X(6);
P13  = X(7); 
P23  = X(8);
P_PV = X(9); % penetrating vein pressure

% CALCULATIONS 
V11 = cf1*((P11 - P_im1)*fact1*C1+V01);
V21 = cf1*((P21 - P_im1)*C2+V02);
R11 = rf1*R01*(V01/V11)^2;
R21 = rf1*R02*(V02/V21)^2;
Rm1 = R0m*(gamma*R11/R01 + (1-gamma)*R21/R02);
Q11 = (P_PA - P11)/R11;
Qm1 = (P11 - P21)/Rm1;
Q21 = (P21 - P_PV)/R21;

V12 = cf2*((P12 - P_im2)*fact2*C1+V01);
V22 = cf2*((P22 - P_im2)*C2+V02);
R12 = rf2*R01*(V01/V12)^2;
R22 = rf2*R02*(V02/V22)^2;
Rm2 = R0m*(gamma*R12/R01 + (1-gamma)*R22/R02);
Q12 = (P_PA - P12)/R12;
Qm2 = (P12 - P22)/Rm2;
Q22 = (P22 - P_PV)/R22;

V13 = (P13 - P_im3)*fact3*C1+V01;
V23 = (P23 - P_im3)*C2+V02;
R13 = R01*(V01/V13)^2;
R23 = R02*(V02/V23)^2;
Rm3 = R0m*(gamma*R13/R01 + (1-gamma)*R23/R02);
Q13 = (P_PA - P13)/R13;
Qm3 = (P13 - P23)/Rm3;
Q23 = (P23 - P_PV)/R23;

Q_ima = Q11 + Q12 + Q13;
Q_imv = Q21 + Q22 + Q23;
Q_out = (P_PV - P_RA)/R_PV;

f(1,:) = (Q_PA - Q_ima)/C_PA; % P_PA
f(2,:) = (P_in - P_PA - Q_PA*R_PA)/L_PA; % Q_PA
f(3,:) = (Q11-Qm1)/(fact1*cf1*C1) + dPim1_dt; % P11
f(4,:) = (Qm1-Q21)/(cf1*C2) + dPim1_dt; % P21
f(5,:) = (Q12-Qm2)/(fact2*cf2*C1) + dPim2_dt; % P12
f(6,:) = (Qm2-Q22)/(cf2*C2) + dPim2_dt; % P22
f(7,:) = (Q13-Qm3)/(fact3*C1) + dPim3_dt; % P13
f(8,:) = (Qm3-Q23)/(C2) + dPim3_dt; % P23
f(9,:) = (Q_imv - Q_out)/C_PV; % P_PV

