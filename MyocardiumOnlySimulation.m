%% Inputting flow and pressure data
% clear;
% data = xlsread('TuneExercisePig','2713 Resting','B9:D5005');
data = xlsread('TuneExercisePig','2713 Exercise Level 2','B9:D5005');
dt    = 1/500;
Pao   = data(:,2);
Plv   = 0.85*(data(:,3)-17);
Plv   = smoothdata(Plv,'gaussian','smoothingfactor',0.015); %smoothing makes the numerics easier
Fdata = data(:,1);
tdata = (0:(length(Pao)-1)).*dt;

% figure(1); clf; axes('position',[0.15 0.15 0.75 0.75]); hold on;
% plot(tdata,Plv,'k-','linewidth',1.5,'color',0.5*[1 1 1]);
% plot(tdata,Pao,'k-','linewidth',1.5);
% set(gca,'Fontsize',14); box on
% xlabel('time (sec)','interpreter','latex','fontsize',16);
% ylabel('Pressure (mmHg)','interpreter','latex','fontsize',16);

% 2-point derivative of Plv
dPlvdt(1,:) = (Plv(2)-Plv(1))./dt;
for i = 2:(length(Plv)-1)
  dPlvdt(i,:) = (Plv(i+1)-Plv(i-1))./(2*dt);
end
dPlvdt(length(Plv),:) = (Plv(length(Plv))-Plv(length(Plv)-1))./dt;

Pdata = [tdata; Pao'; Plv'; dPlvdt'];

%% Initial Conditions and Simulation

% fact1 = 2.35;
% fact2 = 2.35;
% fact3 = 2.35;
fact1 = 2.05;
fact2 = 2.35;
fact3 = 3.20;

% Xo_myo = [60 1 50 50 85 85 120 120 5]'; % for 2713 Resting
Xo_myo = [60 1 10 10 10 10 10 10 5]'; % 2713 Exercise Level 2
[t,X] = ode23s(@dXdT_myocardium,[0 9.99],Xo_myo,[],Pdata,fact1,fact2,fact3);

% PARAMETERS
C_PA = 0.0013/3;  % mL / mmHg
L_PA = 1.0; % ?????
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

P_LV = interp1(Pdata(1,:),Pdata(3,:),t);
P_im1 = 1.2*0.167*P_LV;
P_im2 = 1.2*0.500*P_LV;
P_im3 = 1.2*0.833*P_LV;

% STATE VARIABLES
P_PA = X(:,1); % penetrating artery pressure
Q_PA = X(:,2); % inlet flow penetrating artery
P11  = X(:,3); 
P21  = X(:,4);
P12  = X(:,5); 
P22  = X(:,6);
P13  = X(:,7); 
P23  = X(:,8);
P_PV = X(:,9); % penetrating vein pressure

% CALCULATIONS 
V11 = cf1*((P11 - P_im1)*fact1*C1+V01);
V21 = cf1*((P21 - P_im1)*C2+V02);
R11 = rf1*R01*(V01./V11).^2;
R21 = rf1*R02*(V02./V21).^2;
Rm1 = R0m*(gamma*R11/R01 + (1-gamma)*R21/R02);
Q11 = (P_PA - P11)./R11;
Qm1 = (P11 - P21)./Rm1;
Q21 = (P21 - P_PV)./R21;

V12 = cf2*((P12 - P_im2)*fact2*C1+V01);
V22 = cf2*((P22 - P_im2)*C2+V02);
R12 = rf2*R01*(V01./V12).^2;
R22 = rf2*R02*(V02./V22).^2;
Rm2 = R0m*(gamma*R12/R01 + (1-gamma)*R22/R02);
Q12 = (P_PA - P12)./R12;
Qm2 = (P12 - P22)./Rm2;
Q22 = (P22 - P_PV)./R22;

V13 = (P13 - P_im3)*fact3*C1+V01;
V23 = (P23 - P_im3)*C2+V02;
R13 = R01*(V01./V13).^2;
R23 = R02*(V02./V23).^2;
Rm3 = R0m*(gamma*R13/R01 + (1-gamma)*R23/R02);
Q13 = (P_PA - P13)./R13;
Qm3 = (P13 - P23)./Rm3;
Q23 = (P23 - P_PV)./R23;

figure(3); clf; axes('position',[0.15 0.15 0.75 0.75]);  hold on;
plot(tdata,Fdata,'color',0.5*[1 1 1],'linewidth',3);
plot(t,60*X(:,2),'color','r','linewidth',1.5); 
% plot(t,60*(Q11+Q12+Q13),'color','b','linewidth',1.5); 
l = legend('data','model'); 
set(l,'fontsize',12,'location','southeast');
set(gca,'fontsize',14); box on;
ylabel('Myocardial flow (mL/min)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);
axis([5 10 -50 300])

% figure(5); plot(t,V11+V21,t,V12+V22,t,V13+V23)
% figure(6); plot(t,Q11,t,Q12,t,Q13); grid
% axis([5 10 -0.5 2])

[mean(interp1(t,Q11,4:0.01:9.99)), ...
mean(interp1(t,Q12,4:0.01:9.99)), ...
mean(interp1(t,Q13,4:0.01:9.99))]


