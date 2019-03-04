%% Inputting flow and pressure data
% clear;
tic

data = xlsread('TuneExercisePig','2713 Resting','B9:D5005');
dt    = 1/500;
Pao   = data(:,2);
Plv   = 0.85*(data(:,3)-17); % adjusting PLV data to match aortic pressure systole
Plv   = smoothdata(Plv,'gaussian','smoothingfactor',0.015); %smoothing makes the numerics easier
Fdata = data(:,1);
tdata = (0:(length(Pao)-1)).*dt;

figure(1); clf; axes('position',[0.15 0.15 0.75 0.75]); hold on;
plot(tdata,Plv,'k-','linewidth',3,'color',0.5*[1 1 1]);
plot(tdata,Pao,'k-','linewidth',3);
set(gca,'Fontsize',18); box on
xlabel('time (sec)','interpreter','latex','fontsize',20);
ylabel('Pressure (mmHg)','interpreter','latex','fontsize',20);
axis([5 10 0 180]);
legend('LV','Aorta');

figure(2); clf; axes('position',[0.15 0.15 0.75 0.75]); hold on;
plot(tdata,Plv,'k-','linewidth',2,'color',0.5*[1 1 1]);
plot(tdata,Pao,'k-','linewidth',2);
set(gca,'Fontsize',16); box on
xlabel('time (sec)','interpreter','latex','fontsize',18);
ylabel('Pressure (mmHg)','interpreter','latex','fontsize',18);
axis([0 10 -2 180]);
legend('LV','Aorta');
%%
% 2-point derivative of Plv
dPlvdt(1,:) = (Plv(2)-Plv(1))./dt;
for i = 2:(length(Plv)-1)
  dPlvdt(i,:) = (Plv(i+1)-Plv(i-1))./(2*dt);
end
dPlvdt(length(Plv),:) = (Plv(length(Plv))-Plv(length(Plv)-1))./dt;

Pdata = [tdata; Pao'; Plv'; dPlvdt'];
% figure(5); plot(tdata,dPlvdt);

%%
% Setting up network 
segID       = xlsread('LAD Network.xlsx','Coronary 1D Model Parameters','A2:A36');
[~,segNAME] = xlsread('LAD Network.xlsx','Coronary 1D Model Parameters','B2:B36');
len      = xlsread('LAD Network.xlsx','Coronary 1D Model Parameters','C2:C36');
area        = xlsread('LAD Network.xlsx','Coronary 1D Model Parameters','D2:D36');
Co          = xlsread('LAD Network.xlsx','Coronary 1D Model Parameters','E2:E36');
[~,segP]    = xlsread('LAD Network.xlsx','Coronary 1D Model Parameters','F2:F36');
[~,segD1]   = xlsread('LAD Network.xlsx','Coronary 1D Model Parameters','G2:G36');
[~,segD2]   = xlsread('LAD Network.xlsx','Coronary 1D Model Parameters','H2:H36');

% Segment connectivity
Nseg = 35;
d1 = zeros(Nseg,1);
d2 = zeros(Nseg,1);
PARENT = zeros(Nseg,1);
myocardial_outlets = [];
for i = 1:Nseg
  v = char(segD1{i});
  if strcmp(v(1),'P') % the outlet is a perfusion zone
    d1(i) = 999;
    myocardial_outlets = [myocardial_outlets; i];
  end
  v = char(segD2{i});
  if strcmp(v(1),'P')||strcmp(v,'LVfw')||strcmp(v,'Sep') % the outlet is a perfusion zone
    d2(i) = 999;
    myocardial_outlets = [myocardial_outlets; i];
  end
  for j = 1:Nseg
    if strcmp(segD1{i},segNAME{j})
      d1(i) = j;
    end
    if strcmp(segD2{i},segNAME{j})
      d2(i) = j;
    end
    
    if strcmp(segP{i},segNAME{j})
      PARENT(i) = j;
    end
  end
end

fact1 = 1.0;
fact2 = 1.0;
fact3 = 1.0;

% segment parameters
a = 0.2802;
b = -0.5053; % per mm
c = 0.1325;
d = -0.01114; % per mm
E = 10000000; % value ???
mu = (3e-3)/133.3; % value mmHg sec
rho = 1060; % kg / M^3
f = 0.01; % sec
r = sqrt(area/pi); % radii in mm
h = r.*( a*exp(b*r) + c*exp(d*r) ); % wall thickness in mm
C = (2*pi/E)*(r.^3).*len./h; % compliance in mL/mmHg
R = (8*mu/pi)*len./(r.^4)*1000;
Rv = f./C;
L = rho*len./(pi*r.^2)*1000; % (in kg / M^4)
L = L./(133.3e6); % convert to s^2 mmHg / mL

%% Initial Conditions
Po = 100*ones(Nseg,1);
Qo = zeros(Nseg,1);
Xo = [Po; Qo];
Xo_myo = [60 1 50 50 85 85 120 120 5]'; % for 2713 Resting
for i = 1:28
  Xo = [Xo; Xo_myo];
end

%% Simulation

% Xo = X(end,:)';
[t,X] = ode23s(@dXdT_LAD,[0 9.99],Xo,[],d1,d2,PARENT,myocardial_outlets,Rv,R,C,L,Pdata,fact1,fact2,fact3);

Pv = X(:,1:Nseg);
Q  = X(:,(Nseg+1):(2*Nseg)); 
Qout = zeros(size(Q)); 
P = zeros(size(Q));
myocardial_i = 0;
for i = 1:Nseg
   if d1(i) ~= 999
     Qout(:,i) = Qout(:,i) + Q(:,d1(i));
   elseif d2(i) ~= 999
     myocardial_i = myocardial_i + 1;
     Qout(:,i) = Qout(:,i) + X(:,70 + 2 + 9*(myocardial_i-1) )/35;
   else
     myocardial_i = myocardial_i + 1;
     Qout(:,i) = Qout(:,i) + 2*X(:,70 + 2 + 9*(myocardial_i-1) )/35;
   end

   if (d2(i) > 0) && (d2(i)~=999)
     Qout(:,i) = Qout(:,i) + Q(:,d2(i));
   end
   P(:,i) = Pv(:,i) + (Q(:,i)-Qout(:,i)).*Rv(i); 
end  
% Q(:,i) is the flow into the ith epicarial segment
% Qout(i,:) is the flow out of the ith epicardial segment
% X(:,70 + 2 + 9*(myocardial_i-1) )/35 is the flow into the ith myocardial outlet

% figure(2); clf; axes('position',[0.15 0.15 0.75 0.50]); hold on;
% set(gca,'Fontsize',14); box on
% plot(t,P,'linewidth',1.5,'color','k');
% set(gca,'fontsize',14);
% ylabel('Epicardial pressures (mmHg)');
% xlabel('time (sec)');

figure(3); clf; axes('position',[0.15 0.15 0.75 0.75]);  hold on;
plot(tdata,Fdata,'color',0.5*[1 1 1],'linewidth',3);
% plot(t,60*Xout(:,2),'color','r','linewidth',1.5); 
plot(t,60*Q(:,1),'color','r','linewidth',1.5); 
l = legend('data','model'); 
set(l,'fontsize',12,'location','northeast');
set(gca,'fontsize',14); box on;
ylabel('Myocardial flow (mL/min)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);
axis([5 10 -25 150]); grid

figure(4); clf; axes('position',[0.15 0.15 0.75 0.75]);  hold on;
plot(tdata,Fdata,'color',0.5*[1 1 1],'linewidth',3);
%l = legend('data','model'); 
%set(l,'fontsize',12,'location','northeast');
set(gca,'fontsize',18); box on;
ylabel('Myocardial flow (mL/min)','interpreter','latex','fontsize',20);
xlabel('time (sec)','interpreter','latex','fontsize',20);
axis([5 10 -25 160]); grid

% figure(4); clf; axes('position',[0.15 0.15 0.75 0.75]); hold on;
% plot(t,60*Xout(:,2),'color','r','linewidth',1.5); 
% plot(t,60*Xout(:,9)./2,'color','b','linewidth',1.5); 
% l = legend('LAD','venous'); 
% set(l,'fontsize',12,'location','southeast');
% set(gca,'fontsize',14); box on;
% xlabel('time (sec)','interpreter','latex','fontsize',16);
% ylabel('Myocardial flow (mL/min)','interpreter','latex','fontsize',16);
% % axis([0 10 -50 200])

%% Myocardium

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
P_RA = 0; % right atrial pressure (mmHg)

P_LV = interp1(Pdata(1,:),Pdata(3,:),t);
P_im1 = 1.2*0.167*P_LV;
P_im2 = 1.2*0.500*P_LV;
P_im3 = 1.2*0.833*P_LV;

% Flows to myocardium
Q13_list = [];
Q12_list = [];
Q11_list = [];
Qm3_list = [];
Qm2_list = [];
Qm1_list = [];
Q23_list = [];
Q22_list = [];
Q21_list = [];
Qima_list = [];
Qimv_list = [];
Qpa_list = [];
Qpv_list = [];

% Myocardial vascular volumes
V13_list = [];
V12_list = [];
V11_list = [];
V23_list = [];
V22_list = [];
V21_list = [];
Vpa_list = [];
Vpv_list = [];

myocardial_i = 0;
for i = 1:Nseg
  n_out = 0;
  if (d1(i) == 999) && (d2(i) ~= 999)
    myocardial_i = myocardial_i + 1;
    n_out = 1; % number of outlets represented by myocardial_i
  elseif (d1(i) == 999) && (d2(i) == 999)
    myocardial_i = myocardial_i + 1;
    n_out = 2;
  end
  
  if n_out > 0
    out_i = 70 + 9*(myocardial_i-1) + (1:9); % indices for the myocardial_ith outlet model
    Xout = X(:,out_i);

    % STATE VARIABLES
    P_PA = Xout(:,1); % penetrating artery pressure
    Q_PA = Xout(:,2); % inlet flow penetrating artery
    P11  = Xout(:,3); 
    P21  = Xout(:,4);
    P12  = Xout(:,5); 
    P22  = Xout(:,6);
    P13  = Xout(:,7); 
    P23  = Xout(:,8);
    P_PV = Xout(:,9); % penetrating vein pressure

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
    
    Q_out = (P_PV - P_RA)./R_PV;
    
    Vpa = P_PA.*C_PA;
    Vpv = P_PV.*C_PV;

    % Assembling time series of flow and volume to all perfusion zones. 
    % (We add one column for a single outlet, two columns for an outlet
    % representing a two-outlet bifurcation.)
    Q13_list = [Q13_list, Q13*ones(1,n_out)];
    Q12_list = [Q12_list, Q12*ones(1,n_out)];
    Q11_list = [Q11_list, Q11*ones(1,n_out)];
    Qm3_list = [Qm3_list, Qm3*ones(1,n_out)];
    Qm2_list = [Qm2_list, Qm2*ones(1,n_out)];
    Qm1_list = [Qm1_list, Qm1*ones(1,n_out)];
    Q23_list = [Q23_list, Q23*ones(1,n_out)];
    Q22_list = [Q22_list, Q22*ones(1,n_out)];
    Q21_list = [Q21_list, Q21*ones(1,n_out)];
    Qpa_list = [Qpa_list, Q_PA*ones(1,n_out)];
    Qpv_list = [Qpv_list, Q_out*ones(1,n_out)];
    
    V13_list = [V13_list, V13*ones(1,n_out)];
    V12_list = [V12_list, V12*ones(1,n_out)];
    V11_list = [V11_list, V11*ones(1,n_out)];
    V23_list = [V23_list, V23*ones(1,n_out)];
    V22_list = [V22_list, V22*ones(1,n_out)];
    V21_list = [V21_list, V21*ones(1,n_out)];
    Vpa_list = [Vpa_list, Vpa*ones(1,n_out)];
    Vpv_list = [Vpv_list, Vpv*ones(1,n_out)];
    
  end
    
end

% Compute average epi, mid, and endo flow values
Q_13 = mean(Q13_list')';
Q_12 = mean(Q12_list')';
Q_11 = mean(Q11_list')';
Q_m3 = mean(Qm3_list')';
Q_m2 = mean(Qm2_list')';
Q_m1 = mean(Qm1_list')';
Q_23 = mean(Q23_list')';
Q_22 = mean(Q22_list')';
Q_21 = mean(Q21_list')';
Q_pa = mean(Qpa_list')';
Q_pv = mean(Qpv_list')';

% Compute average volumes
V_13 = mean(V13_list')';
V_12 = mean(V12_list')';
V_11 = mean(V11_list')';
V_23 = mean(V23_list')';
V_22 = mean(V22_list')';
V_21 = mean(V21_list')';
V_pv = mean(Vpv_list')';
V_pa = mean(Vpa_list')';

[mean(interp1(t,Q_11,4:0.01:9.99)) , mean(interp1(t,Q_12,4:0.01:9.99)), mean(interp1(t,Q_13,4:0.01:9.99))]

% Endo/Epi flow ratio:
EndoEpiRatio = mean(interp1(t,Q_13,4:0.01:9.99)) / mean(interp1(t,Q_11,4:0.01:9.99))

% Checking conservation of flow. Mean LAD flow should be equal to mean of endo+mid_epi
mean(interp1(t,Q(:,1),4:0.01:9.99))
mean(interp1(t,Q_11,4:0.01:9.99)) + mean(interp1(t,Q_12,4:0.01:9.99)) + mean(interp1(t,Q_13,4:0.01:9.99))

figure(5); clf; axes('position',[0.15 0.15 0.75 0.75]);  hold on;
plot(t,60*Q_13,'color','r','linewidth',1.5); 
plot(t,60*Q_12,'color','g','linewidth',1.5); 
plot(t,60*Q_11,'color','b','linewidth',1.5); 
l = legend('endo','mid','epi'); 
set(l,'fontsize',12,'location','northeast');
set(gca,'fontsize',14); box on;
ylabel('Myocardial flow (mL/min)','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);
axis([5 10 -40 100]); grid

toc

