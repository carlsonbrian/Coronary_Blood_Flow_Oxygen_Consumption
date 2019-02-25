function f = dXdT_LAD(t,X,d1,d2,PARENT,myocardial_outlets,Rv,R,C,L,Pdata,alpha,beta1,beta2)

Nseg = 35; % number of epicarial segments

% INPUT Pressures
P_Ao = interp1(Pdata(1,:),Pdata(2,:),t);
P_LV = interp1(Pdata(1,:),Pdata(3,:),t);
dPlvdt = interp1(Pdata(1,:),Pdata(4,:),t);

% STATE VARIABLES
%  index 1-35: the Voigt body pressure in the segments
%  index 36-70: the flow (Q) in the segment
%  index 71-322: the 28 9-variable 0-d networks
Pv = X(1:Nseg);
Q  = X((Nseg+1):(2*Nseg));

Qout = zeros(Nseg,1);
myocardial_i = 0;
for i = 1:Nseg
   if d1(i) ~= 999 % no myocardial outlet
     Qout(i) = Qout(i) + Q(d1(i));
   elseif d2(i) ~= 999 % in this case, 1 myocardial outlet
     myocardial_i = myocardial_i + 1;
     Qout(i) = Qout(i) + X(70 + 2 + 9*(myocardial_i-1) )/35;
   else
     myocardial_i = myocardial_i + 1; % two myocardial outlets
     Qout(i) = Qout(i) + 2*X(70 + 2 + 9*(myocardial_i-1) )/35;
   end

   if (d2(i) > 0) && (d2(i)~=999)
     Qout(i) = Qout(i) + Q(d2(i));
   end
end  

P = Pv + (Q-Qout).*Rv; 
Pin(1,:) = P_Ao;
Pin(2:Nseg,:) = P(PARENT(2:Nseg));

f(1:Nseg,:) = (Q-Qout)./C;
f((Nseg+1):(2*Nseg),:) = (Pin - P - (Q.*R))./L;

for i = 1:28
  f( 70 + 9*(i-1) + (1:9), :) = dXdT_LV_myocardium(t,X(70 + 9*(i-1) + (1:9)),P(myocardial_outlets(i)),P_LV,dPlvdt,alpha,beta1,beta2);
end

