function f = dCdT_oxygen(t,x,T,V_pa,V_11,V_12,V_21,V_22,V_13, ...
    V_23,V_pv,Q_pa,Q_11,Q_m1,Q_12,Q_21,Q_m2,Q_22,Q_13,Q_m3,Q_23,Q_pv)

% volumes
v_pa = interp1(T,V_pa,t);
v_11 = interp1(T,V_11,t);
v_12 = interp1(T,V_12,t);
v_21 = interp1(T,V_21,t);
v_22 = interp1(T,V_22,t);
v_13 = interp1(T,V_13,t);
v_23 = interp1(T,V_23,t);
v_pv = interp1(T,V_pv,t) + 0.1; % add offset to v_pv to keep it stays positive

%flows
q_pa = interp1(T,Q_pa,t);
q_11 = interp1(T,Q_11,t);
q_m1 = interp1(T,Q_m1,t);
q_12 = interp1(T,Q_12,t);
q_21 = interp1(T,Q_21,t);
q_m2 = interp1(T,Q_m2,t);
q_22 = interp1(T,Q_22,t);
q_13 = interp1(T,Q_13,t);
q_m3 = interp1(T,Q_m3,t);
q_23 = interp1(T,Q_23,t);
q_pv = interp1(T,Q_pv,t);

% masses
m_pa = x(1);
m_11 = x(2);
m_21 = x(3);
m_12 = x(4);
m_22 = x(5);
m_13 = x(6);
m_23 = x(7);
m_pv = x(8);
c_t1 = x(9);
c_t2 = x(10);
c_t3 = x(11);

% concentrations
c_pa = m_pa/v_pa;
c_11 = m_11/v_11;
c_21 = m_21/v_21;
c_12 = m_12/v_12;
c_22 = m_22/v_22;
c_13 = m_13/v_13;
c_23 = m_23/v_23;
c_pv = m_pv/v_pv;

% parameters
CA = 1; % input concentration
MVO2 = 1.4449*0.9228/3;                 % Factor 1.00 Rest, 1.4449 Exercise
PS = 100;
v_t = 2.5;

f(1,:) = q_pa*CA*(q_pa>0) + q_pa*c_pa*(q_pa<0) - ...        % m_pa
    q_11*c_pa*(q_11>0) - q_11*c_11*(q_11<0) - ...
    q_12*c_pa*(q_12>0) - q_12*c_12*(q_12<0) - ...
    q_13*c_pa*(q_13>0) - q_13*c_13*(q_13<0); 
f(2,:) = q_11*c_pa*(q_11>0) + q_11*c_11*(q_11<0) - ...      % m_11
    q_m1*c_11*(q_m1>0) - q_m1*c_21*(q_m1<0); 
f(3,:) = q_m1*c_11*(q_m1>0) + q_m1*c_21*(q_m1<0) - ...      % m_21
    q_21*c_21*(q_21>0) - q_21*c_pv*(q_21<0) ...
    - PS*(c_21 - c_t1); 
f(4,:) = q_12*c_pa*(q_12>0) + q_12*c_12*(q_12<0) - ...      % m_12
    q_m2*c_12*(q_m2>0) - q_m2*c_22*(q_m2<0); 
f(5,:) = q_m2*c_12*(q_m2>0) + q_m2*c_22*(q_m2<0) - ...      % m_22
    q_22*c_22*(q_22>0) - q_22*c_pv*(q_22<0) ...
    - PS*(c_22 - c_t2); 
f(6,:) = q_13*c_pa*(q_13>0) + q_13*c_13*(q_13<0) - ...      % m_13
    q_m3*c_13*(q_m3>0) - q_m3*c_23*(q_m3<0); 
f(7,:) = q_m3*c_13*(q_m3>0) + q_m3*c_23*(q_m3<0) - ...      % m_23
    q_23*c_23*(q_23>0) - q_23*c_pv*(q_23<0) ...
    - PS*(c_23 - c_t3);
f(8,:) = q_21*c_21*(q_21>0) + q_21*c_pv*(q_21<0) + ...      % m_pv
         q_22*c_22*(q_22>0) + q_22*c_pv*(q_22<0) + ...
         q_23*c_23*(q_23>0) + q_23*c_pv*(q_23<0) - q_pv*c_pv; 
f(9,:)  = + (PS/v_t)*(c_21 - c_t1) - MVO2/v_t;
f(10,:) = + (PS/v_t)*(c_22 - c_t2) - MVO2/v_t;
f(11,:) = + (PS/v_t)*(c_23 - c_t3) - MVO2/v_t;

