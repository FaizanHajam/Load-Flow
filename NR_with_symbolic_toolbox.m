%%NR program in polar form using the symbolic toolbox
clc
clear
BMVA=input('Input the base MVA \n')
tic
%% A five-bus system with Q limits specified at the PV buses (We are considering the nominal-pi model)
%Line data of the system
%           FB TB  R      X    B/2
linedata = [1  2  0.042 0.168 0.041;
            1  5  0.031 0.126 0.031;
            2  3  0.031 0.126 0.031;
            3  4  0.024 0.136 0.082;
            3  5  0.053 0.210 0.051;
            4  5  0.063 0.252 0.061]; 
%%Bus data of the system

% Bus type 1 slack bus, 2 PV, 3 PQ       
%                 No.  Type |v|  delta Pg    Qg    Pl      Ql      Qmin Qmax
busdata=         [1       1       1    0     0     0       0  0    0    0;
                  2       2       1    0     50    0       0  0   -50  50;
                  3       2       1    0    100    0       0  0   -50  50;
                  4       3       1    0     0     0      115 60   0    0;
                  5       3       1    0     0     0       85 40   0    0]; 




%% Formulation of Y_bus and extracting the information
[nbus, nbranch, Ybus,  Yb, theta, G, B]=Y_bus(linedata);
[bus, type, V, del, P, Q, Qmin, Qmax]=busdatafn(busdata,BMVA);
npvs_initial=length(find(type==1|type==2)); %Determine the no. of PV buses including the slack bus
%% Define the symbolic variables
syms P_sym Q_sym [nbus,1]
syms PQ_gen [2*nbus,1]
syms del_sym V_sym [nbus,1] 
syms theta_sym Y_sym [nbus,nbus]
assume(Y_sym,'real')
assume(theta_sym,'real')
Psp=P;  %Specified power
Qsp=Q;
[P_sym, Q_sym]=power_symbolic_expression(nbus, V_sym, del_sym, Y_sym, theta_sym); % Symbolic expression for power at all buses
[J_gen,PQ_gen]=generalizedjacobianNR(P_sym, Q_sym,V_sym, del_sym); % Expression for the jacobian for all the buses, that is, from bus 1 to nbus; note that the jacobian is calculated with respect to all deltas and voltages
iter=1; %Set the iteration counter
error=1;  %Tolerance

while error>1e-10  % tolerance
%% Check for violation of Q limit
[Pcal,Qcal]=Powercalculation(P_sym, Q_sym, nbus, Yb, theta, V, del, V_sym, del_sym, Y_sym, theta_sym); % Numerical values of the power calculated at all buses
% Note: We check for the limit violation in each iteration. If in some
% iteration, Q at a PV bus is within the specified limits, switch the bus
% back to PV bus. 
for m=2:npvs_initial    
                 if Qcal(m)>Qmax(m) % if calculated power at the PV bus exceeds the maximum reactive power
                        Qsp(m)=Qmax(m); % then substitute specified reactive power at bus m as maximum reactive power that is already defined
                        type(m)=3; % moreover, set the bus type as PQ bus
             
                 elseif Qcal(m)<Qmin(m) % Vice versa
                        Qsp(m)=Qmin(m);
                        type(m)=3;  %In this case also, set type to PQ
                 else 
                     type(m)=2;
                 end
end

%% If voltage limit is violated at a PQ bus
% for m=npvs_initial+1:nbus
%     if V(m)>Vmax(m)
%         V(m)=Vmax(m);
%         type(m)=2;
%     elseif V(m)<Vmin(m)
%         V(m)=Vmin(m);
%         type(m)=2;
%     else 
%         type(m)=3;
%     end
% end
PVbus=find(type==1|type==2);  % All the PV buses, including the slack bus
PQbus= find(type==3); % PQ buses
npvs=length(PVbus);  % No. of slack plus PV buses
npq=length(PQbus);   % No. of PQ buses
%% Determine the jacobian
[PQ_cal, J_cal]=jacobianNR(J_gen, PQ_gen, nbus,npvs, Yb, theta, V, del, V_sym, del_sym, Y_sym, theta_sym); % Calculate Jacobian
dpa=Psp(2:nbus)-PQ_cal(1:nbus-1); % Delta P at 2:nbus
dqa=Qsp(npvs+1:nbus)-PQ_cal(nbus:end); % Delta Q at the all the PQ buses
mismatch1=[dpa;dqa]; %updated mismatches dP and dQ with slack removed
error=max(abs(mismatch1)); %updating the error 
if error < 1e-10
    break 
else
X1=J_cal\(mismatch1);
Ddel=X1(1:nbus-1);
DV=X1(nbus:end);
del(2:nbus)=del(2:nbus)+Ddel;
V(npvs+1:nbus)=V(npvs+1:nbus)+DV;
iter=iter+1;
end
end
%DISPLAYING RESULTS:
disp('RESULTS')
disp('The total no of iterations is')
iter
disp('Final V and del are:')
Voltages=V 
Delta=del
toc
%% Power expressions
function [P_sym, Q_sym]=power_symbolic_expression(nbus, V_sym, del_sym, Y_sym, theta_sym)
for i=1:nbus
P_sym(i)=sum(V_sym(i).*V_sym.*(Y_sym(i,:))'.*cos(del_sym(i)-del_sym-(theta_sym(i,:))'),1);
Q_sym(i)=sum(V_sym(i).*V_sym.*(Y_sym(i,:))'.*sin(del_sym(i)-del_sym-(theta_sym(i,:))'),1);
end
end
%% Power evaluation
function [Pcal,Qcal]=Powercalculation(P_sym, Q_sym, nbus, Yb, theta, V, del, V_sym, del_sym, Y_sym, theta_sym)
Pcal=double(subs(P_sym,[del_sym, V_sym, theta_sym,Y_sym],[del, V, theta, Yb]));
Qcal=double(subs(Q_sym,[del_sym, V_sym, theta_sym,Y_sym],[del, V, theta, Yb]));
end
%% Determine generalized Jacobian
function [J_gen,PQ_gen]=generalizedjacobianNR(P_sym, Q_sym,V_sym, del_sym)
PQ_gen=[P_sym'; Q_sym'];
J_gen=jacobian(PQ_gen,[del_sym;V_sym]);
end
%% Jacobian evaluation
function [PQ_cal, J_cal]=jacobianNR(J_gen, PQ_gen, nbus,npvs, Yb, theta, V, del, V_sym, del_sym, Y_sym, theta_sym)
PQ=PQ_gen([2:nbus,nbus+npvs+1:2*nbus]);
J=J_gen([2:nbus,nbus+npvs+1:2*nbus],[2:nbus,nbus+npvs+1:2*nbus]);
J_cal=double(subs(J,[del_sym, V_sym, theta_sym,Y_sym],[del, V, theta, Yb]));
PQ_cal=double(subs(PQ,[del_sym, V_sym, theta_sym,Y_sym],[del, V, theta, Yb]));
end