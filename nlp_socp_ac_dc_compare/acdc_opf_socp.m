clc; clear; close all;

%% S1.1 Load AC network parameters & variables 

% -------------------------------------------------------------------------
% 1.1.1 AC grid parameters
% -------------------------------------------------------------------------
mpc_ac = case9;                       
bus_ac     = mpc_ac.bus;
gen_ac     = mpc_ac.gen;
branch_ac  = mpc_ac.branch;
baseMVA_ac = mpc_ac.baseMVA;
gencost_ac = mpc_ac.gencost;

[Ybus_ac, ~, ~] = makeYbus(baseMVA_ac, bus_ac, branch_ac);
Ybus_ac = full(Ybus_ac);

G_ac = real(Ybus_ac);
B_ac = imag(Ybus_ac);

nbus_ac    = size(bus_ac, 1);
ngen_ac    = size(gen_ac, 1);
nbranch_ac = size(branch_ac, 1);

pd_ac = bus_ac(:, 3) / baseMVA_ac;
qd_ac = bus_ac(:, 4) / baseMVA_ac;

refBus_ac = find(bus_ac(:, 2) == 3);    
genBus_ac = gen_ac(:, 1);  
a2vBus_ac = [4; 6; 8];

% -------------------------------------------------------------------------
% 1.1.2 AC grid variables
% -------------------------------------------------------------------------
pn_ac       = sdpvar(nbus_ac, 1);
qn_ac       = sdpvar(nbus_ac, 1);
pgen_ac     = sdpvar(nbus_ac, 1);
qgen_ac     = sdpvar(nbus_ac, 1);
pa2v_ac     = sdpvar(nbus_ac, 1);
qa2v_ac     = sdpvar(nbus_ac, 1);
pl_ac       = sdpvar(nbus_ac, nbus_ac, 'full');
ql_ac       = sdpvar(nbus_ac, nbus_ac, 'full');
cc_ac       = sdpvar(nbus_ac, nbus_ac, 'full');
ss_ac       = sdpvar(nbus_ac, nbus_ac, 'full');


%% S1.2 Add SOCP-relaxed AC network constraints

% "Strong SOCP relaxations for the optimal power flow problem". 
% Operations Research, 2016, 64(6): 1177-1196.

con_ac = [];

% -------------------------------------------------------------------------
% 1.2.1 Branch active/reactive power flow equations
%   Smax: apparent power limits derived from MATPOWER branch ratings (rateA)
%   Polygonal approximation is used to enforce 
%   p_ij^2 + q_ij^2 <= Smax(i,j)^2
%
%   "Linear decision rules for control of reactive power by 
%   distributed photovoltaic generators[J]. "
%   IEEE Transactions on Power Systems, 2017, 33(2): 2165-2174."
% -------------------------------------------------------------------------
Smax_ac = zeros(nbus_ac, nbus_ac);

for k = 1:nbranch_ac
    i = branch_ac(k, 1);
    j = branch_ac(k, 2);

    rateA = branch_ac(k, 6);

    if rateA == 0
        Smax_val = 1e3;   
    else
        Smax_val = rateA / baseMVA_ac;
    end

    Smax_ac(i, j) = Smax_val;
    Smax_ac(j, i) = Smax_val;   
end

N = 16;
for k = 1:nbranch_ac
    i = branch_ac(k, 1);
    j = branch_ac(k, 2);

    con_ac = [con_ac, pl_ac(i, j) == -G_ac(i, j)*(cc_ac(i, i) ...
        -cc_ac(i, j))-B_ac(i, j)*ss_ac(i, j)];
    con_ac = [con_ac, ql_ac(i, j) == B_ac(i, j)*(cc_ac(i, i) ...
        -cc_ac(i, j))-G_ac(i, j)*ss_ac(i, j)];
    % ij and ji
    con_ac = [con_ac, pl_ac(j, i) == -G_ac(j, i)*(cc_ac(j, j) ...
        -cc_ac(j, i))-B_ac(j, i)*ss_ac(j, i)];
    con_ac = [con_ac, ql_ac(j, i) == B_ac(j, i)*(cc_ac(j, j) ...
        -cc_ac(j, i))-G_ac(j, i)*ss_ac(j, i)];
    
    for k = 1: N
        con_ac = [con_ac, cos(k/N*pi)*pl_ac(i, j) ...
            + sin(k/N*pi)*ql_ac(i, j) <= Smax_ac(i, j)];
        con_ac = [con_ac; -(cos(k/N*pi)*pl_ac(i, j) ...
            + sin(k/N*pi)*ql_ac(i, j)) <= Smax_ac(i, j)];
        % ij and ji
        con_ac = [con_ac, cos(k/N*pi)*pl_ac(j, i) ...
            + sin(k/N*pi)*ql_ac(j, i) <= Smax_ac(j, i)];
        con_ac = [con_ac; -(cos(k/N*pi)*pl_ac(j, i) ...
            + sin(k/N*pi)*ql_ac(j, i)) <= Smax_ac(j, i)];
    end
    
end

% -------------------------------------------------------------------------
% 1.2.2 Nodal active/reactive power flow equations
% -------------------------------------------------------------------------
for i = 1:nbus_ac
    con_ac = [con_ac, pn_ac(i) == ss_ac(i, i)*B_ac(i, i) ...
        + sum(cc_ac(i, :).*G_ac(i, :) - ss_ac(i, :).*B_ac(i, :))];
    con_ac = [con_ac, qn_ac(i) == ss_ac(i, i)*G_ac(i, i) ...
        - sum(cc_ac(i, :).*B_ac(i, :) + ss_ac(i, :).*G_ac(i, :))];
end

% -------------------------------------------------------------------------
% 1.2.3 Nodal power balance
% -------------------------------------------------------------------------
for i = 1:nbus_ac
    % if ~ismember(i, refBus_ac)
        p_expr = -pd_ac(i);
        q_expr = -qd_ac(i);

        if ismember(i, genBus_ac)
            p_expr = p_expr + pgen_ac(i);
            q_expr = q_expr + qgen_ac(i);
        end

        if ismember(i, a2vBus_ac)
            p_expr = p_expr - pa2v_ac(i);
            q_expr = q_expr - qa2v_ac(i);
        end

        con_ac = [con_ac, pn_ac(i) == p_expr];
        con_ac = [con_ac, qn_ac(i) == q_expr];

    % end
end

% -------------------------------------------------------------------------
% 1.2.4 SOCP relaxation constraints
%     Symmetry constraints
% -------------------------------------------------------------------------
for k = 1:nbranch_ac
    i = branch_ac(k, 1);
    j = branch_ac(k, 2);
    con_ac = [con_ac, cc_ac(i, j)^2 + ss_ac(i, j)^2 ...
        <= cc_ac(i, i)*cc_ac(j, j)];
end

for i = 1:nbus_ac
    con_ac = [con_ac, ss_ac(i, i) == 0];  
    con_ac=[con_ac, cc_ac(i, i)>=0];
    for j = i + 1:nbus_ac
        con_ac = [con_ac, cc_ac(i, j) == cc_ac(j, i)];
        con_ac = [con_ac, ss_ac(i, j) + ss_ac(j, i) == 0];
    end
end


% -------------------------------------------------------------------------
% 1.2.5 Generator active/reactive power limits
% -------------------------------------------------------------------------
pgenMin_ac = zeros(nbus_ac, 1);
pgenMax_ac = zeros(nbus_ac, 1);
qgenMin_ac = zeros(nbus_ac, 1);
qgenMax_ac = zeros(nbus_ac, 1);

for i = 1:nbus_ac
    gen_idx = find(gen_ac(:,1) == i, 1, 'first');
    if ~isempty(gen_idx)
        qgenMax_ac(i) = gen_ac(gen_idx, 4)  / baseMVA_ac;   
        qgenMin_ac(i) = gen_ac(gen_idx, 5)  / baseMVA_ac;   
        pgenMax_ac(i) = gen_ac(gen_idx, 9)  / baseMVA_ac;   
        pgenMin_ac(i) = gen_ac(gen_idx, 10) / baseMVA_ac;   
    end
end

for i = 1:nbus_ac
    if ismember(i, genBus_ac)
        con_ac = [con_ac, pgenMin_ac(i) <= pgen_ac(i) <= pgenMax_ac(i)];
        con_ac = [con_ac, qgenMin_ac(i) <= qgen_ac(i) <= qgenMax_ac(i)];
    else
        con_ac = [con_ac, pgen_ac(i) == 0];
        con_ac = [con_ac, qgen_ac(i) == 0];
    end
end

% -------------------------------------------------------------------------
% 1.2.6 Voltage magnitude bounds
% -------------------------------------------------------------------------
for i = 1:nbus_ac
    Vmax_ac = bus_ac(i, 12);
    Vmin_ac = bus_ac(i, 13);
    con_ac = [con_ac, Vmin_ac^2 <= cc_ac(i, i) <= Vmax_ac^2];
end

% -------------------------------------------------------------------------
% 1.2.7 Voltage references
% -------------------------------------------------------------------------
r = refBus_ac.'
gen_idx = find(gen_ac(:,1) == r, 1, 'first');
Vref_ac = gen_ac(gen_idx,6);
con_ac = [con_ac, cc_ac(r, r) == Vref_ac^2];

% -------------------------------------------------------------------------
% 1.2.8 Obj: Minizing generation costs
% -------------------------------------------------------------------------

c2_ac = zeros(nbus_ac, 1);
c1_ac = zeros(nbus_ac, 1);
c0_ac = zeros(nbus_ac, 1);

for i = 1:nbus_ac
    gen_idx = find(gen_ac(:,1) == i, 1, 'first');
    if ~isempty(gen_idx)
        ncost = gencost_ac(gen_idx, 4);

        if ncost == 3
            c2_ac(i) = gencost_ac(gen_idx, end-2);
            c1_ac(i) = gencost_ac(gen_idx, end-1);
            c0_ac(i) = gencost_ac(gen_idx, end);
        elseif ncost == 2
            c2_ac(i) = 0;
            c1_ac(i) = gencost_ac(gen_idx, end-1);
            c0_ac(i) = gencost_ac(gen_idx, end);
        end
    end
end

obj_ac = 0;

for i = 1:nbus_ac
    if ismember(i, genBus_ac)
        Pg_MW = pgen_ac(i) * baseMVA_ac;
        obj_ac = obj_ac + c2_ac(i) * Pg_MW^2 + c1_ac(i) * Pg_MW + c0_ac(i);
    end
end

% -------------------------------------------------------------------------
% 1.2.9 Solve
% -------------------------------------------------------------------------
% opt = sdpsettings('solver','gurobi', ...
%                   'verbose',1, ...
%                   'gurobi.QCPDual',1);
% 
% sol_ac = optimize(con_ac, obj_ac, opt);
% 
% if sol_ac.problem ~= 0
%     disp('Optimization failed.');
%     disp(sol_ac.info);
%     yalmiperror(sol_ac.problem);
%     return;
% end
% 
% disp('AC Optimization succeeded.');


%% S2.1 Load DC network parameters & variables 

% -------------------------------------------------------------------------
% 2.1.1 DC grid parameters
% -------------------------------------------------------------------------
mpc_dc      = case5_dc_5branch;                       
results_dc  = runpf(mpc_dc);
bus_dc      = results_dc.bus;
gen_dc      = results_dc.gen;
branch_dc   = results_dc.branch;
baseMW_dc   = results_dc.baseMVA;

[Ybus_dc, ~, ~] = makeYbus(baseMW_dc, bus_dc, branch_dc);
Ybus_dc = full(Ybus_dc);

nbus_dc         = size(bus_dc, 1);
ngen_dc         = size(gen_dc, 1);
nbranch_dc      = size(branch_dc, 1);

zij_dc = 1 ./ abs(Ybus_dc);
zij_dc = zij_dc - diag(diag(zij_dc));
zij_dc(isinf(zij_dc)) = 1e5;

bus_dc(:, 3) = bus_dc(:, 3) / baseMW_dc; 
pd_dc = bus_dc(:, 3);

fbus_dc = branch_dc(:, 1);
tbus_dc = branch_dc(:, 2);

refBus_dc = find(bus_dc(:, 2) == 3);     

% -------------------------------------------------------------------------
% 2.1.2 DC grid parameters
% -------------------------------------------------------------------------
u_dc        = sdpvar(nbus_dc, 1);              
lij_dc      = sdpvar(nbus_dc, nbus_dc, 'full');   
pn_dc       = sdpvar(nbus_dc, 1);              
pij_dc      = sdpvar(nbus_dc, nbus_dc, 'full');   

% switching transmission
alpha_dc    = binvar(nbus_dc, nbus_dc, 'full'); 
beta_dc     = binvar(nbus_dc, nbus_dc, 'full');
ep_dc       = binvar(nbus_dc, nbus_dc, 'full');
bij_dc      = sdpvar(nbus_dc, nbus_dc, 'full');
tij_dc      = sdpvar(nbus_dc, nbus_dc, 'full');


%% S2.2 Add SOCP-relaxed DC network constraints

% "Optimalpowerflowindirectcurrentnetworks," 
% IEEE Transactions on Power Systems,vol.29,no.6,pp.2892–2904, 2014.

con_dc=[]; 
M = 1e5;

% -------------------------------------------------------------------------
% 2.2.1 Nodal power balance
% -------------------------------------------------------------------------
for i = 1:nbus_dc
    expr = 0;
    
    for k = 1:nbranch_dc
        if fbus_dc(k) == i
            j = tbus_dc(k);
            expr = expr + pij_dc(i, j);
        elseif tbus_dc(k) == i
            j = fbus_dc(k);
            expr = expr + pij_dc(i, j);
        end
    end
    
    con_dc = [con_dc; pn_dc(i) == expr];
end

% -------------------------------------------------------------------------
% 2.2.2 Branch equations
% -------------------------------------------------------------------------
for k = 1:nbranch_dc
    i = fbus_dc(k);
    j = tbus_dc(k);
    rateA = branch_dc(k, 6);

    if rateA == 0
        Pmax_val = 1e3;   
    else
        Pmax_val = rateA / baseMW_dc;
    end
    Pmax_val = 100 / baseMW_dc;

    % loss
    con_dc = [con_dc;
        pij_dc(i, j) + pij_dc(j, i) == zij_dc(i, j) * lij_dc(i, j)];

    % volt drop
    con_dc = [con_dc;
       u_dc(i) - u_dc(j) == zij_dc(i, j) * (pij_dc(i, j) - pij_dc(j, i))];

    % SOCP
    con_dc = [con_dc;
        pij_dc(i, j)^2 <= lij_dc(i, j) * u_dc(i);
        pij_dc(j, i)^2 <= lij_dc(j, i) * u_dc(j)];

    con_dc = [con_dc;
        lij_dc(i, j) >= 0
        lij_dc(j, i) >= 0];

    con_dc = [con_dc; 
        u_dc(i) >= 0];

    con_dc = [con_dc;
        lij_dc(j, i) == lij_dc(i, j)];

    % Bounds
    con_dc = [con_dc;
       pij_dc(i, j) <= Pmax_val
       pij_dc(j, i) <= Pmax_val];

end

% ---------------------------------------------------------------------
% 2.2.3 Power Injection
% -------------------------------------------------------------------------
% for i = 1:nbus_dc
%     if ~ismember(i, refBus_dc)
%         con_dc = [con_dc; pn_dc(i) == -pd_dc(i)];
%     end
% end

% -------------------------------------------------------------------------
% 2.2.4 Slack bus voltage
% -------------------------------------------------------------------------
% r = refBus_dc.'
% gen_idx = find(gen_dc(:,1) == r, 1, 'first');
% Vref_dc = gen_dc(gen_idx,6);
% con_dc = [con_dc, u_dc(r) == Vref_dc^2];

% -------------------------------------------------------------------------
% 2.2.5 Voltage magnitude bounds
% -------------------------------------------------------------------------
for i = 1:nbus_dc
    Vmax_dc = bus_dc(i,12);
    Vmin_dc = bus_dc(i,13);
    con_dc = [con_dc, Vmin_dc^2 <= u_dc(i) <= Vmax_dc^2];
end

% -------------------------------------------------------------------------
% 2.2.8 Solve
% -------------------------------------------------------------------------
% opt = sdpsettings('solver','gurobi','verbose',1);
% sol_dc = optimize(con_dc, obj_dc, opt);
% 
% if sol_dc.problem ~= 0
%     disp('Optimization failed.');
%     disp(sol_dc.info);
%     return;
% end
% 
% disp('DC Optimization succeeded.');



%% S3.1 Load RES parameters & variables 

% -------------------------------------------------------------------------
% 3.1.1 RES parameters 
% -------------------------------------------------------------------------
baseMVA_res1 = 100;
poutMax_res1 = 80 / baseMVA_res1;
SoutMax_res1 = 100 / baseMVA_res1;


baseMVA_res2 = 100;
poutMax_res2 = 70 / baseMVA_res2;
SoutMax_res2 = 100 / baseMVA_res2;

% -------------------------------------------------------------------------
% 3.1.2 RES variables
% -------------------------------------------------------------------------
pout_res1    = sdpvar(1, 1);
qout_res1    = sdpvar(1, 1);

pr2v_res1    = sdpvar(1, 1);
qr2v_res1    = sdpvar(1, 1);


pout_res2    = sdpvar(1, 1);
qout_res2    = sdpvar(1, 1);

pr2v_res2    = sdpvar(1, 1);
qr2v_res2    = sdpvar(1, 1);


%% S3.2 RES constraints 
con_res1 = [];
con_res2 = [];

% -------------------------------------------------------------------------
% 3.1.1 Capacity constraints
%   Polygonal approximation is used to enforce 
%   p_out^2 + q_out^2 <= SoutMax(i,j)^2
% -------------------------------------------------------------------------
N = 16;

for k = 1: N
        con_res1 = [con_res1; cos(k/N*pi)*pout_res1 ...
            + sin(k/N*pi)*qout_res1 <= SoutMax_res1];
        con_res1 = [con_res1; -(cos(k/N*pi)*pout_res1 ...
            + sin(k/N*pi)*qout_res1) <= SoutMax_res1];
end


for k = 1: N
        con_res2 = [con_res2; cos(k/N*pi)*pout_res2 ...
            + sin(k/N*pi)*qout_res2 <= SoutMax_res2];
        con_res2 = [con_res2; -(cos(k/N*pi)*pout_res2 ...
            + sin(k/N*pi)*qout_res2) <= SoutMax_res2];
end

% -------------------------------------------------------------------------
% 3.1.2 Maximum power bound
% -------------------------------------------------------------------------
con_res1 = [con_res1; 0 <= pout_res1 <= poutMax_res1];

con_res2 = [con_res2; 0 <= pout_res2 <= poutMax_res2];


% -------------------------------------------------------------------------
% 3.1.3 Speicify power
% -------------------------------------------------------------------------

con_res1 = [con_res1; pout_res1 == pr2v_res1;
            qout_res1 == qr2v_res1];

con_res2 = [con_res2; pout_res2 == pr2v_res2;
            qout_res2 == qr2v_res2];


% -------------------------------------------------------------------------
% 3.1.3 Speicify power
% -------------------------------------------------------------------------
% obj_res1 = sum(pout_res1);
% obj_res2 = sum(pout_res2);

% -------------------------------------------------------------------------
% 3.1.4 Solve
% -------------------------------------------------------------------------
% opt = sdpsettings('solver','gurobi','verbose',1);
% sol_res = optimize([con_res1; con_res2], obj_res1+obj_res2, opt);
% 
% if sol_res.problem ~= 0
%     disp('Optimization failed.');
%     disp(sol_res.info);
%     return;
% end
% 
% disp('RES Optimization succeeded.');


%% S4.1 Load VSC network parameters & variables 

% -------------------------------------------------------------------------
% 4.1.1 VSC parameters
% -------------------------------------------------------------------------
mpc_vsc1 = case2_vsc;                       
bus_vsc1     = mpc_vsc1.bus;
branch_vsc1  = mpc_vsc1.branch;
baseMVA_vsc1 = mpc_vsc1.baseMVA;

[Ybus_vsc1, ~, ~] = makeYbus(baseMVA_vsc1, bus_vsc1, branch_vsc1);
Ybus_vsc1 = full(Ybus_vsc1);

G_vsc1 = real(Ybus_vsc1);
B_vsc1 = imag(Ybus_vsc1);

nbus_vsc1    = size(bus_vsc1, 1);
nbranch_vsc1 = size(branch_vsc1, 1);

pd_vsc1 = bus_vsc1(:, 3) / baseMVA_vsc1;
qd_vsc1 = bus_vsc1(:, 4) / baseMVA_vsc1;

refBus_vsc1 = find(bus_vsc1(:, 2) == 3);    


mpc_vsc2 = case2_vsc;                       
bus_vsc2     = mpc_vsc2.bus;
branch_vsc2  = mpc_vsc2.branch;
baseMVA_vsc2 = mpc_vsc2.baseMVA;

[Ybus_vsc2, ~, ~] = makeYbus(baseMVA_vsc2, bus_vsc2, branch_vsc2);
Ybus_vsc2 = full(Ybus_vsc2);

G_vsc2 = real(Ybus_vsc2);
B_vsc2 = imag(Ybus_vsc2);

nbus_vsc2    = size(bus_vsc2, 1);
nbranch_vsc2 = size(branch_vsc2, 1);

pd_vsc2 = bus_vsc2(:, 3) / baseMVA_vsc2;
qd_vsc2 = bus_vsc2(:, 4) / baseMVA_vsc2;

refBus_vsc2 = find(bus_vsc2(:, 2) == 3);   


mpc_vsc3 = case2_vsc;                       
bus_vsc3     = mpc_vsc3.bus;
branch_vsc3  = mpc_vsc3.branch;
baseMVA_vsc3 = mpc_vsc3.baseMVA;

[Ybus_vsc3, ~, ~] = makeYbus(baseMVA_vsc3, bus_vsc3, branch_vsc3);
Ybus_vsc3 = full(Ybus_vsc3);

G_vsc3 = real(Ybus_vsc3);
B_vsc3 = imag(Ybus_vsc3);

nbus_vsc3    = size(bus_vsc3, 1);
nbranch_vsc3 = size(branch_vsc3, 1);

pd_vsc3 = bus_vsc2(:, 3) / baseMVA_vsc3;
qd_vsc3 = bus_vsc2(:, 4) / baseMVA_vsc3;

refBus_vsc3 = find(bus_vsc3(:, 2) == 3);


mpc_vsc4 = case2_vsc;                       
bus_vsc4     = mpc_vsc4.bus;
branch_vsc4  = mpc_vsc4.branch;
baseMVA_vsc4 = mpc_vsc4.baseMVA;

[Ybus_vsc4, ~, ~] = makeYbus(baseMVA_vsc4, bus_vsc4, branch_vsc4);
Ybus_vsc4 = full(Ybus_vsc4);

G_vsc4 = real(Ybus_vsc4);
B_vsc4 = imag(Ybus_vsc4);

nbus_vsc4    = size(bus_vsc4, 1);
nbranch_vsc4 = size(branch_vsc4, 1);

pd_vsc4 = bus_vsc2(:, 3) / baseMVA_vsc4;
qd_vsc4 = bus_vsc2(:, 4) / baseMVA_vsc4;

refBus_vsc4 = find(bus_vsc4(:, 2) == 3);


mpc_vsc5 = case2_vsc;                       
bus_vsc5     = mpc_vsc5.bus;
branch_vsc5  = mpc_vsc5.branch;
baseMVA_vsc5 = mpc_vsc5.baseMVA;

[Ybus_vsc5, ~, ~] = makeYbus(baseMVA_vsc5, bus_vsc5, branch_vsc5);
Ybus_vsc5 = full(Ybus_vsc5);

G_vsc5 = real(Ybus_vsc5);
B_vsc5 = imag(Ybus_vsc5);

nbus_vsc5    = size(bus_vsc5, 1);
nbranch_vsc5 = size(branch_vsc5, 1);

pd_vsc5 = bus_vsc2(:, 3) / baseMVA_vsc5;
qd_vsc5 = bus_vsc2(:, 4) / baseMVA_vsc5;

refBus_vsc5 = find(bus_vsc5(:, 2) == 3);

aloss_vsc=0.0011; bloss_vsc=0.003; closs_vsc=0.0046;
Imax_vsc = 1;
sMax_vsc1 = 100 /baseMVA_vsc1 ;
sMax_vsc2 = 100 /baseMVA_vsc2 ;
sMax_vsc3 = 100 /baseMVA_vsc3 ;
sMax_vsc4 = 100 /baseMVA_vsc4 ;
sMax_vsc5 = 100 /baseMVA_vsc5 ;

% -------------------------------------------------------------------------
% 4.1.2 VSC variables
% -------------------------------------------------------------------------
pn_vsc1       = sdpvar(nbus_vsc1, 1);
qn_vsc1       = sdpvar(nbus_vsc1, 1);
pgen_vsc1     = sdpvar(nbus_vsc1, 1);
qgen_vsc1     = sdpvar(nbus_vsc1, 1);
pa2v_vsc1     = sdpvar(nbus_vsc1, 1);
qa2v_vsc1     = sdpvar(nbus_vsc1, 1);
pl_vsc1       = sdpvar(nbus_vsc1, nbus_vsc1, 'full');
ql_vsc1       = sdpvar(nbus_vsc1, nbus_vsc1, 'full');
cc_vsc1       = sdpvar(nbus_vsc1, nbus_vsc1, 'full');
ss_vsc1       = sdpvar(nbus_vsc1, nbus_vsc1, 'full');
pv2d_vsc1     = sdpvar(1, 1);
pv2a_vsc1     = sdpvar(1, 1);
qv2a_vsc1     = sdpvar(1, 1);
I_vsc1        = sdpvar(1, 1);
l_vsc1        = sdpvar(1, 1);
ploss_vsc1    = sdpvar(1, 1);
vmin_vsc1     = 0.9;
vmax_vsc1     = 1.1;

pn_vsc2       = sdpvar(nbus_vsc2, 1);
qn_vsc2       = sdpvar(nbus_vsc2, 1);
pgen_vsc2     = sdpvar(nbus_vsc2, 1);
qgen_vsc2     = sdpvar(nbus_vsc2, 1);
pa2v_vsc2     = sdpvar(nbus_vsc2, 1);
qa2v_vsc2     = sdpvar(nbus_vsc2, 1);
pl_vsc2       = sdpvar(nbus_vsc2, nbus_vsc2, 'full');
ql_vsc2       = sdpvar(nbus_vsc2, nbus_vsc2, 'full');
cc_vsc2       = sdpvar(nbus_vsc2, nbus_vsc2, 'full');
ss_vsc2       = sdpvar(nbus_vsc2, nbus_vsc2, 'full');
pv2d_vsc2     = sdpvar(1, 1);
pv2a_vsc2     = sdpvar(1, 1);
qv2a_vsc2     = sdpvar(1, 1);
I_vsc2        = sdpvar(1, 1);
l_vsc2        = sdpvar(1, 1);
ploss_vsc2    = sdpvar(1, 1);
vmin_vsc2     = 0.9;
vmax_vsc2     = 1.1;

pn_vsc3       = sdpvar(nbus_vsc3, 1);
qn_vsc3       = sdpvar(nbus_vsc3, 1);
pgen_vsc3     = sdpvar(nbus_vsc3, 1);
qgen_vsc3     = sdpvar(nbus_vsc3, 1);
pa2v_vsc3     = sdpvar(nbus_vsc3, 1);
qa2v_vsc3     = sdpvar(nbus_vsc3, 1);
pl_vsc3       = sdpvar(nbus_vsc3, nbus_vsc3, 'full');
ql_vsc3       = sdpvar(nbus_vsc3, nbus_vsc3, 'full');
cc_vsc3       = sdpvar(nbus_vsc3, nbus_vsc3, 'full');
ss_vsc3       = sdpvar(nbus_vsc3, nbus_vsc3, 'full');
pv2d_vsc3     = sdpvar(1, 1);
pv2a_vsc3     = sdpvar(1, 1);
qv2a_vsc3     = sdpvar(1, 1);
I_vsc3        = sdpvar(1, 1);
l_vsc3        = sdpvar(1, 1);
ploss_vsc3    = sdpvar(1, 1);
vmin_vsc3     = 0.9;
vmax_vsc3     = 1.1;

pn_vsc4       = sdpvar(nbus_vsc4, 1);
qn_vsc4       = sdpvar(nbus_vsc4, 1);
pgen_vsc4     = sdpvar(nbus_vsc4, 1);
qgen_vsc4     = sdpvar(nbus_vsc4, 1);
pa2v_vsc4     = sdpvar(nbus_vsc4, 1);
qa2v_vsc4     = sdpvar(nbus_vsc4, 1);
pl_vsc4       = sdpvar(nbus_vsc4, nbus_vsc4, 'full');
ql_vsc4       = sdpvar(nbus_vsc4, nbus_vsc4, 'full');
cc_vsc4       = sdpvar(nbus_vsc4, nbus_vsc4, 'full');
ss_vsc4       = sdpvar(nbus_vsc4, nbus_vsc4, 'full');
pv2d_vsc4     = sdpvar(1, 1);
pv2a_vsc4     = sdpvar(1, 1);
qv2a_vsc4     = sdpvar(1, 1);
I_vsc4        = sdpvar(1, 1);
l_vsc4        = sdpvar(1, 1);
ploss_vsc4    = sdpvar(1, 1);
vmin_vsc4     = 0.9;
vmax_vsc4     = 1.1;

pn_vsc5       = sdpvar(nbus_vsc5, 1);
qn_vsc5       = sdpvar(nbus_vsc5, 1);
pgen_vsc5     = sdpvar(nbus_vsc5, 1);
qgen_vsc5     = sdpvar(nbus_vsc5, 1);
pa2v_vsc5     = sdpvar(nbus_vsc5, 1);
qa2v_vsc5     = sdpvar(nbus_vsc5, 1);
pl_vsc5       = sdpvar(nbus_vsc5, nbus_vsc5, 'full');
ql_vsc5       = sdpvar(nbus_vsc5, nbus_vsc5, 'full');
cc_vsc5       = sdpvar(nbus_vsc5, nbus_vsc5, 'full');
ss_vsc5       = sdpvar(nbus_vsc5, nbus_vsc5, 'full');
pv2d_vsc5     = sdpvar(1, 1);
pv2a_vsc5     = sdpvar(1, 1);
qv2a_vsc5     = sdpvar(1, 1);
I_vsc5        = sdpvar(1, 1);
l_vsc5        = sdpvar(1, 1);
ploss_vsc5    = sdpvar(1, 1);
vmin_vsc5     = 0.9;
vmax_vsc5     = 1.1;

%% S4.2 VSC constraints

con_vsc1 = [];
con_vsc2 = [];
con_vsc3 = [];
con_vsc4 = [];
con_vsc5 = [];

% -------------------------------------------------------------------------
% 4.2.1 Nodal active/reactive power flow equations
% -------------------------------------------------------------------------
for i = 1:nbus_vsc1
    con_vsc1 = [con_vsc1; pn_vsc1(i) == ss_vsc1(i, i)*B_vsc1(i, i) ...
        + sum(cc_vsc1(i, :).*G_vsc1(i, :) - ss_vsc1(i, :).*B_vsc1(i, :))];
    con_vsc1 = [con_vsc1; qn_vsc1(i) == ss_vsc1(i, i)*G_vsc1(i, i) ...
        - sum(cc_vsc1(i, :).*B_vsc1(i, :) + ss_vsc1(i, :).*G_vsc1(i, :))];
    con_vsc1 = [con_vsc1; vmin_vsc1^2^2 <= cc_vsc1(i, i) <= vmax_vsc1^2];
end

for i = 1:nbus_vsc2
    con_vsc2 = [con_vsc2; pn_vsc2(i) == ss_vsc2(i, i)*B_vsc2(i, i) ...
        + sum(cc_vsc2(i, :).*G_vsc2(i, :) - ss_vsc2(i, :).*B_vsc2(i, :))];
    con_vsc2 = [con_vsc2; qn_vsc2(i) == ss_vsc2(i, i)*G_vsc2(i, i) ...
        - sum(cc_vsc2(i, :).*B_vsc2(i, :) + ss_vsc2(i, :).*G_vsc2(i, :))];
    con_vsc2 = [con_vsc2; vmin_vsc2^2 <= cc_vsc2(i, i) <= vmax_vsc2^2];
end

for i = 1:nbus_vsc3
    con_vsc3 = [con_vsc3; pn_vsc3(i) == ss_vsc3(i, i)*B_vsc3(i, i) ...
        + sum(cc_vsc3(i, :).*G_vsc3(i, :) - ss_vsc3(i, :).*B_vsc3(i, :))];
    con_vsc3 = [con_vsc3; qn_vsc3(i) == ss_vsc3(i, i)*G_vsc3(i, i) ...
        - sum(cc_vsc3(i, :).*B_vsc3(i, :) + ss_vsc3(i, :).*G_vsc3(i, :))];
    con_vsc3 = [con_vsc3; vmin_vsc3^2 <= cc_vsc3(i, i) <= vmax_vsc3^2];
end

for i = 1:nbus_vsc4
    con_vsc4 = [con_vsc4; pn_vsc4(i) == ss_vsc4(i, i)*B_vsc4(i, i) ...
        + sum(cc_vsc4(i, :).*G_vsc4(i, :) - ss_vsc4(i, :).*B_vsc4(i, :))];
    con_vsc4 = [con_vsc4; qn_vsc4(i) == ss_vsc4(i, i)*G_vsc4(i, i) ...
        - sum(cc_vsc4(i, :).*B_vsc4(i, :) + ss_vsc4(i, :).*G_vsc4(i, :))];
    con_vsc4 = [con_vsc4; vmin_vsc4^2 <= cc_vsc4(i, i) <= vmax_vsc4^2];
end

for i = 1:nbus_vsc5
    con_vsc5 = [con_vsc5; pn_vsc5(i) == ss_vsc5(i, i)*B_vsc5(i, i) ...
        + sum(cc_vsc5(i, :).*G_vsc5(i, :) - ss_vsc5(i, :).*B_vsc5(i, :))];
    con_vsc5 = [con_vsc5; qn_vsc5(i) == ss_vsc5(i, i)*G_vsc5(i, i) ...
        - sum(cc_vsc5(i, :).*B_vsc5(i, :) + ss_vsc5(i, :).*G_vsc5(i, :))];
    con_vsc5 = [con_vsc5; vmin_vsc5^2 <= cc_vsc5(i, i) <= vmax_vsc5^2];
end


% -------------------------------------------------------------------------
% 4.2.2 PCC side power balance
% -------------------------------------------------------------------------

con_vsc1 = [con_vsc1; pn_vsc1(1) == -pv2a_vsc1];
con_vsc1 = [con_vsc1; qn_vsc1(1) == -qv2a_vsc1];    

con_vsc2 = [con_vsc2; pn_vsc2(1) == -pv2a_vsc2];
con_vsc2 = [con_vsc2; qn_vsc2(1) == -qv2a_vsc2];    

con_vsc3 = [con_vsc3; pn_vsc3(1) == -pv2a_vsc3];
con_vsc3 = [con_vsc3; qn_vsc3(1) == -qv2a_vsc3];    

con_vsc4 = [con_vsc4; pn_vsc4(1) == -pv2a_vsc4];
con_vsc4 = [con_vsc4; qn_vsc4(1) == -qv2a_vsc4];    

con_vsc5 = [con_vsc5; pn_vsc5(1) == -pv2a_vsc5];
con_vsc5 = [con_vsc5; qn_vsc5(1) == -qv2a_vsc5];    

% -------------------------------------------------------------------------
% 4.2.3 Capacity cosntraint
%   Polygonal approximation is used to enforce 
%   pn_vsc^2 + qn_vsc^2 <= SMax_vsc^2
% -------------------------------------------------------------------------

N = 16;

for k = 1: N
        con_vsc1 = [con_vsc1; cos(k/N*pi)*pn_vsc1(2) ...
            + sin(k/N*pi)*qn_vsc1(2) <= sMax_vsc1];
        con_vsc1 = [con_vsc1; -(cos(k/N*pi)*pn_vsc1(2) ...
            + sin(k/N*pi)*qn_vsc1(2)) <= sMax_vsc1];

        con_vsc2 = [con_vsc2; cos(k/N*pi)*pn_vsc2(2) ...
            + sin(k/N*pi)*qn_vsc2(2) <= sMax_vsc2];
        con_vsc2 = [con_vsc2; -(cos(k/N*pi)*pn_vsc2(2) ...
            + sin(k/N*pi)*qn_vsc2(2)) <= sMax_vsc2];


        con_vsc3 = [con_vsc3; cos(k/N*pi)*pn_vsc3(2) ...
            + sin(k/N*pi)*qn_vsc3(2) <= sMax_vsc3];
        con_vsc3 = [con_vsc3; -(cos(k/N*pi)*pn_vsc3(2) ...
            + sin(k/N*pi)*qn_vsc3(2)) <= sMax_vsc3];

        con_vsc4 = [con_vsc4; cos(k/N*pi)*pn_vsc4(2) ...
            + sin(k/N*pi)*qn_vsc4(2) <= sMax_vsc4];
        con_vsc4 = [con_vsc4; -(cos(k/N*pi)*pn_vsc4(2) ...
            + sin(k/N*pi)*qn_vsc4(2)) <= sMax_vsc4];

        con_vsc5 = [con_vsc5; cos(k/N*pi)*pn_vsc5(2) ...
            + sin(k/N*pi)*qn_vsc5(2) <= sMax_vsc5];
        con_vsc5 = [con_vsc5; -(cos(k/N*pi)*pn_vsc5(2) ...
            + sin(k/N*pi)*qn_vsc5(2)) <= sMax_vsc5];
end


% -------------------------------------------------------------------------
% 4.2.4 VSC power loss
% -------------------------------------------------------------------------
con_vsc1 = [con_vsc1; pn_vsc1(2) + ploss_vsc1 + pv2d_vsc1 == 0;
    ploss_vsc1 == aloss_vsc + bloss_vsc*I_vsc1 +closs_vsc*l_vsc1;
    0 <= I_vsc1 <= Imax_vsc;
    0 <= l_vsc1 <= Imax_vsc^2
    ];

con_vsc2 = [con_vsc2; pn_vsc2(2) + ploss_vsc2 + pv2d_vsc2 == 0;
    ploss_vsc2 == aloss_vsc + bloss_vsc*I_vsc2 +closs_vsc*l_vsc2;
    0 <= I_vsc2 <= Imax_vsc;
    0 <= l_vsc2 <= Imax_vsc^2
 ];

con_vsc3 = [con_vsc3; pn_vsc3(2) + ploss_vsc3 + pv2d_vsc3 == 0;
    ploss_vsc3 == aloss_vsc + bloss_vsc*I_vsc3 +closs_vsc*l_vsc3;
    0 <= I_vsc3 <= Imax_vsc;
    0 <= l_vsc3 <= Imax_vsc^2
 ];
  
con_vsc4 = [con_vsc4; pn_vsc4(2) + ploss_vsc4 + pv2d_vsc4 == 0;
    ploss_vsc4 == aloss_vsc + bloss_vsc*I_vsc4 +closs_vsc*l_vsc4;
    0 <= I_vsc4 <= Imax_vsc;
    0 <= l_vsc4 <= Imax_vsc^2
 ];

con_vsc5 = [con_vsc5; pn_vsc5(2) + ploss_vsc5 + pv2d_vsc5 == 0;
    ploss_vsc5 == aloss_vsc + bloss_vsc*I_vsc5 +closs_vsc*l_vsc5;
    0 <= I_vsc5 <= Imax_vsc;
    0 <= l_vsc5 <= Imax_vsc^2
 ];


% -------------------------------------------------------------------------
% 4.2.5 SOCP relaxation
% -------------------------------------------------------------------------
con_vsc1 = [con_vsc1; I_vsc1^2 <= l_vsc1;
   pn_vsc1(2)^2 + qn_vsc1(2)^2 <= cc_vsc1(2, 2) * l_vsc1];

for k = 1:nbranch_vsc1
    i = branch_vsc1(k, 1);
    j = branch_vsc1(k, 2);
    con_vsc1 = [con_vsc1, cc_vsc1(i, j)^2 + ss_vsc1(i, j)^2 ...
        <= cc_vsc1(i, i)*cc_vsc1(j, j)];
end

for i = 1:nbus_vsc1
    con_vsc1 = [con_vsc1, ss_vsc1(i, i) == 0];  
    con_vsc1=[con_vsc1, cc_vsc1(i, i)>=0];
    for j = i + 1:nbus_vsc1
        con_vsc1 = [con_vsc1, cc_vsc1(i, j) == cc_vsc1(j, i)];
        con_vsc1 = [con_vsc1, ss_vsc1(i, j) + ss_vsc1(j, i) == 0];
    end
end

con_vsc2 = [con_vsc2; I_vsc2^2 <= l_vsc2;
   pn_vsc2(2)^2 + qn_vsc2(2)^2 <= cc_vsc2(2, 2) * l_vsc2];

for k = 1:nbranch_vsc2
    i = branch_vsc2(k, 1);
    j = branch_vsc2(k, 2);
    con_vsc2 = [con_vsc2, cc_vsc2(i, j)^2 + ss_vsc2(i, j)^2 ...
        <= cc_vsc2(i, i)*cc_vsc2(j, j)];
end

for i = 1:nbus_vsc2
    con_vsc2 = [con_vsc2, ss_vsc2(i, i) == 0];  
    con_vsc2=[con_vsc2, cc_vsc2(i, i)>=0];
    for j = i + 1:nbus_vsc2
        con_vsc2 = [con_vsc2, cc_vsc2(i, j) == cc_vsc2(j, i)];
        con_vsc2 = [con_vsc2, ss_vsc2(i, j) + ss_vsc2(j, i) == 0];
    end
end

con_vsc3 = [con_vsc3; I_vsc3^2 <= l_vsc3;
   pn_vsc3(2)^2 + qn_vsc3(2)^2 <= cc_vsc3(2, 2) * l_vsc3];

for k = 1:nbranch_vsc3
    i = branch_vsc3(k, 1);
    j = branch_vsc3(k, 2);
    con_vsc3 = [con_vsc3, cc_vsc3(i, j)^2 + ss_vsc3(i, j)^2 ...
        <= cc_vsc3(i, i)*cc_vsc3(j, j)];
end

for i = 1:nbus_vsc3
    con_vsc3 = [con_vsc3, ss_vsc3(i, i) == 0];  
    con_vsc3=[con_vsc3, cc_vsc3(i, i)>=0];
    for j = i + 1:nbus_vsc3
        con_vsc3 = [con_vsc3, cc_vsc3(i, j) == cc_vsc3(j, i)];
        con_vsc3 = [con_vsc3, ss_vsc3(i, j) + ss_vsc3(j, i) == 0];
    end
end

con_vsc4 = [con_vsc4; I_vsc4^2 <= l_vsc4;
   pn_vsc4(2)^2 + qn_vsc4(2)^2 <= cc_vsc4(2, 2) * l_vsc4];

for k = 1:nbranch_vsc4
    i = branch_vsc4(k, 1);
    j = branch_vsc4(k, 2);
    con_vsc4 = [con_vsc4, cc_vsc4(i, j)^2 + ss_vsc4(i, j)^2 ...
        <= cc_vsc4(i, i)*cc_vsc4(j, j)];
end

for i = 1:nbus_vsc4
    con_vsc4 = [con_vsc4, ss_vsc4(i, i) == 0];  
    con_vsc4 =[con_vsc4, cc_vsc4(i, i)>=0];
    for j = i + 1:nbus_vsc4
        con_vsc4 = [con_vsc4, cc_vsc4(i, j) == cc_vsc4(j, i)];
        con_vsc4 = [con_vsc4, ss_vsc4(i, j) + ss_vsc4(j, i) == 0];
    end
end

con_vsc5 = [con_vsc5; I_vsc5^2 <= l_vsc5;
   pn_vsc5(2)^2 + qn_vsc5(2)^2 <= cc_vsc5(2, 2) * l_vsc5];

for k = 1:nbranch_vsc5
    i = branch_vsc5(k, 1);
    j = branch_vsc5(k, 2);
    con_vsc5 = [con_vsc5, cc_vsc5(i, j)^2 + ss_vsc5(i, j)^2 ...
        <= cc_vsc5(i, i)*cc_vsc5(j, j)];
end

for i = 1:nbus_vsc5
    con_vsc5 = [con_vsc5, ss_vsc5(i, i) == 0];  
    con_vsc5=[con_vsc5, cc_vsc5(i, i)>=0];
    for j = i + 1:nbus_vsc5
        con_vsc5 = [con_vsc5, cc_vsc5(i, j) == cc_vsc5(j, i)];
        con_vsc5 = [con_vsc5, ss_vsc5(i, j) + ss_vsc5(j, i) == 0];
    end
end

% % -------------------------------------------------------------------------
% % 4.2.6 Obj
% % -------------------------------------------------------------------------
% obj_vsc = sum(pn_vsc1.^2) + sum(pn_vsc2.^2) + sum(pn_vsc3.^2) ...
%     + sum(pn_vsc4.^2) + sum(pn_vsc5.^2)


%% S5 Coupling constraints
con_cp = [];

con_cp = [con_cp; pa2v_ac(4) + pv2a_vsc1 == 0;
     qa2v_ac(4) + qv2a_vsc1 == 0;
     cc_ac(4, 4) == cc_vsc1(1, 1)];

con_cp = [con_cp; pa2v_ac(6) + pv2a_vsc2 == 0;
     qa2v_ac(6) + qv2a_vsc2 == 0;
     cc_ac(6, 6) == cc_vsc2(1, 1)];

con_cp = [con_cp; pa2v_ac(8) + pv2a_vsc3 == 0;
     qa2v_ac(8) + qv2a_vsc3 == 0;
     cc_ac(8, 8) == cc_vsc3(1, 1)];

con_cp = [con_cp; pr2v_res1 + pv2a_vsc4 == 0;
     qr2v_res1 + qv2a_vsc4 == 0];
 
con_cp = [con_cp; pr2v_res2 + pv2a_vsc5 == 0;
     qr2v_res2 + qv2a_vsc5 == 0];

con_cp = [con_cp; pn_dc(1) - pv2d_vsc1 == 0;
    pn_dc(2) - pv2d_vsc2 == 0;
    pn_dc(3) - pv2d_vsc3 == 0;
    pn_dc(4) - pv2d_vsc4 == 0;
    pn_dc(5) - pv2d_vsc5 == 0];


%% S6 Merge
con_vsc     = [con_vsc1; con_vsc2; con_vsc3; con_vsc4; con_vsc5];
con_res     = [con_res1; con_res2];
con         = [con_ac; con_dc; con_vsc; con_res; con_cp];
obj         = obj_ac;

opt = sdpsettings('solver','gurobi', 'verbose',1);
sol = optimize(con, obj, opt);

if sol.problem ~= 0
    disp('Optimization failed.');
    disp(sol.info);
    yalmiperror(sol.problem);
    return;
end

fprintf('\n============================================================\n');
fprintf('Optimization succeeded.\n');
fprintf('Total generation cost: %.6f\n', value(obj));
fprintf('Solver time: %.6f s\n', sol.solvertime);
fprintf('============================================================\n');

fprintf('\nAC generator outputs:\n');
fprintf('  %-8s %-12s %-12s\n', 'Gen', 'Pg(MW)', 'Qg(Mvar)');

for g = 1:ngen_ac
    bus_id = gen_ac(g, 1);

    Pg_val = value(pgen_ac(bus_id)) * baseMVA_ac;
    Qg_val = value(qgen_ac(bus_id)) * baseMVA_ac;

    fprintf('  %-8d %-12.6f %-12.6f\n', ...
        g, Pg_val, Qg_val);
end