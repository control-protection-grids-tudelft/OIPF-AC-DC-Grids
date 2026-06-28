clc; clear; close all;

%% ========================================================================
%% SOCP-relaxed DC power flow reconstruction for comparison with MATPOWER
%% ========================================================================

%% S1. Load DC network parameters
mpc = case5_dc;                       
results = runpf(mpc);
bus     = results.bus;
gen     = results.gen;
branch  = results.branch;
baseMVA = results.baseMVA;

[Ybus, ~, ~] = makeYbus(baseMVA, bus, branch);
Ybus = full(Ybus);

n_bus    = size(bus, 1);
n_gen    = size(gen, 1);
n_branch = size(branch, 1);

zij = 1 ./ abs(Ybus);
zij = zij - diag(diag(zij));
zij(isinf(zij)) = 1e4;

bus(:, 3) = bus(:, 3) / baseMVA; 
p_d = bus(:, 3);

n_bus    = length(bus(:, 1));
n_branch = size(branch, 1);

fbus = branch(:, 1);
tbus = branch(:, 2);

ref_bus = find(bus(:, 2) == 3);     

%% S2. Define optimization variables (keep original names)
u     = sdpvar(n_bus, 1);              % voltage square
lij   = sdpvar(n_bus, n_bus, 'full');   % current square
pn = sdpvar(n_bus, 1);       
pij   = sdpvar(n_bus, n_bus, 'full');  

%% S3. Construct SOCP model
con=[];

% -------------------------------------------------------------------------
% 3.0 Nodal power balance
% -------------------------------------------------------------------------
for i = 1:n_bus
    expr = 0;
    for k = 1:n_branch
        if fbus(k) == i
            j = tbus(k);
            expr = expr + pij(i, j);
        elseif tbus(k) == i
            j = fbus(k);
            expr = expr + pij(i, j);
        end
    end
    con = [con; pn(i) == expr];
end

% -------------------------------------------------------------------------
% 3.1 Branch equations
% -------------------------------------------------------------------------
for k = 1:n_branch
    i = fbus(k);
    j = tbus(k);
    rateA = branch(k, 6);

    if rateA == 0
        Smax_val = 1e3;   
    else
        Smax_val = rateA / baseMVA;
    end

    % loss
    con = [con;
        pij(i, j) + pij(j, i) == zij(i, j) * lij(i, j)];

    % volt drop
    con = [con;
        u(i) - u(j) == zij(i, j) * (pij(i, j) - pij(j, i))];

    % SOCP
    con = [con;
        pij(i, j)^2 <= lij(i, j) * u(i);
        pij(j, i)^2 <= lij(j, i) * u(i)];

    con = [con;
        lij(i, j) >= 0
        lij(j, i) >= 0];

    con = [con; 
        u(i) >= 0]

    con = [con;
        lij(j, i) == lij(i, j)];

    % Bounds
    con = [con;
        pij(i, j) <= Smax_val
        pij(j, i) <= Smax_val];
end

% -------------------------------------------------------------------------
% 3.2 Power Injection
% -------------------------------------------------------------------------
for i = 1:n_bus
    if ~ismember(i, ref_bus)
        con = [con; pn(i) == -p_d(i)];
    end
end

% -------------------------------------------------------------------------
% 3.3 Slack bus voltage
% -------------------------------------------------------------------------
r = ref_bus.'
gen_idx = find(gen(:,1) == r, 1, 'first');
Vref = gen(gen_idx,6);
con = [con, u(r) == Vref^2];

% -------------------------------------------------------------------------
% 3.4 Voltage magnitude bounds
% -------------------------------------------------------------------------
for i = 1:n_bus
    Vmax = bus(i,12);
    Vmin = bus(i,13);
    con = [con, Vmin^2 <= u(i) <= Vmax^2];
end

%% S4. Objective
% Only a mild auxiliary objective 
obj = (pn);

%% S5. Solve
opt = sdpsettings('solver','cplex','verbose',1);
sol = optimize(con, obj, opt);

if sol.problem ~= 0
    disp('Optimization failed.');
    disp(sol.info);
    return;
end

%% S6. Extract SOCP results
u_val   = value(u);
p_val   = value(pn);
pij_val = value(pij);
lij_val = value(lij);

v_opt = sqrt(u_val);

%% S7. MATPOWER reference
[~,bus_pf,~,branch_pf,~] = runpf('case5_dc');

v_matpower = bus_pf(:,8);

%% S8. Voltage comparison
disp('--------------------------------------------------');
disp('Voltage magnitude comparison (p.u.):');
disp(table((1:n_bus)', v_matpower, v_opt, v_opt - v_matpower, ...
    'VariableNames', {'Bus','Vm_MATPOWER','Vm_SOCP','Error'}));

max_vm_err = max(abs(v_opt - v_matpower));
disp(['Maximum voltage error = ', num2str(max_vm_err)]);
disp('--------------------------------------------------');

%% S9. Branch flow comparison
Pf_pf = branch_pf(:, 14) / baseMVA;   
Pt_pf = branch_pf(:, 16) / baseMVA;  

Pij_socp = zeros(n_branch, 1);
Pji_socp = zeros(n_branch, 1);

for k = 1:n_branch
    i = fbus(k);
    j = tbus(k);

    Pij_socp(k) = pij_val(i, j);
    Pji_socp(k) = pij_val(j, i);
end

Pf_err = Pij_socp - Pf_pf;
Pt_err = Pji_socp - Pt_pf;

disp('--------------------------------------------------');
disp('Branch power comparison (from-end):');
disp(table((1:n_branch)', fbus, tbus, ...
    Pf_pf, Pij_socp, Pf_err, ...
    'VariableNames', {'Branch','From','To','Pf_MP','Pf_SOCP','Error'}));

disp('--------------------------------------------------');
disp('Branch power comparison (to-end):');
disp(table((1:n_branch)', fbus, tbus, ...
    Pt_pf, Pji_socp, Pt_err, ...
    'VariableNames', {'Branch','From','To','Pt_MP','Pt_SOCP','Error'}));

max_pf_err = max(abs(Pf_err));
max_pt_err = max(abs(Pt_err));

fprintf('Max Pf error = %.4e\n', max_pf_err);
fprintf('Max Pt error = %.4e\n', max_pt_err);

%% S10. SOCP relaxation gap
gap = zeros(n_branch,1);

for k = 1:n_branch
    i = fbus(k);
    j = tbus(k);

    lhs = pij_val(i, j)^2;
    rhs = lij_val(i, j) * u_val(i);

    gap(k) = rhs - lhs;
end

disp(['Max SOCP gap = ', num2str(max(abs(gap)))]);
disp(['Avg SOCP gap = ', num2str(mean(abs(gap)))]);

%% S11. Plot voltage
figure;
plot(1:n_bus, v_matpower, '-o'); hold on;
plot(1:n_bus, v_opt, '-xr');
xlabel('Bus index');
ylabel('Voltage magnitude error (p.u.)');
legend('MATPOWER','SOCP');
title('Voltage comparison');
grid on;