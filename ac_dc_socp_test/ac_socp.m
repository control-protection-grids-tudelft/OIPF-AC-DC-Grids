 clc; clear; close all;

%% ========================================================================
%% SOCP-relaxed AC power flow reconstruction for comparison with MATPOWER
%% ========================================================================

%% S1. Load AC network parameters
mpc = case14s;                       
results = runpf(mpc);
bus     = results.bus;
gen     = results.gen;
branch  = results.branch;
baseMVA = results.baseMVA;

[Ybus, ~, ~] = makeYbus(baseMVA, bus, branch);
Ybus = full(Ybus);

G = real(Ybus);
B = imag(Ybus);

n_bus    = size(bus, 1);
n_gen    = size(gen, 1);
n_branch = size(branch, 1);

pd = bus(:, 3) / baseMVA;
qd = bus(:, 4) / baseMVA;

%% S2. Define AC variables
ref_bus = find(bus(:, 2) == 3);     
pv_bus  = find(bus(:, 2) == 2);    
pq_bus  = find(bus(:, 2) == 1);     
gen_bus = gen(:,1);               

Pg_bus = zeros(n_bus, 1);
Qg_bus = zeros(n_bus, 1);

for k = 1:n_gen
    bidx = gen(k, 1);   
    Pg_bus(bidx) = Pg_bus(bidx) + gen(k,2) / baseMVA;
    Qg_bus(bidx) = Qg_bus(bidx) + gen(k,3) / baseMVA;
end

%% S7. Construct fixed net injections
Pinj = Pg_bus - pd;
Qinj = Qg_bus - qd;

disp('Net active power injections Pinj (p.u.):');
disp(Pinj.');

disp('Net reactive power injections Qinj (p.u.):');
disp(Qinj.');

disp(['Sum of Pinj over all buses (p.u.): ', num2str(sum(Pinj))]);
disp(['Sum of Qinj over all buses (p.u.): ', num2str(sum(Qinj))]);

%% S8. Define optimization variables
% p, q   : nodal injections
% cc(i,j): e_i*e_j + f_i*f_j
% ss(i,j): e_i*f_j - e_j*f_i
pn  = sdpvar(n_bus, 1);
qn  = sdpvar(n_bus, 1);
pl = sdpvar(n_bus, n_bus, 'full');
ql = sdpvar(n_bus, n_bus, 'full');
cc = sdpvar(n_bus, n_bus, 'full');
ss = sdpvar(n_bus, n_bus, 'full');

%% S9. Construct SOCP-relaxed AC power flow model
con = [];

% -------------------------------------------------------------------------
%   Branch active/reactive power flow equations
%   Polygonal approximation is used to enforce |S_ij| <= Smax(i,j),
%   i.e., p_ij^2 + q_ij^2 <= Smax(i,j)^2
% -------------------------------------------------------------------------
Smax = zeros(n_bus, n_bus);

for k = 1:n_branch
    i = branch(k, 1);
    j = branch(k, 2);

    rateA = branch(k, 6);

    if rateA == 0
        Smax_val = 1e3;   
    else
        Smax_val = rateA / baseMVA;
    end

    Smax(i, j) = Smax_val;
    Smax(j, i) = Smax_val;   
end

N = 16;
for k = 1:n_branch
    i = branch(k, 1);
    j = branch(k, 2);

    con = [con, pl(i, j) == -G(i, j)*(cc(i, i)-cc(i, j))-B(i, j)*ss(i, j)];
    con = [con, ql(i, j) == B(i, j)*(cc(i, i)-cc(i, j))-G(i, j)*ss(i, j)];
    % ij and ji
    con = [con, pl(j, i) == -G(j, i)*(cc(j, j)-cc(j, i))-B(j, i)*ss(j, i)];
    con = [con, ql(j, i) == B(j, i)*(cc(j, j)-cc(j, i))-G(j, i)*ss(j, i)];
    
    for k = 1: N
        con = [con, cos(k/N*pi)*pl(i, j) + sin(k/N*pi)*ql(i, j) <= Smax(i, j)];
        con = [con; -(cos(k/N*pi)*pl(i, j) + sin(k/N*pi)*ql(i, j)) <= Smax(i, j)];
        % ij and ji
        con = [con, cos(k/N*pi)*pl(j, i) + sin(k/N*pi)*ql(j, i) <= Smax(j, i)];
        con = [con; -(cos(k/N*pi)*pl(j, i) + sin(k/N*pi)*ql(j, i)) <= Smax(j, i)];
    end
    
end

% -------------------------------------------------------------------------
% 9.1 Nodal active/reactive power flow equations
% -------------------------------------------------------------------------
for i = 1:n_bus
    con = [con, pn(i) == ss(i, i)*B(i, i)+sum(cc(i, :).*G(i, :) ...
        - ss(i, :).*B(i, :))];
    con = [con, qn(i) == ss(i, i)*G(i, i)-sum(cc(i, :).*B(i, :) ...
        + ss(i, :).*G(i, :))];
end


% -------------------------------------------------------------------------
% 9.2 Fix injections of all non-slack buses
% -------------------------------------------------------------------------
for i = 1:n_bus
    if ~ismember(i, ref_bus)
        con = [con, pn(i) == Pinj(i)];
        con = [con, qn(i) == Qinj(i)];
    end
end

% -------------------------------------------------------------------------
% 9.3 Symmetry constraints
% -------------------------------------------------------------------------
for i = 1:n_bus
    con = [con, ss(i, i) == 0];  
    con=[con, cc(i, i)>=0];
    for j = i+1:n_bus
        con = [con, cc(i, j) == cc(j, i)];
        con = [con, ss(i, j) + ss(j, i) == 0];
    end
end

% -------------------------------------------------------------------------
% 9.4 SOCP relaxation constraints
%     Add only on existing branches 
% -------------------------------------------------------------------------
for k = 1:n_branch
    i = branch(k, 1);
    j = branch(k, 2);
    con = [con, cc(i, j)^2 + ss(i, j)^2 <= cc(i, i)*cc(j, j)];
end

% -------------------------------------------------------------------------
% 9.5 Voltage magnitude bounds
% -------------------------------------------------------------------------
for i = 1:n_bus
    Vmax = bus(i,12);
    Vmin = bus(i,13);
    con = [con, Vmin^2 <= cc(i,i) <= Vmax^2];
end

% -------------------------------------------------------------------------
% 9.6 Fix voltage magnitude at all generator buses that are slack or PV buses
% -------------------------------------------------------------------------
r = ref_bus.'
gen_idx = find(gen(:,1) == r, 1, 'first');
Vref = gen(gen_idx,6);
con = [con, cc(r, r) == Vref^2];

%% S10. Objective function
% Only a mild auxiliary objective 
obj = 0;
for r = ref_bus.'
    obj = obj + pn(r)^2 + qn(r)^2;
end

%% S11. Solver settings
opt=sdpsettings('savesolveroutput',1,'solver','gurobi');
op=optimize(con,obj,opt);

%% S12. Solve
sol = optimize(con, obj, opt);

%% S13. Check results
if sol.problem ~= 0
    disp('Optimization failed.');
    disp(sol.info);
    yalmiperror(sol.problem);
    return;
end

disp('Optimization succeeded.');

%% S14. Extract SOCP solution
cc_val = value(cc);
ss_val = value(ss);
p_val  = value(pn);
q_val  = value(qn);
pl_val = value(pl);
ql_val = value(ql);

Vm_socp = sqrt(diag(cc_val)); 
Vm_pf = bus(:,8); 

%% S15. Display nodal voltage comparison
disp('--------------------------------------------------');
disp('Voltage magnitude comparison (p.u.):');
disp(table((1:n_bus)', Vm_pf, Vm_socp, Vm_socp - Vm_pf, ...
    'VariableNames', {'Bus','Vm_MATPOWER','Vm_SOCP','Error'}));
disp('--------------------------------------------------');

max_vm_err = max(abs(Vm_socp - Vm_pf));
disp(['Maximum voltage magnitude error = ', num2str(max_vm_err)]);

%% S16. Display branch flow comparison
Pf_pf = branch(:,14) / baseMVA;   
Qf_pf = branch(:,15) / baseMVA;  
Pt_pf = branch(:,16) / baseMVA;   
Qt_pf = branch(:,17) / baseMVA;  

Pij_socp = zeros(n_branch,1);
Qij_socp = zeros(n_branch,1);
Pji_socp = zeros(n_branch,1);
Qji_socp = zeros(n_branch,1);

for k = 1:n_branch
    i = branch(k,1);
    j = branch(k,2);

    % SOCP reconstructed from-end flow: i -> j
    Pij_socp(k) = pl_val(i, j)

    Qij_socp(k) = ql_val(i, j)

    % SOCP reconstructed to-end flow: j -> i
    % Note: ss(j,i) = -ss(i,j), cc(j,i) = cc(i,j)
    Pji_socp(k) = pl_val(j, i);

    Qji_socp(k) = ql_val(j, i);
end

Pf_err = Pij_socp - Pf_pf;
Qf_err = Qij_socp - Qf_pf;
Pt_err = Pji_socp - Pt_pf;
Qt_err = Qji_socp - Qt_pf;

disp('--------------------------------------------------');
disp('Branch power flow comparison (from-end, p.u.):');
branch_table_from = table((1:n_branch)', branch(:,1), branch(:,2), ...
    Pf_pf, Pij_socp, Pf_err, ...
    Qf_pf, Qij_socp, Qf_err, ...
    'VariableNames', {'Branch','FromBus','ToBus', ...
    'Pf_MATPOWER','Pf_SOCP','Pf_Error', ...
    'Qf_MATPOWER','Qf_SOCP','Qf_Error'});
disp(branch_table_from);

disp('--------------------------------------------------');
disp('Branch power flow comparison (to-end, p.u.):');
branch_table_to = table((1:n_branch)', branch(:,1), branch(:,2), ...
    Pt_pf, Pji_socp, Pt_err, ...
    Qt_pf, Qji_socp, Qt_err, ...
    'VariableNames', {'Branch','FromBus','ToBus', ...
    'Pt_MATPOWER','Pt_SOCP','Pt_Error', ...
    'Qt_MATPOWER','Qt_SOCP','Qt_Error'});
disp(branch_table_to);
disp('--------------------------------------------------');

max_pf_err = max(abs(Pf_err));
max_qf_err = max(abs(Qf_err));
max_pt_err = max(abs(Pt_err));
max_qt_err = max(abs(Qt_err));

avg_pf_err = mean(abs(Pf_err));
avg_qf_err = mean(abs(Qf_err));
avg_pt_err = mean(abs(Pt_err));
avg_qt_err = mean(abs(Qt_err));

disp('Branch flow error summary (p.u.):');
fprintf('From-end active power   max error = %.8e, avg error = %.8e\n', max_pf_err, avg_pf_err);
fprintf('From-end reactive power max error = %.8e, avg error = %.8e\n', max_qf_err, avg_qf_err);
fprintf('To-end active power     max error = %.8e, avg error = %.8e\n', max_pt_err, avg_pt_err);
fprintf('To-end reactive power   max error = %.8e, avg error = %.8e\n', max_qt_err, avg_qt_err);

%% S17. Check SOCP relaxation gap on branches
relax_gap = zeros(n_branch, 1);

for k = 1:n_branch
    i = branch(k, 1);
    j = branch(k, 2);

    lhs = cc_val(i, j)^2 + ss_val(i, j)^2;
    rhs = cc_val(i, i)*cc_val(j, j);

    relax_gap(k) = rhs - lhs;
end

disp(['Maximum SOCP relaxation gap on branches = ', num2str(max(abs(relax_gap)))]);
disp(['Average SOCP relaxation gap on branches = ', num2str(mean(abs(relax_gap)))]);

%% S18. Voltage comparison
figure;
plot(1:n_bus, Vm_pf, '-o', 'LineWidth', 1.2); hold on;
plot(1:n_bus, Vm_socp, '-xr', 'LineWidth', 1.2);
xlabel('Bus index');
ylabel('Voltage magnitude (p.u.)');
legend('MATPOWER runpf', 'SOCP reconstruction', 'Location', 'best');
title('Voltage magnitude comparison');
grid on;

%% S19. Voltage error
figure;
stem(1:n_bus, Vm_socp - Vm_pf, 'LineWidth', 1.2);
xlabel('Bus index');
ylabel('Voltage magnitude error (p.u.)');
title('SOCP&MATPOWER voltage magnitude error');
grid on;