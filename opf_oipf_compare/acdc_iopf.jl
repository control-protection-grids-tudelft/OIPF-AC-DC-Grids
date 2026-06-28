"""
Reference:https://arxiv.org/abs/2606.09436
"""

# =============================================================================
# AC/DC IOPF model 
# =============================================================================

using JuMP
using Gurobi
using LinearAlgebra
using SparseArrays
using Printf
include("parse_mpc.jl")
include("makeYbus.jl")

# -----------------------------------------------------------------------------
# Additional functions
# -----------------------------------------------------------------------------

"""
Add polygonal approximation constraints for apparent power capacity:
    |p + jq| <= Smax
"""
function add_polygon_capacity!(model::Model, p, q, Smax, N::Int)
    for kk in 1:N
        ck = cos(kk / N * pi)
        sk = sin(kk / N * pi)
        @constraint(model,  ck * p + sk * q <= Smax)
        @constraint(model, -ck * p - sk * q <= Smax)
    end
end


"""
Describe the coupling logic between information flow and power flow
"""
function add_couple_logic_constraints!(model::Model, e, x, x0)
    @constraint(model, e >= x - x0)
    @constraint(model, e >= x0 - x)
    @constraint(model, e <= x + x0)
    @constraint(model, e <= 2 - (x0 + x))
end

function add_vsc_nodal_constraints!(
    model::Model,
    nbus_vsc::Int,
    G_vsc,
    B_vsc,
    pn_vsc,
    qn_vsc,
    cc_vsc,
    ss_vsc;
    Vmin::Float64 = 0.9,
    Vmax::Float64 = 1.1,
)
    for i in 1:nbus_vsc
        @constraint(model,
            pn_vsc[i] == ss_vsc[i, i] * B_vsc[i, i]
                + sum(cc_vsc[i, j] * G_vsc[i, j] - ss_vsc[i, j] * B_vsc[i, j] for j in 1:nbus_vsc)
        )

        @constraint(model,
            qn_vsc[i] == ss_vsc[i, i] * G_vsc[i, i]
                - sum(cc_vsc[i, j] * B_vsc[i, j] + ss_vsc[i, j] * G_vsc[i, j] for j in 1:nbus_vsc)
        )

        @constraint(model, Vmin^2 <= cc_vsc[i, i] <= Vmax^2)
    end
end

function add_vsc_socp_constraints!(
    model::Model,
    nbus_vsc::Int,
    nbranch_vsc::Int,
    branch_vsc,
    cc_vsc,
    ss_vsc,
)
    for kk in 1:nbranch_vsc
        i = Int(branch_vsc[kk, 1])
        j = Int(branch_vsc[kk, 2])
        @constraint(model, cc_vsc[i, j]^2 + ss_vsc[i, j]^2 <= cc_vsc[i, i] * cc_vsc[j, j])
    end

    for i in 1:nbus_vsc
        @constraint(model, ss_vsc[i, i] == 0)
        @constraint(model, cc_vsc[i, i] >= 0)
        for j in (i + 1):nbus_vsc
            @constraint(model, cc_vsc[i, j] == cc_vsc[j, i])
            @constraint(model, ss_vsc[i, j] + ss_vsc[j, i] == 0)
        end
    end
end

function main()
    # =============================================================================
    # Build model
    # =============================================================================

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 1)

    # -----------------------------------------------------------------------------
    # 1.1.1 AC grid parameters
    # -----------------------------------------------------------------------------
    mpc_ac = parse_mpc(joinpath(@__DIR__, "case9.m"))
    bus_ac     = mpc_ac["bus"]
    gen_ac     = mpc_ac["gen"]
    branch_ac  = mpc_ac["branch"]
    baseMVA_ac = mpc_ac["baseMVA"]
    gencost_ac = mpc_ac["gencost"]

    Ybus_ac = makeYbus(baseMVA_ac, bus_ac, branch_ac)
    Ybus_ac = Matrix(Ybus_ac)

    G_ac = real.(Ybus_ac)
    B_ac = imag.(Ybus_ac)

    nbus_ac    = size(bus_ac, 1)
    ngen_ac    = size(gen_ac, 1)
    nbranch_ac = size(branch_ac, 1)

    pd_ac = bus_ac[:, 3] ./ baseMVA_ac
    qd_ac = bus_ac[:, 4] ./ baseMVA_ac

    refBus_ac = findall(bus_ac[:, 2] .== 3)
    genBus_ac = Int.(gen_ac[:, 1])
    a2vBus_ac = [4, 6, 8]

    # -----------------------------------------------------------------------------
    # 1.1.2 AC grid variables
    # -----------------------------------------------------------------------------
    @variable(model, pn_ac[1:nbus_ac])
    @variable(model, qn_ac[1:nbus_ac])
    @variable(model, pgen_ac[1:nbus_ac])
    @variable(model, qgen_ac[1:nbus_ac])
    @variable(model, pa2v_ac[1:nbus_ac])
    @variable(model, qa2v_ac[1:nbus_ac])
    @variable(model, pl_ac[1:nbus_ac, 1:nbus_ac])
    @variable(model, ql_ac[1:nbus_ac, 1:nbus_ac])
    @variable(model, cc_ac[1:nbus_ac, 1:nbus_ac])
    @variable(model, ss_ac[1:nbus_ac, 1:nbus_ac])

    # =============================================================================
    # S1.2 Add SOCP-relaxed AC network constraints
    # =============================================================================

    # -----------------------------------------------------------------------------
    # 1.2.1 Branch active/reactive power flow equations
    # -----------------------------------------------------------------------------
    Smax_ac = zeros(nbus_ac, nbus_ac)

    for kk in 1:nbranch_ac
        i = Int(branch_ac[kk, 1])
        j = Int(branch_ac[kk, 2])
        rateA = branch_ac[kk, 6]
        Smax_val = rateA == 0 ? 1e3 : rateA / baseMVA_ac
        Smax_ac[i, j] = Smax_val
        Smax_ac[j, i] = Smax_val
    end

    N = 16
    for kk in 1:nbranch_ac
        i = Int(branch_ac[kk, 1])
        j = Int(branch_ac[kk, 2])

        @constraint(model, pl_ac[i, j] == -G_ac[i, j] * (cc_ac[i, i] - cc_ac[i, j]) - B_ac[i, j] * ss_ac[i, j])
        @constraint(model, ql_ac[i, j] ==  B_ac[i, j] * (cc_ac[i, i] - cc_ac[i, j]) - G_ac[i, j] * ss_ac[i, j])

        @constraint(model, pl_ac[j, i] == -G_ac[j, i] * (cc_ac[j, j] - cc_ac[j, i]) - B_ac[j, i] * ss_ac[j, i])
        @constraint(model, ql_ac[j, i] ==  B_ac[j, i] * (cc_ac[j, j] - cc_ac[j, i]) - G_ac[j, i] * ss_ac[j, i])

        add_polygon_capacity!(model, pl_ac[i, j], ql_ac[i, j], Smax_ac[i, j], N)
        add_polygon_capacity!(model, pl_ac[j, i], ql_ac[j, i], Smax_ac[j, i], N)
    end

    # -----------------------------------------------------------------------------
    # 1.2.2 Nodal active/reactive power flow equations
    # -----------------------------------------------------------------------------
    for i in 1:nbus_ac
        @constraint(model,
            pn_ac[i] == ss_ac[i, i] * B_ac[i, i]
                + sum(cc_ac[i, j] * G_ac[i, j] - ss_ac[i, j] * B_ac[i, j] for j in 1:nbus_ac)
        )
        @constraint(model,
            qn_ac[i] == ss_ac[i, i] * G_ac[i, i]
                - sum(cc_ac[i, j] * B_ac[i, j] + ss_ac[i, j] * G_ac[i, j] for j in 1:nbus_ac)
        )
    end

    # -----------------------------------------------------------------------------
    # 1.2.3 Fix injections
    # -----------------------------------------------------------------------------
    for i in 1:nbus_ac
        p_expr = -pd_ac[i]
        q_expr = -qd_ac[i]

        if i in genBus_ac
            p_expr += pgen_ac[i]
            q_expr += qgen_ac[i]
        end

        if i in a2vBus_ac
            p_expr -= pa2v_ac[i]
            q_expr -= qa2v_ac[i]
        end

        @constraint(model, pn_ac[i] == p_expr)
        @constraint(model, qn_ac[i] == q_expr)
    end

    # -----------------------------------------------------------------------------
    # 1.2.4 SOCP relaxation constraints and symmetry constraints
    # -----------------------------------------------------------------------------
    for kk in 1:nbranch_ac
        i = Int(branch_ac[kk, 1])
        j = Int(branch_ac[kk, 2])
        @constraint(model, cc_ac[i, j]^2 + ss_ac[i, j]^2 <= cc_ac[i, i] * cc_ac[j, j])
    end

    for i in 1:nbus_ac
        @constraint(model, ss_ac[i, i] == 0)
        @constraint(model, cc_ac[i, i] >= 0)
        for j in (i + 1):nbus_ac
            @constraint(model, cc_ac[i, j] == cc_ac[j, i])
            @constraint(model, ss_ac[i, j] + ss_ac[j, i] == 0)
        end
    end

    # -----------------------------------------------------------------------------
    # 1.2.5 Generator active/reactive power limits
    # -----------------------------------------------------------------------------
    pgenMin_ac = zeros(nbus_ac)
    pgenMax_ac = zeros(nbus_ac)
    qgenMin_ac = zeros(nbus_ac)
    qgenMax_ac = zeros(nbus_ac)

    for i in 1:nbus_ac
        gen_idx = findfirst(gen_ac[:, 1] .== i)
        if gen_idx !== nothing
            qgenMax_ac[i] = gen_ac[gen_idx, 4]  / baseMVA_ac
            qgenMin_ac[i] = gen_ac[gen_idx, 5]  / baseMVA_ac
            pgenMax_ac[i] = gen_ac[gen_idx, 9]  / baseMVA_ac
            pgenMin_ac[i] = gen_ac[gen_idx, 10] / baseMVA_ac
        end
    end

    for i in 1:nbus_ac
        if i in genBus_ac
            @constraint(model, pgenMin_ac[i] <= pgen_ac[i] <= pgenMax_ac[i])
            @constraint(model, qgenMin_ac[i] <= qgen_ac[i] <= qgenMax_ac[i])
        else
            @constraint(model, pgen_ac[i] == 0)
            @constraint(model, qgen_ac[i] == 0)
        end
    end

    # -----------------------------------------------------------------------------
    # 1.2.6 Voltage magnitude bounds
    # -----------------------------------------------------------------------------
    for i in 1:nbus_ac
        Vmax_ac = bus_ac[i, 12]
        Vmin_ac = bus_ac[i, 13]
        @constraint(model, Vmin_ac^2 <= cc_ac[i, i] <= Vmax_ac^2)
    end

    # -----------------------------------------------------------------------------
    # 1.2.7 Voltage references
    # -----------------------------------------------------------------------------
    r = refBus_ac[1]
    gen_idx = findfirst(gen_ac[:, 1] .== r)
    Vref_ac = gen_ac[gen_idx, 6]
    @constraint(model, cc_ac[r, r] == Vref_ac^2)

    # -----------------------------------------------------------------------------
    # 1.2.8 Objective: generation costs
    # -----------------------------------------------------------------------------
    c2_ac = zeros(nbus_ac)
    c1_ac = zeros(nbus_ac)
    c0_ac = zeros(nbus_ac)

    for i in 1:nbus_ac
        gen_idx = findfirst(gen_ac[:, 1] .== i)
        if gen_idx !== nothing
            ncost = Int(gencost_ac[gen_idx, 4])
            if ncost == 3
                c2_ac[i] = gencost_ac[gen_idx, end - 2]
                c1_ac[i] = gencost_ac[gen_idx, end - 1]
                c0_ac[i] = gencost_ac[gen_idx, end]
            elseif ncost == 2
                c2_ac[i] = 0
                c1_ac[i] = gencost_ac[gen_idx, end - 1]
                c0_ac[i] = gencost_ac[gen_idx, end]
            end
        end
    end

    obj_ac = @expression(model,
        sum((i in genBus_ac) ?
            (c2_ac[i] * (pgen_ac[i] * baseMVA_ac)^2 + c1_ac[i] * (pgen_ac[i] * baseMVA_ac) + c0_ac[i]) :
            0
            for i in 1:nbus_ac)
    )

    # =============================================================================
    # S2.1 Load DC network parameters & variables
    # =============================================================================

    mpc_dc = parse_mpc(joinpath(@__DIR__, "case5_dc.m"))
    results_dc = mpc_dc
    bus_dc      = results_dc["bus"]
    gen_dc      = results_dc["gen"]
    branch_dc   = results_dc["branch"]
    baseMW_dc   = results_dc["baseMVA"]

    Ybus_dc = makeYbus(baseMW_dc, bus_dc, branch_dc)
    Ybus_dc = Matrix(Ybus_dc)

    nbus_dc    = size(bus_dc, 1)
    ngen_dc    = size(gen_dc, 1)
    nbranch_dc = size(branch_dc, 1)

    zij_dc = 1.0 ./ abs.(Ybus_dc)
    zij_dc = zij_dc - Diagonal(diag(zij_dc))
    zij_dc[isinf.(zij_dc)] .= 1e5

    bus_dc[:, 3] ./= baseMW_dc
    pd_dc = bus_dc[:, 3]

    fbus_dc = Int.(branch_dc[:, 1])
    tbus_dc = Int.(branch_dc[:, 2])

    refBus_dc = findall(bus_dc[:, 2] .== 3)

    @variable(model, u_dc[1:nbus_dc])
    @variable(model, lij_dc[1:nbus_dc, 1:nbus_dc])
    @variable(model, pn_dc[1:nbus_dc])
    @variable(model, pij_dc[1:nbus_dc, 1:nbus_dc])

    @variable(model, alpha_dc[1:nbus_dc, 1:nbus_dc], Bin)
    @variable(model, beta_dc[1:nbus_dc, 1:nbus_dc], Bin)
    @variable(model, ep_dc[1:nbus_dc, 1:nbus_dc], Bin)
    @variable(model, bij_dc[1:nbus_dc, 1:nbus_dc])
    @variable(model, tij_dc[1:nbus_dc, 1:nbus_dc])

    # =============================================================================
    # S2.2 Add DC network constraints
    # =============================================================================
    M = 1e5

    # -----------------------------------------------------------------------------
    # 2.2.1 Nodal power balance
    # -----------------------------------------------------------------------------
    for i in 1:nbus_dc
        @constraint(model,
            pn_dc[i] == sum(fbus_dc[kk] == i ? pij_dc[i, tbus_dc[kk]] :
                            tbus_dc[kk] == i ? pij_dc[i, fbus_dc[kk]] : 0
                            for kk in 1:nbranch_dc)
        )
    end

    # -----------------------------------------------------------------------------
    # 2.2.2 Branch equations
    # -----------------------------------------------------------------------------
    for kk in 1:nbranch_dc
        i = fbus_dc[kk]
        j = tbus_dc[kk]

        rateA = branch_dc[kk, 6]
        Pmax_val = rateA == 0 ? 1e3 : rateA / baseMW_dc
        Pmax_val = 100 / baseMW_dc

        @constraint(model, pij_dc[i, j] + pij_dc[j, i] == zij_dc[i, j] * lij_dc[i, j])
        @constraint(model, pij_dc[i, j]^2 <= lij_dc[i, j] * u_dc[i])
        @constraint(model, pij_dc[j, i]^2 <= lij_dc[j, i] * u_dc[j])
        @constraint(model, lij_dc[i, j] >= 0)
        @constraint(model, lij_dc[j, i] >= 0)
        @constraint(model, u_dc[i] >= 0)
        @constraint(model, lij_dc[j, i] == lij_dc[i, j])
        @constraint(model, pij_dc[i, j] <= Pmax_val)
        @constraint(model, pij_dc[j, i] <= Pmax_val)
    end

    # -----------------------------------------------------------------------------
    # 2.2.5 Voltage magnitude bounds
    # -----------------------------------------------------------------------------
    for i in 1:nbus_dc
        Vmax_dc = bus_dc[i, 12]
        Vmin_dc = bus_dc[i, 13]
        @constraint(model, Vmin_dc^2 <= u_dc[i] <= Vmax_dc^2)
    end

    # -----------------------------------------------------------------------------
    # 2.2.6 Transmission switch
    # -----------------------------------------------------------------------------
    for kk in 1:nbranch_dc
        i = fbus_dc[kk]
        j = tbus_dc[kk]

        @constraint(model, pij_dc[i, j] >= -M * alpha_dc[i, j])
        @constraint(model, pij_dc[i, j] <=  M * alpha_dc[i, j])
        @constraint(model, pij_dc[j, i] >= -M * alpha_dc[j, i])
        @constraint(model, pij_dc[j, i] <=  M * alpha_dc[j, i])

        @constraint(model,
            (u_dc[i] - bij_dc[i, j]) - (u_dc[j] - tij_dc[i, j]) ==
            zij_dc[i, j] * (pij_dc[i, j] - pij_dc[j, i])
        )

        Vmax_dc = max(bus_dc[i, 12], bus_dc[j, 12])
        Vmin_dc = min(bus_dc[i, 13], bus_dc[j, 13])

        @constraint(model, bij_dc[i, j] >= Vmin_dc^2 * (1 - alpha_dc[i, j]))
        @constraint(model, bij_dc[i, j] <= Vmax_dc^2 * (1 - alpha_dc[i, j]))

        @constraint(model, tij_dc[i, j] >= Vmin_dc^2 * (1 - alpha_dc[i, j]))
        @constraint(model, tij_dc[i, j] <= Vmax_dc^2 * (1 - alpha_dc[i, j]))

        @constraint(model, u_dc[i] - bij_dc[i, j] >= Vmin_dc^2 * alpha_dc[i, j])
        @constraint(model, u_dc[i] - bij_dc[i, j] <= Vmax_dc^2 * alpha_dc[i, j])

        @constraint(model, u_dc[j] - tij_dc[i, j] >= Vmin_dc^2 * alpha_dc[i, j])
        @constraint(model, u_dc[j] - tij_dc[i, j] <= Vmax_dc^2 * alpha_dc[i, j])

        @constraint(model, alpha_dc[i, j] - alpha_dc[j, i] == 0)
        @constraint(model, tij_dc[i, j] - tij_dc[j, i] == 0)
        @constraint(model, bij_dc[i, j] - bij_dc[j, i] == 0)

        @constraint(model, alpha_dc[1, 4] == 0)
        @constraint(model, alpha_dc[3, 5] == 0)

        @constraint(model, alpha_dc[i, j] <= beta_dc[i, j])
        @constraint(model, alpha_dc[i, j] <= beta_dc[j, i])
        @constraint(model, alpha_dc[i, j] >= beta_dc[i, j] + beta_dc[j, i] - 1)
    end

    # -----------------------------------------------------------------------------
    # 2.2.7 Objective: switching costs
    # -----------------------------------------------------------------------------
    betak_dc = zeros(nbus_dc, nbus_dc)
    lam_dc   = zeros(nbus_dc, nbus_dc)

    betak_dc[1, 2] = 1; lam_dc[1, 2] = 150
    betak_dc[2, 1] = 1; lam_dc[2, 1] = 100
    betak_dc[1, 4] = 1; lam_dc[1, 4] = 250
    betak_dc[4, 1] = 1; lam_dc[4, 1] = 300
    betak_dc[1, 5] = 0; lam_dc[1, 5] = 100
    betak_dc[5, 1] = 0; lam_dc[5, 1] = 100
    betak_dc[2, 3] = 1; lam_dc[2, 3] = 200
    betak_dc[3, 2] = 1; lam_dc[3, 2] = 300
    betak_dc[2, 4] = 0; lam_dc[2, 4] = 350
    betak_dc[4, 2] = 0; lam_dc[4, 2] = 250
    betak_dc[2, 5] = 0; lam_dc[2, 5] = 150
    betak_dc[5, 2] = 0; lam_dc[5, 2] = 100
    betak_dc[3, 4] = 0; lam_dc[3, 4] = 350
    betak_dc[4, 3] = 0; lam_dc[4, 3] = 450
    betak_dc[3, 5] = 1; lam_dc[3, 5] = 250
    betak_dc[5, 3] = 1; lam_dc[5, 3] = 350
    betak_dc[4, 5] = 1; lam_dc[4, 5] = 450
    betak_dc[5, 4] = 1; lam_dc[5, 4] = 400

    for kk in 1:nbranch_dc
        i = fbus_dc[kk]
        j = tbus_dc[kk]
        add_couple_logic_constraints!(model, ep_dc[i, j], beta_dc[i, j], betak_dc[i, j])
        add_couple_logic_constraints!(model, ep_dc[j, i], beta_dc[j, i], betak_dc[j, i])
    end

    obj_dc = @expression(model,
        sum(lam_dc[fbus_dc[kk], tbus_dc[kk]] * ep_dc[fbus_dc[kk], tbus_dc[kk]] +
            lam_dc[tbus_dc[kk], fbus_dc[kk]] * ep_dc[tbus_dc[kk], fbus_dc[kk]]
            for kk in 1:nbranch_dc)
    )

    # =============================================================================
    # S3.1 Load RES parameters & variables
    # =============================================================================

    baseMVA_res1 = 100
    poutMax_res1 = 80 / baseMVA_res1
    SoutMax_res1 = 100 / baseMVA_res1

    baseMVA_res2 = 100
    poutMax_res2 = 70 / baseMVA_res2
    SoutMax_res2 = 100 / baseMVA_res2

    @variable(model, pout_res1)
    @variable(model, qout_res1)
    @variable(model, pr2v_res1)
    @variable(model, qr2v_res1)

    @variable(model, pout_res2)
    @variable(model, qout_res2)
    @variable(model, pr2v_res2)
    @variable(model, qr2v_res2)

    # =============================================================================
    # S3.2 RES constraints
    # =============================================================================

    con_res1 = JuMP.ConstraintRef[]
    con_res2 = JuMP.ConstraintRef[]

    N = 16
    add_polygon_capacity!(model, pout_res1, qout_res1, SoutMax_res1, N)
    add_polygon_capacity!(model, pout_res2, qout_res2, SoutMax_res2, N)

    @constraint(model, 0 <= pout_res1 <= poutMax_res1)
    @constraint(model, 0 <= pout_res2 <= poutMax_res2)

    @constraint(model, pout_res1 == pr2v_res1)
    @constraint(model, qout_res1 == qr2v_res1)
    @constraint(model, pout_res2 == pr2v_res2)
    @constraint(model, qout_res2 == qr2v_res2)

    # =============================================================================
    # S4.1 Load VSC network parameters & variables
    # =============================================================================

    # -----------------------------------------------------------------------------
    # 4.1.1 VSC parameters
    # -----------------------------------------------------------------------------
    mpc_vsc1     = parse_mpc(joinpath(@__DIR__, "case2_vsc.m"))
    bus_vsc1     = mpc_vsc1["bus"]
    branch_vsc1  = mpc_vsc1["branch"]
    baseMVA_vsc1 = mpc_vsc1["baseMVA"]
    Ybus_vsc1 = makeYbus(baseMVA_vsc1, bus_vsc1, branch_vsc1)
    Ybus_vsc1 = Matrix(Ybus_vsc1)
    G_vsc1 = real.(Ybus_vsc1)
    B_vsc1 = imag.(Ybus_vsc1)
    nbus_vsc1    = size(bus_vsc1, 1)
    nbranch_vsc1 = size(branch_vsc1, 1)
    pd_vsc1 = bus_vsc1[:, 3] ./ baseMVA_vsc1
    qd_vsc1 = bus_vsc1[:, 4] ./ baseMVA_vsc1
    refBus_vsc1 = findall(bus_vsc1[:, 2] .== 3)

    mpc_vsc2 = parse_mpc(joinpath(@__DIR__, "case2_vsc.m"))
    bus_vsc2     = mpc_vsc2["bus"]
    branch_vsc2  = mpc_vsc2["branch"]
    baseMVA_vsc2 = mpc_vsc2["baseMVA"]
    Ybus_vsc2 = makeYbus(baseMVA_vsc2, bus_vsc2, branch_vsc2)
    Ybus_vsc2 = Matrix(Ybus_vsc2)
    G_vsc2 = real.(Ybus_vsc2)
    B_vsc2 = imag.(Ybus_vsc2)
    nbus_vsc2    = size(bus_vsc2, 1)
    nbranch_vsc2 = size(branch_vsc2, 1)
    pd_vsc2 = bus_vsc2[:, 3] ./ baseMVA_vsc2
    qd_vsc2 = bus_vsc2[:, 4] ./ baseMVA_vsc2
    refBus_vsc2 = findall(bus_vsc2[:, 2] .== 3)

    mpc_vsc3 = parse_mpc(joinpath(@__DIR__, "case2_vsc.m"))
    bus_vsc3     = mpc_vsc3["bus"]
    branch_vsc3  = mpc_vsc3["branch"]
    baseMVA_vsc3 = mpc_vsc3["baseMVA"]
    Ybus_vsc3 = makeYbus(baseMVA_vsc3, bus_vsc3, branch_vsc3)
    Ybus_vsc3 = Matrix(Ybus_vsc3)
    G_vsc3 = real.(Ybus_vsc3)
    B_vsc3 = imag.(Ybus_vsc3)
    nbus_vsc3    = size(bus_vsc3, 1)
    nbranch_vsc3 = size(branch_vsc3, 1)
    pd_vsc3 = bus_vsc2[:, 3] ./ baseMVA_vsc3
    qd_vsc3 = bus_vsc2[:, 4] ./ baseMVA_vsc3
    refBus_vsc3 = findall(bus_vsc3[:, 2] .== 3)

    mpc_vsc4 = parse_mpc(joinpath(@__DIR__, "case2_vsc.m"))
    bus_vsc4     = mpc_vsc4["bus"]
    branch_vsc4  = mpc_vsc4["branch"]
    baseMVA_vsc4 = mpc_vsc4["baseMVA"]
    Ybus_vsc4 = makeYbus(baseMVA_vsc4, bus_vsc4, branch_vsc4)
    Ybus_vsc4 = Matrix(Ybus_vsc4)
    G_vsc4 = real.(Ybus_vsc4)
    B_vsc4 = imag.(Ybus_vsc4)
    nbus_vsc4    = size(bus_vsc4, 1)
    nbranch_vsc4 = size(branch_vsc4, 1)
    pd_vsc4 = bus_vsc2[:, 3] ./ baseMVA_vsc4
    qd_vsc4 = bus_vsc2[:, 4] ./ baseMVA_vsc4
    refBus_vsc4 = findall(bus_vsc4[:, 2] .== 3)

    mpc_vsc5 = parse_mpc(joinpath(@__DIR__, "case2_vsc.m"))
    bus_vsc5     = mpc_vsc5["bus"]
    branch_vsc5  = mpc_vsc5["branch"]
    baseMVA_vsc5 = mpc_vsc5["baseMVA"]
    Ybus_vsc5 = makeYbus(baseMVA_vsc5, bus_vsc5, branch_vsc5)
    Ybus_vsc5 = Matrix(Ybus_vsc5)
    G_vsc5 = real.(Ybus_vsc5)
    B_vsc5 = imag.(Ybus_vsc5)
    nbus_vsc5    = size(bus_vsc5, 1)
    nbranch_vsc5 = size(branch_vsc5, 1)
    pd_vsc5 = bus_vsc2[:, 3] ./ baseMVA_vsc5
    qd_vsc5 = bus_vsc2[:, 4] ./ baseMVA_vsc5
    refBus_vsc5 = findall(bus_vsc5[:, 2] .== 3)

    aloss_vsc = 0.0011
    bloss_vsc = 0.003
    closs_vsc = 0.0046
    Imax_vsc = 1
    sMax_vsc1 = 100 / baseMVA_vsc1
    sMax_vsc2 = 100 / baseMVA_vsc2
    sMax_vsc3 = 100 / baseMVA_vsc3
    sMax_vsc4 = 100 / baseMVA_vsc4
    sMax_vsc5 = 100 / baseMVA_vsc5

    # -----------------------------------------------------------------------------
    # 4.1.2 VSC variables
    # -----------------------------------------------------------------------------
    @variable(model, pn_vsc1[1:nbus_vsc1])
    @variable(model, qn_vsc1[1:nbus_vsc1])
    @variable(model, pgen_vsc1[1:nbus_vsc1])
    @variable(model, qgen_vsc1[1:nbus_vsc1])
    @variable(model, pa2v_vsc1[1:nbus_vsc1])
    @variable(model, qa2v_vsc1[1:nbus_vsc1])
    @variable(model, pl_vsc1[1:nbus_vsc1, 1:nbus_vsc1])
    @variable(model, ql_vsc1[1:nbus_vsc1, 1:nbus_vsc1])
    @variable(model, cc_vsc1[1:nbus_vsc1, 1:nbus_vsc1])
    @variable(model, ss_vsc1[1:nbus_vsc1, 1:nbus_vsc1])
    @variable(model, pv2d_vsc1)
    @variable(model, pv2a_vsc1)
    @variable(model, qv2a_vsc1)
    @variable(model, I_vsc1)
    @variable(model, l_vsc1)
    @variable(model, ploss_vsc1)

    @variable(model, pn_vsc2[1:nbus_vsc2])
    @variable(model, qn_vsc2[1:nbus_vsc2])
    @variable(model, pgen_vsc2[1:nbus_vsc2])
    @variable(model, qgen_vsc2[1:nbus_vsc2])
    @variable(model, pa2v_vsc2[1:nbus_vsc2])
    @variable(model, qa2v_vsc2[1:nbus_vsc2])
    @variable(model, pl_vsc2[1:nbus_vsc2, 1:nbus_vsc2])
    @variable(model, ql_vsc2[1:nbus_vsc2, 1:nbus_vsc2])
    @variable(model, cc_vsc2[1:nbus_vsc2, 1:nbus_vsc2])
    @variable(model, ss_vsc2[1:nbus_vsc2, 1:nbus_vsc2])
    @variable(model, pv2d_vsc2)
    @variable(model, pv2a_vsc2)
    @variable(model, qv2a_vsc2)
    @variable(model, I_vsc2)
    @variable(model, l_vsc2)
    @variable(model, ploss_vsc2)

    @variable(model, pn_vsc3[1:nbus_vsc3])
    @variable(model, qn_vsc3[1:nbus_vsc3])
    @variable(model, pgen_vsc3[1:nbus_vsc3])
    @variable(model, qgen_vsc3[1:nbus_vsc3])
    @variable(model, pa2v_vsc3[1:nbus_vsc3])
    @variable(model, qa2v_vsc3[1:nbus_vsc3])
    @variable(model, pl_vsc3[1:nbus_vsc3, 1:nbus_vsc3])
    @variable(model, ql_vsc3[1:nbus_vsc3, 1:nbus_vsc3])
    @variable(model, cc_vsc3[1:nbus_vsc3, 1:nbus_vsc3])
    @variable(model, ss_vsc3[1:nbus_vsc3, 1:nbus_vsc3])
    @variable(model, pv2d_vsc3)
    @variable(model, pv2a_vsc3)
    @variable(model, qv2a_vsc3)
    @variable(model, I_vsc3)
    @variable(model, l_vsc3)
    @variable(model, ploss_vsc3)

    @variable(model, pn_vsc4[1:nbus_vsc4])
    @variable(model, qn_vsc4[1:nbus_vsc4])
    @variable(model, pgen_vsc4[1:nbus_vsc4])
    @variable(model, qgen_vsc4[1:nbus_vsc4])
    @variable(model, pa2v_vsc4[1:nbus_vsc4])
    @variable(model, qa2v_vsc4[1:nbus_vsc4])
    @variable(model, pl_vsc4[1:nbus_vsc4, 1:nbus_vsc4])
    @variable(model, ql_vsc4[1:nbus_vsc4, 1:nbus_vsc4])
    @variable(model, cc_vsc4[1:nbus_vsc4, 1:nbus_vsc4])
    @variable(model, ss_vsc4[1:nbus_vsc4, 1:nbus_vsc4])
    @variable(model, pv2d_vsc4)
    @variable(model, pv2a_vsc4)
    @variable(model, qv2a_vsc4)
    @variable(model, I_vsc4)
    @variable(model, l_vsc4)
    @variable(model, ploss_vsc4)

    @variable(model, pn_vsc5[1:nbus_vsc5])
    @variable(model, qn_vsc5[1:nbus_vsc5])
    @variable(model, pgen_vsc5[1:nbus_vsc5])
    @variable(model, qgen_vsc5[1:nbus_vsc5])
    @variable(model, pa2v_vsc5[1:nbus_vsc5])
    @variable(model, qa2v_vsc5[1:nbus_vsc5])
    @variable(model, pl_vsc5[1:nbus_vsc5, 1:nbus_vsc5])
    @variable(model, ql_vsc5[1:nbus_vsc5, 1:nbus_vsc5])
    @variable(model, cc_vsc5[1:nbus_vsc5, 1:nbus_vsc5])
    @variable(model, ss_vsc5[1:nbus_vsc5, 1:nbus_vsc5])
    @variable(model, pv2d_vsc5)
    @variable(model, pv2a_vsc5)
    @variable(model, qv2a_vsc5)
    @variable(model, I_vsc5)
    @variable(model, l_vsc5)
    @variable(model, ploss_vsc5)

    # =============================================================================
    # S4.2 VSC constraints
    # =============================================================================

    con_vsc1 = JuMP.ConstraintRef[]
    con_vsc2 = JuMP.ConstraintRef[]
    con_vsc3 = JuMP.ConstraintRef[]
    con_vsc4 = JuMP.ConstraintRef[]
    con_vsc5 = JuMP.ConstraintRef[]

    # -----------------------------------------------------------------------------
    # 4.2.1 Nodal active/reactive power flow equations
    # -----------------------------------------------------------------------------
    add_vsc_nodal_constraints!(model, nbus_vsc1, G_vsc1, B_vsc1, pn_vsc1, qn_vsc1, cc_vsc1, ss_vsc1)
    add_vsc_nodal_constraints!(model, nbus_vsc2, G_vsc2, B_vsc2, pn_vsc2, qn_vsc2, cc_vsc2, ss_vsc2)
    add_vsc_nodal_constraints!(model, nbus_vsc3, G_vsc3, B_vsc3, pn_vsc3, qn_vsc3, cc_vsc3, ss_vsc3)
    add_vsc_nodal_constraints!(model, nbus_vsc4, G_vsc4, B_vsc4, pn_vsc4, qn_vsc4, cc_vsc4, ss_vsc4)
    add_vsc_nodal_constraints!(model, nbus_vsc5, G_vsc5, B_vsc5, pn_vsc5, qn_vsc5, cc_vsc5, ss_vsc5)

    # -----------------------------------------------------------------------------
    # 4.2.2 PCC-side power balance
    # -----------------------------------------------------------------------------
    
    acSideBus_vsc = 1
    dcSideBus_vsc = 2

    @constraint(model, pn_vsc1[acSideBus_vsc] == -pv2a_vsc1)
    @constraint(model, qn_vsc1[dcSideBus_vsc] == -qv2a_vsc1)
    @constraint(model, pn_vsc2[acSideBus_vsc] == -pv2a_vsc2)
    @constraint(model, qn_vsc2[dcSideBus_vsc] == -qv2a_vsc2)
    @constraint(model, pn_vsc3[acSideBus_vsc] == -pv2a_vsc3)
    @constraint(model, qn_vsc3[dcSideBus_vsc] == -qv2a_vsc3)
    @constraint(model, pn_vsc4[acSideBus_vsc] == -pv2a_vsc4)
    @constraint(model, qn_vsc4[dcSideBus_vsc] == -qv2a_vsc4)
    @constraint(model, pn_vsc5[acSideBus_vsc] == -pv2a_vsc5)
    @constraint(model, qn_vsc5[dcSideBus_vsc] == -qv2a_vsc5)

    # -----------------------------------------------------------------------------
    # 4.2.3 Capacity constraint
    # -----------------------------------------------------------------------------
    N = 16
    add_polygon_capacity!(model, pn_vsc1[dcSideBus_vsc], qn_vsc1[dcSideBus_vsc], sMax_vsc1, N)
    add_polygon_capacity!(model, pn_vsc2[dcSideBus_vsc], qn_vsc2[dcSideBus_vsc], sMax_vsc2, N)
    add_polygon_capacity!(model, pn_vsc3[dcSideBus_vsc], qn_vsc3[dcSideBus_vsc], sMax_vsc3, N)
    add_polygon_capacity!(model, pn_vsc4[dcSideBus_vsc], qn_vsc4[dcSideBus_vsc], sMax_vsc4, N)
    add_polygon_capacity!(model, pn_vsc5[dcSideBus_vsc], qn_vsc5[dcSideBus_vsc], sMax_vsc5, N)

    # -----------------------------------------------------------------------------
    # 4.2.4 VSC power loss
    # -----------------------------------------------------------------------------
    @constraint(model, pn_vsc1[dcSideBus_vsc] + ploss_vsc1 + pv2d_vsc1 == 0)
    @constraint(model, ploss_vsc1 == aloss_vsc + bloss_vsc * I_vsc1 + closs_vsc * l_vsc1)
    @constraint(model, 0 <= I_vsc1 <= Imax_vsc)
    @constraint(model, 0 <= l_vsc1 <= Imax_vsc^2)

    @constraint(model, pn_vsc2[dcSideBus_vsc] + ploss_vsc2 + pv2d_vsc2 == 0)
    @constraint(model, ploss_vsc2 == aloss_vsc + bloss_vsc * I_vsc2 + closs_vsc * l_vsc2)
    @constraint(model, 0 <= I_vsc2 <= Imax_vsc)
    @constraint(model, 0 <= l_vsc2 <= Imax_vsc^2)

    @constraint(model, pn_vsc3[dcSideBus_vsc] + ploss_vsc3 + pv2d_vsc3 == 0)
    @constraint(model, ploss_vsc3 == aloss_vsc + bloss_vsc * I_vsc3 + closs_vsc * l_vsc3)
    @constraint(model, 0 <= I_vsc3 <= Imax_vsc)
    @constraint(model, 0 <= l_vsc3 <= Imax_vsc^2)

    @constraint(model, pn_vsc4[dcSideBus_vsc] + ploss_vsc4 + pv2d_vsc4 == 0)
    @constraint(model, ploss_vsc4 == aloss_vsc + bloss_vsc * I_vsc4 + closs_vsc * l_vsc4)
    @constraint(model, 0 <= I_vsc4 <= Imax_vsc)
    @constraint(model, 0 <= l_vsc4 <= Imax_vsc^2)

    @constraint(model, pn_vsc5[dcSideBus_vsc] + ploss_vsc5 + pv2d_vsc5 == 0)
    @constraint(model, ploss_vsc5 == aloss_vsc + bloss_vsc * I_vsc5 + closs_vsc * l_vsc5)
    @constraint(model, 0 <= I_vsc5 <= Imax_vsc)
    @constraint(model, 0 <= l_vsc5 <= Imax_vsc^2)

    # -----------------------------------------------------------------------------
    # 4.2.5 SOCP relaxation
    # -----------------------------------------------------------------------------
    @constraint(model, I_vsc1^2 <= l_vsc1)
    @constraint(model, pn_vsc1[dcSideBus_vsc]^2 + qn_vsc1[dcSideBus_vsc]^2 <= cc_vsc1[dcSideBus_vsc, dcSideBus_vsc] * l_vsc1)
    add_vsc_socp_constraints!(model, nbus_vsc1, nbranch_vsc1, branch_vsc1, cc_vsc1, ss_vsc1)

    @constraint(model, I_vsc2^2 <= l_vsc2)
    @constraint(model, pn_vsc2[dcSideBus_vsc]^2 + qn_vsc2[dcSideBus_vsc]^2 <= cc_vsc2[dcSideBus_vsc, dcSideBus_vsc] * l_vsc2)
    add_vsc_socp_constraints!(model, nbus_vsc2, nbranch_vsc2, branch_vsc2, cc_vsc2, ss_vsc2)

    @constraint(model, I_vsc3^2 <= l_vsc3)
    @constraint(model, pn_vsc3[dcSideBus_vsc]^2 + qn_vsc3[dcSideBus_vsc]^2 <= cc_vsc3[dcSideBus_vsc, dcSideBus_vsc] * l_vsc3)
    add_vsc_socp_constraints!(model, nbus_vsc3, nbranch_vsc3, branch_vsc3, cc_vsc3, ss_vsc3)

    @constraint(model, I_vsc4^2 <= l_vsc4)
    @constraint(model, pn_vsc4[dcSideBus_vsc]^2 + qn_vsc4[dcSideBus_vsc]^2 <= cc_vsc4[dcSideBus_vsc, dcSideBus_vsc] * l_vsc4)
    add_vsc_socp_constraints!(model, nbus_vsc4, nbranch_vsc4, branch_vsc4, cc_vsc4, ss_vsc4)

    @constraint(model, I_vsc5^2 <= l_vsc5)
    @constraint(model, pn_vsc5[dcSideBus_vsc]^2 + qn_vsc5[dcSideBus_vsc]^2 <= cc_vsc5[dcSideBus_vsc, dcSideBus_vsc] * l_vsc5)
    add_vsc_socp_constraints!(model, nbus_vsc5, nbranch_vsc5, branch_vsc5, cc_vsc5, ss_vsc5)

    # -----------------------------------------------------------------------------
    # 4.2.6 Objective term
    # -----------------------------------------------------------------------------
    obj_vsc = @expression(model,
        sum(pn_vsc1[i]^2 for i in 1:nbus_vsc1) +
        sum(pn_vsc2[i]^2 for i in 1:nbus_vsc2) +
        sum(pn_vsc3[i]^2 for i in 1:nbus_vsc3) +
        sum(pn_vsc4[i]^2 for i in 1:nbus_vsc4) +
        sum(pn_vsc5[i]^2 for i in 1:nbus_vsc5)
    )

    # =============================================================================
    # S5.1 Load Communication network parameters & variables
    # =============================================================================

    mpc_info = parse_mpc(joinpath(@__DIR__, "case5_info.m"))
    bus_info     = mpc_info["bus"]
    gen_info     = mpc_info["gen"]
    branch_info  = mpc_info["branch"]

    nbus_info    = size(bus_info, 1)
    ngen_info    = size(gen_info, 1)
    nbranch_info = size(branch_info, 1)

    refBus_info = findall(bus_info[:, 2] .== 3)
    fbus_info   = Int.(branch_info[:, 1])
    tbus_info   = Int.(branch_info[:, 2])

    dij_info = zeros(nbus_info, nbus_info)
    for kk in 1:nbranch_info
        i = fbus_info[kk]
        j = tbus_info[kk]
        dij_info[i, j] = branch_info[kk, 3]
        dij_info[j, i] = dij_info[i, j] - 5
    end

    @variable(model, s_info[1:nbus_info])
    @variable(model, d_info[1:nbus_info])
    @variable(model, fij_info[1:nbus_info, 1:nbus_info])
    @variable(model, fposij_info[1:nbus_info, 1:nbus_info])
    @variable(model, gama_info[1:nbus_info, 1:nbus_info], Bin)

    # =============================================================================
    # S5.2 Communication network constraints
    # =============================================================================

    # -----------------------------------------------------------------------------
    # 5.2.1 Information balance
    # -----------------------------------------------------------------------------
    for k in 1:nbus_info
        idxFrom_info = findall(fbus_info .== k)
        idxTo_info   = findall(tbus_info .== k)

        if k in refBus_info
            @constraint(model,
                s_info[k] == sum(fij_info[k, tbus_info[ii]] for ii in idxFrom_info) +
                             sum(fij_info[k, fbus_info[ii]] for ii in idxTo_info)
            )
        else
            @constraint(model,
                d_info[k] == -(
                    sum(fij_info[k, tbus_info[ii]] for ii in idxFrom_info) +
                    sum(fij_info[k, fbus_info[ii]] for ii in idxTo_info)
                )
            )
            @constraint(model,
                d_info[k] == -(
                    sum(dij_info[k, tbus_info[ii]] * gama_info[k, tbus_info[ii]] for ii in idxFrom_info) +
                    sum(dij_info[k, fbus_info[ii]] * gama_info[k, fbus_info[ii]] for ii in idxTo_info)
                )
            )
        end
    end

    for kk in 1:nbranch_info
        i = fbus_info[kk]
        j = tbus_info[kk]
        @constraint(model, fij_info[i, j] == -fij_info[j, i])
    end

    # -----------------------------------------------------------------------------
    # 5.2.2 Coupling with power-flow switching decisions
    # -----------------------------------------------------------------------------
    lam_info = zeros(nbus_info, nbus_info)

    betak_dc[1, 2] = 1; lam_info[1, 2] = 15
    betak_dc[2, 1] = 1; lam_info[2, 1] = 15
    betak_dc[1, 4] = 1; lam_info[1, 4] = 25
    betak_dc[4, 1] = 1; lam_info[4, 1] = 30
    betak_dc[1, 5] = 0; lam_info[1, 5] = 1e6  # means info link #1-#5 break
    betak_dc[5, 1] = 0; lam_info[5, 1] = 1e6
    betak_dc[2, 3] = 1; lam_info[2, 3] = 20
    betak_dc[3, 2] = 1; lam_info[3, 2] = 20
    betak_dc[2, 4] = 0; lam_info[2, 4] = 1e6
    betak_dc[4, 2] = 0; lam_info[4, 2] = 1e6
    betak_dc[2, 5] = 0; lam_info[2, 5] = 15
    betak_dc[5, 2] = 0; lam_info[5, 2] = 10
    betak_dc[3, 4] = 0; lam_info[3, 4] = 1e6
    betak_dc[4, 3] = 0; lam_info[4, 3] = 1e6
    betak_dc[3, 5] = 1; lam_info[3, 5] = 25
    betak_dc[5, 3] = 1; lam_info[5, 3] = 30
    betak_dc[4, 5] = 1; lam_info[4, 5] = 40
    betak_dc[5, 4] = 1; lam_info[5, 4] = 40

    for kk in 1:nbranch_dc
        i = fbus_dc[kk]
        j = tbus_dc[kk]
        add_couple_logic_constraints!(model, gama_info[i, j], beta_dc[i, j], betak_dc[i, j])
        add_couple_logic_constraints!(model, gama_info[j, i], beta_dc[j, i], betak_dc[j, i])
    end

    # -----------------------------------------------------------------------------
    # 5.2.3 Communication objective
    # -----------------------------------------------------------------------------
    for kk in 1:nbranch_info
        i = fbus_info[kk]
        j = tbus_info[kk]
        @constraint(model, fposij_info[i, j] >= fij_info[i, j])
        @constraint(model, fposij_info[i, j] >= 0)
        @constraint(model, fposij_info[j, i] >= fij_info[j, i])
        @constraint(model, fposij_info[j, i] >= 0)
    end

    obj_info = @expression(model,
        sum(lam_info[fbus_info[kk], tbus_info[kk]] * fposij_info[fbus_info[kk], tbus_info[kk]] +
            lam_info[tbus_info[kk], fbus_info[kk]] * fposij_info[tbus_info[kk], fbus_info[kk]]
            for kk in 1:nbranch_info)
    )

    # =============================================================================
    # S6 Coupling constraints
    # =============================================================================

    @constraint(model, pa2v_ac[4] + pv2a_vsc1 == 0)
    @constraint(model, qa2v_ac[4] + qv2a_vsc1 == 0)
    @constraint(model, cc_ac[4, 4] == cc_vsc1[1, 1])

    @constraint(model, pa2v_ac[6] + pv2a_vsc2 == 0)
    @constraint(model, qa2v_ac[6] + qv2a_vsc2 == 0)
    @constraint(model, cc_ac[6, 6] == cc_vsc2[1, 1])

    @constraint(model, pa2v_ac[8] + pv2a_vsc3 == 0)
    @constraint(model, qa2v_ac[8] + qv2a_vsc3 == 0)
    @constraint(model, cc_ac[8, 8] == cc_vsc3[1, 1])

    @constraint(model, pr2v_res1 + pv2a_vsc4 == 0)
    @constraint(model, qr2v_res1 + qv2a_vsc4 == 0)

    @constraint(model, pr2v_res2 + pv2a_vsc5 == 0)
    @constraint(model, qr2v_res2 + qv2a_vsc5 == 0)

    @constraint(model, pn_dc[1] - pv2d_vsc1 == 0)
    @constraint(model, pn_dc[2] - pv2d_vsc2 == 0)
    @constraint(model, pn_dc[3] - pv2d_vsc3 == 0)
    @constraint(model, pn_dc[4] - pv2d_vsc4 == 0)
    @constraint(model, pn_dc[5] - pv2d_vsc5 == 0)

    # =============================================================================
    # S7 Objective and solve
    # =============================================================================

    obj = @expression(model, obj_ac + obj_dc + obj_info)
    @objective(model, Min, obj)

    optimize!(model)

    term_status = termination_status(model)
    primal_status_val = primal_status(model)

    if term_status != MOI.OPTIMAL && term_status != MOI.LOCALLY_SOLVED
        @warn "Optimization did not report an optimal solution." term_status primal_status_val
    end

    P_sym = [
        (abs(value(pij_dc[i, j])) + abs(value(pij_dc[j, i]))) / 2
        for i in 1:nbus_dc, j in 1:nbus_dc
    ]

    println("Termination status: ", term_status)
    println("Objective value: ", objective_value(model))

    println("\n================ beta_dc matrix =================")

    beta_val = round.(Int, value.(beta_dc))

    display(beta_val)
    
    println("\n================ P_sym matrix (MW) =================")
    display(P_sym .* baseMW_dc)

    println("\n================ fij_info matrix =================")
    display(value.(fij_info))

    println("\n================ Communication Flow =================")
    @printf("%-10s %-10s %14s %14s\n",
        "From Node", "To Node", "Fij(MB/s)", "Fji(MB/s)")

    for kk in 1:nbranch_info
        i = fbus_info[kk]
        j = tbus_info[kk]

        lam_ij = lam_info[i, j]
        lam_ji = lam_info[j, i]
        if lam_ij >= 1e6 || lam_ji >= 1e6
            continue
        end

        Fij = value(fij_info[i, j])
        Fji = value(fij_info[j, i])

        @printf("%-10d %-10d %14.2f %14.2f\n",
            i, j, Fij, Fji)
    end

    println("\n================ DC Breaker Status =================")
    @printf("%-10s %-10s %-14s %-14s\n", "From Bus", "To Bus", "Breaker@From", "Breaker@To")

    for kk in 1:nbranch_dc
        i = fbus_dc[kk]
        j = tbus_dc[kk]

        status_i = round(Int, value(beta_dc[i,j])) == 1 ? "Closed" : "Open"
        status_j = round(Int, value(beta_dc[j,i])) == 1 ? "Closed" : "Open"

        @printf("%-10d %-10d %-14s %-14s\n", i, j, status_i, status_j)
    end

    println("\n================ DC Branch Power Flow =================")
    @printf("%-10s %-10s %12s %12s %12s\n", "From Bus", "To Bus", "Pij(MW)", "Pji(MW)", "Psym(MW)")

    for kk in 1:nbranch_dc
        i = fbus_dc[kk]
        j = tbus_dc[kk]

        Pij  = value(pij_dc[i,j]) * baseMW_dc
        Pji  = value(pij_dc[j,i]) * baseMW_dc
        Psym = P_sym[i,j] * baseMW_dc

        @printf("%-10d %-10d %12.2f %12.2f %12.2f\n", i, j, Pij, Pji, Psym)
    end
end

main()
