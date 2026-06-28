# EEEIC2026

Code for the paper [Leveraging Optimal Information-Power Flow for Transmission Switching in AC/MTDC Grids](https://arxiv.org/abs/2606.09436), which will be presented in 26th International Conference on Environment and Electrical
Engineering (EEEIC2026)

This code contains optimization models for studying transmission switching in AC/multi-terminal DC grids. The main idea is to compare conventional optimal power flow (OPF) with optimal information-power flow (OIPF), where communication network constraints and costs are considered together with power flow switching decisions.

## Comparison Experiments

Run `nlp_socp_ac_dc_compare/acdc_opf_nlp.m` and `nlp_socp_ac_dc_compare/acdc_opf_socp.m` to compare the nonlinear AC/DC OPF model with its SOCP-relaxed formulation. The outputs show the optimization status, total generation cost, solver time, and AC generator active/reactive power outputs, which can be used to observe the accuracy of the SOCP-relaxed AC/DC OPF model.

Run `opf_oipf_compare/acdc_opf.jl` and `opf_oipf_compare/acdc_iopf.jl` to compare conventional OPF with proposed OIPF. The outputs show the objective value, DC breaker switching status, DC branch power flows, and, for OIPF, the communication link flows. This comparison shows how information-network constraints and costs affect transmission-switching decisions.
