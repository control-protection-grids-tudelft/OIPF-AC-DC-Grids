# OIPF for Transmission Switching in AC/MTDC Grids

<table width="100%">
  <tr>
    <td align="left"><img src="IEPG_logo.jpg" alt="Intelligent Electrical Power Grids (IEPG)" height="80"></td>
    <td align="right"><img src="cresym.png" alt="CRESYM" height="80"></td>
  </tr>
</table>

This repository contains the source code accompanying the paper
**"Leveraging Optimal Information-Power Flow for Transmission Switching in AC/MTDC Grids"**, presented at the **26th International Conference on Environment and Electrical Engineering (IEEE EEEIC 2026)**.

## Paper and presentation

- [Paper on arXiv](https://arxiv.org/abs/2606.09436)
- [EEEIC 2026 presentation slides (PDF)](h.li_169_eeeic2026.pdf)
- arXiv DOI: [10.48550/arXiv.2606.09436](https://doi.org/10.48550/arXiv.2606.09436)
- Zenodo DOI: pending archival release

## Overview

This code contains optimization models for studying transmission switching in AC/multi-terminal DC grids. The main idea is to compare conventional optimal power flow (OPF) with optimal information-power flow (OIPF), where communication network constraints and costs are considered together with power flow switching decisions.

## Comparison Experiments

Run `nlp_socp_ac_dc_compare/acdc_opf_nlp.m` and `nlp_socp_ac_dc_compare/acdc_opf_socp.m` to compare the nonlinear AC/DC OPF model with its SOCP-relaxed formulation. The outputs show the optimization status, total generation cost, solver time, and AC generator active/reactive power outputs, which can be used to observe the accuracy of the SOCP-relaxed AC/DC OPF model.

Run `opf_oipf_compare/acdc_opf.jl` and `opf_oipf_compare/acdc_iopf.jl` to compare conventional OPF with proposed OIPF. The outputs show the objective value, DC breaker switching status, DC branch power flows, and, for OIPF, the communication link flows. This comparison shows how information-network constraints and costs affect transmission-switching decisions.

## Citation

If you use this code, please cite:

```bibtex
@misc{li2026oipf,
  title         = {Leveraging Optimal Information-Power Flow for Transmission Switching in AC/MTDC Grids},
  author        = {Li, Haixiao and Leki\'{c}, Aleksandra},
  year          = {2026},
  eprint        = {2606.09436},
  archivePrefix = {arXiv},
  primaryClass  = {eess.SY},
  doi           = {10.48550/arXiv.2606.09436},
  url           = {https://arxiv.org/abs/2606.09436}
}
```
