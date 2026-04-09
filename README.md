# Framework TET-CVTL and NV-Spintronic Chip v2.0

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19485796.svg)](https://doi.org/10.5281/zenodo.19485796)
[![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

**Stabilizzazione Negentropica Retrocausale tramite Braiding Anyonico e termine β = 1.61803^{-2}**  
**QuTiP Simulations of the RENASCENT-Q Protocol with Renascent Entanglement and Weak Value Modulation**

**Autore:** Simon Soliman  
**Visual Artist & Independent Researcher** — TET Collective, Roma, Italy  
**ORCID:** 0009-0002-3533-3772  
**Sito:** [tetcollective.org](https://tetcollective.org)  
**Email:** tetcollective@proton.me  

**Data:** Aprile 2026  
**Zenodo DOI:** [10.5281/zenodo.19485796](https://doi.org/10.5281/zenodo.19485796)

---

## Abstract

All’interno del **Framework TET-CVTL**, questo lavoro presenta un’estensione retrocausale negentropica della master equation Floquet-Lindblad per sistemi quantistici aperti, finalizzata alla stabilizzazione coerente a temperatura ambiente (T ≈ 300 K).

Viene introdotto un termine di impulso custom parametrizzato dal **rapporto aureo inverso β = 1.61803^{-2} ≈ 0.382**, ispirato a un’estensione embodied della teoria Orchestrated Objective Reduction (Orch-OR) di Penrose-Hameroff e alle condizioni al contorno del Two-State Vector Formalism (TSVF).  

Utilizzando simulazioni ad alta fedeltà con il framework open-source **QuTiP** su proxy realistici a 2 e 4 qubit (rappresentativi di dimeri di tubulina nei microtubuli o centri NV in diamante) soggetti a modulazione strain periodica tramite Surface Acoustic Waves (SAW), si dimostra una significativa riduzione dell’entropia di von Neumann, un marcato incremento della concurrence e un prolungamento sostanziale del tempo di coerenza, ben oltre i limiti dissipativi standard.

I risultati rafforzano la fattibilità di reti anyoniche stabili all’interno del **Chip NV-Spintronic TET-CVTL v2.0** e del protocollo **RENASCENT-Q**.

---

## Contenuti del Repository

- `Framework_TET_CVTL_e_Chip_NV_Spintronic_v2_0___Simulazioni_QuTiP_del_Protocollo_RENASCENT_Q.pdf` → Documento principale (89 pagine)
- File sorgente LaTeX (`.tex`, `.bib`)
- Tutti gli script Python per le simulazioni QuTiP
- Figure e grafici delle simulazioni
- `CITATION.cff`, `LICENSE.md`, `requirements.txt`

---

## How to Cite

**APA Style**  
Soliman, S. (2026). *Framework TET-CVTL and NV-Spintronic Chip v2.0: Negentropic Retrocausal Stabilization via Anyonic Braiding and the β = 1.61803^{-2} Term*. Zenodo. https://doi.org/10.5281/zenodo.19485796

**BibTeX**
```bibtex
@misc{soliman2026framework,
  author       = {Soliman, Simon},
  title        = {Framework TET-CVTL and NV-Spintronic Chip v2.0: Negentropic Retrocausal Stabilization via Anyonic Braiding and the β = 1.61803^{-2} Term. QuTiP Simulations of the RENASCENT-Q Protocol with Renascent Entanglement and Weak Value Modulation},
  year         = {2026},
  month        = {March},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.19485796},
  url          = {https://doi.org/10.5281/zenodo.19485796},
  note         = {ORCID: 0009-0002-3533-3772}
}
- `renascent_q_measurements.py`  
  → Funzioni per il calcolo di entropia di von Neumann, purezza, concurrence, logarithmic negativity e weak value \(\langle A \rangle_w\).

- `protocollo_tsvf_4qubit_anyonic_braiding.py`  
  → Simulazione completa a 4 qubit con braiding anyonico esplicito (statistiche non-Abeliane, fase \(\theta = \pi/4\)), weak value con imprint topologico e protezione modulata da \(\beta = \phi^{-2}\). Produce i risultati di concurrence media finale 0.5360 e oscillazioni del weak value presentati nel paper.

- `run_simulation.py`  
  → Script principale per lanciare le simulazioni (proxy a 2 qubit e modello a 4 qubit).

- `plot_entropy_vs_time_beta.py`  
  → Genera le figure principali del paper (entropia vs tempo per diversi \(\beta\), concurrence renascent e weak value evolution).

## Come eseguire

```bash
# Esecuzione della simulazione principale (proxy + modello a 4 qubit)
python run_simulation.py

