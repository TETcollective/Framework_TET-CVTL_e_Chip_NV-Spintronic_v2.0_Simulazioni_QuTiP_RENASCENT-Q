# Code Folder - RENASCENT-Q Proxy Simulation

Tutti i codici di simulazione sviluppati per questo lavoro sono contenuti nella cartella `code/` del repository del TET Collective. Gli script sono modulari, ben commentati e direttamente collegati alle sezioni teoriche e numeriche del paper.

## File principali

- `system_setup.py`  
  → Definizione dei parametri del sistema (standard e modified), Hamiltoniana, operatori dissipativi e condizioni iniziali.

- `retro_negentropic_term.py`  
  → Implementazione del termine negentropico custom \(\mathcal{L}_{\text{neg}}(\beta)\) parametrizzato dal rapporto aureo inverso \(\beta = \phi^{-2}\), con approccio anti-decoerenza mirata basato sul weak value.

- `floquet_d_beta_qobevo.py`  
  → Costruzione del Liouvillian completo con QobjEvo, inclusa la modulazione periodica SAW (Surface Acoustic Waves) e il termine dissipativo retrocausale.

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

