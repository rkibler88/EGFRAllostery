# Structural panels

PyMOL sessions for the structural panels of Figures 1, 2 and 5, and the scripts that write
them. These are the author's originals; they are not part of `make figures` because they need
PyMOL and the coordinate files rather than the committed tables.

| file | what it does |
|---|---|
| `write_alphafold_pml.py` | reads the community table, emits `viz_alphafold_communities.pml` |
| `write_hubs_bridges_pml.py` | emits `viz_hubs_and_bridges.pml` (top five by node strength and by betweenness) |
| `write_1m17_pml.py` | emits `viz_1m17_comparison.pml`, applying the +24 offset between 1M17 and UniProt P00533 numbering |
| `viz_*.pml` | the generated sessions, as used |
| `map_rmsf_to_bfactor.py` | writes per-residue RMSF into the B-factor column of the AlphaFold model for the rigidity map |
| `pymol_switch_handshake.py` | renders the switch/handshake panel on 1M17 (erlotinib-bound) |

Structures: AlphaFold2 model of UniProt P00533 and PDB 1M17. Superposition uses CEALIGN on
the C-lobe scaffold, as described in Methods.
