"""cofold-bench: Boltz-2 vs AlphaFold2-Multimer on post-cutoff protein heterodimers.

Sub-modules
-----------
config       load ``config.yaml`` and resolve repo-relative paths
rcsb         RCSB Search / Data / Files API helpers
select_targets  CLI: query RCSB, filter, verify with Biotite, write targets.csv
download     CLI: fetch mmCIF references and canonical FASTA sequences
make_inputs  CLI: write Boltz-2 YAML and ColabFold FASTA inputs
structure    Biotite structure utilities (chains, CB coordinates, contacts)
pdockq       pDockQ (Bryant et al. 2022)
confidence   parse Boltz-2 / ColabFold confidence JSON files
dockq_runner thin wrapper over the DockQ package
score        CLI: aggregate everything into results/results.csv
plot         CLI: figures
"""

__version__ = "0.1.0"
