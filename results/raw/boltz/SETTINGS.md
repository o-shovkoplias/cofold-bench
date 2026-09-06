# Boltz-2 run settings (cofold-bench, RTX 4060 Laptop 8 GB, --no_kernels)

Uniform benchmark run (this directory): `--recycling_steps 3 --diffusion_samples 5 --sampling_steps 200 --max_parallel_samples 1
--max_msa_seqs 512 --num_subsampled_msa 128 --no_kernels --seed 42 --use_msa_server --msa_pairing_strategy greedy`.
Reason: with the Boltz defaults (8192 / 1024 MSA sequences) 15 of 22 targets ran out of GPU memory on 8 GB; the memory of the MSA
module scales with the number of loaded MSA sequences, and 512/128 peaked at 3.7 GiB (3819 MiB, nvidia-smi) on 9ASS (326 residues).
Targets that still fail at 512/128 are retried with `--max_msa_seqs 256 --num_subsampled_msa 64` (listed in FALLBACK.txt if any).

Side comparison (`../boltz_default/`): the 7 targets that completed with Boltz defaults (8192 / 1024) before the memory limit was
understood: 8CM0 9E9G 9FWR 9HMX 9J9W 9O24 9Q0B. Used only to estimate the effect of MSA truncation on ipTM/DockQ.
