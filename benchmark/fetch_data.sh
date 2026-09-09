#!/bin/sh
# Fetch the benchmark inputs. Neither is vendored: the labels are third-party
# CC-licensed data and the trajectories are ~700 MB.
set -e
cd "$(dirname "$0")"

# CryptoBench cryptic-site labels (Skrhak et al., Bioinformatics 2025)
curl -sSL -o cb_labels.json \
  https://raw.githubusercontent.com/skrhakv/CryptoBench/master/src/F-statistics/conservation/label_dataset.json

# The 14 CryptoBench entries that have ATLAS trajectories. Regenerate this list
# by probing every CryptoBench entry against the ATLAS API on the chain that
# carries its labels; it is a lower bound, since only that chain is probed.
mkdir -p atlas_cohort && cd atlas_cohort
for t in 1kx9_A 1esw_A 1hp1_A 2fp1_A 1y6i_A 2h7g_X 2jlq_A 2po4_A \
         3b49_A 3ikw_A 3vjz_A 4uc8_A 5op0_B 6irx_A; do
  [ -f "$t.zip" ] || curl -sSL --max-time 900 -o "$t.zip" \
    "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/$t"
  echo "$t"
done
