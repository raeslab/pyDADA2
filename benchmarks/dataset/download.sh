#!/usr/bin/env bash
# Fetch the phase-11 benchmark data (spec §9.2) into benchmarks/dataset/data/ (git-ignored).
#
#   * MiSeq SOP (Kozich et al. 2013): 20 samples (19 mice + mock), 16S V4, MiSeq 2x250
#   * SILVA 138.2 DADA2-formatted training set + species set (Zenodo 14169026)
#
# Every file is SHA-256 checked. Usage: benchmarks/dataset/download.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
data="$here/data"
mkdir -p "$data"
cd "$data"

fetch() {  # url file sha256
  if [[ ! -f "$2" ]] || ! echo "$3  $2" | sha256sum -c --quiet - 2>/dev/null; then
    echo "downloading $2"
    curl -fL --retry 3 -o "$2.part" "$1"
    mv "$2.part" "$2"
  fi
  echo "$3  $2" | sha256sum -c -
}

fetch https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip miseqsopdata.zip \
  728f7544824f8ca0223618b18bb60fce6a5dc20ae3d476267dbd5cf5594a35d4
fetch https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz \
  silva_nr99_v138.2_toGenus_trainset.fa.gz e0af8246b900a922ea3ea68e7cb1a5e70ceed4db84246b6ea9865eb132452093
fetch https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz \
  silva_v138.2_assignSpecies.fa.gz 45c7b0d2197e481f60c3bc3d53882f7bc0f3ad3a947248cd14780c3a27a305e6

if [[ ! -d MiSeq_SOP ]]; then
  unzip -q miseqsopdata.zip 'MiSeq_SOP/*.fastq' -d .
  rm -rf __MACOSX
fi
echo "FASTQ files: $(ls MiSeq_SOP/*_R1_001.fastq | wc -l) forward, $(ls MiSeq_SOP/*_R2_001.fastq | wc -l) reverse"
