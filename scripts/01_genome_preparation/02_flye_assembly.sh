#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NANOPORE_ENV="$PROJECT_DIR/.micromamba/envs/gridion_nanopore"
METADATA="$PROJECT_DIR/metadata/barcode_muestras.tsv"
THREADS="${THREADS:-4}"
ONLY_BARCODE="${ONLY_BARCODE:-}"

export PATH="$NANOPORE_ENV/bin:$PATH"

mkdir -p "$PROJECT_DIR/results/flye" "$PROJECT_DIR/logs"

tail -n +2 "$METADATA" | while IFS=$'\t' read -r barcode_num barcode_id muestra_id fastq estado notas; do
  if [[ "$estado" != "asignado" ]]; then
    echo "Saltando $barcode_id: $estado"
    continue
  fi

  if [[ -n "$ONLY_BARCODE" && "$barcode_id" != "$ONLY_BARCODE" ]]; then
    continue
  fi

  sample="cepa_${muestra_id}_${barcode_id}"
  input="$PROJECT_DIR/results/filtered_fastq/${sample}.filtered.fastq.gz"
  outdir="$PROJECT_DIR/results/flye/$sample"

  if [[ ! -s "$input" ]]; then
    echo "No existe input filtrado para $sample: $input" >&2
    exit 1
  fi

  echo "Ensamblando $sample"
  flye \
    --nano-raw "$input" \
    --out-dir "$outdir" \
    --threads "$THREADS" \
    > "$PROJECT_DIR/logs/${sample}.flye.log" 2>&1
done
