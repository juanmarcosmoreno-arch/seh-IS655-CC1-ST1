#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NANOPORE_ENV="$PROJECT_DIR/.micromamba/envs/gridion_nanopore"
METADATA="$PROJECT_DIR/metadata/barcode_muestras.tsv"
MIN_LENGTH="${MIN_LENGTH:-1000}"
KEEP_PERCENT="${KEEP_PERCENT:-90}"
TARGET_BASES="${TARGET_BASES:-500000000}"
ONLY_BARCODE="${ONLY_BARCODE:-}"

export PATH="$NANOPORE_ENV/bin:$PATH"

mkdir -p "$PROJECT_DIR/results/filtered_fastq" "$PROJECT_DIR/logs"

resolve_fastq() {
  local requested="$1"
  local direct="$PROJECT_DIR/$requested"
  local gz="${direct}.gz"

  if [[ -s "$direct" ]]; then
    printf '%s\n' "$direct"
  elif [[ -s "$gz" ]]; then
    printf '%s\n' "$gz"
  else
    return 1
  fi
}

tail -n +2 "$METADATA" | while IFS=$'\t' read -r barcode_num barcode_id muestra_id fastq estado notas; do
  if [[ "$estado" != "asignado" ]]; then
    echo "Saltando $barcode_id: $estado"
    continue
  fi

  if [[ -n "$ONLY_BARCODE" && "$barcode_id" != "$ONLY_BARCODE" ]]; then
    continue
  fi

  input="$(resolve_fastq "$fastq")" || {
    echo "No existe FASTQ para $barcode_id: $PROJECT_DIR/$fastq o $PROJECT_DIR/$fastq.gz" >&2
    exit 1
  }
  sample="cepa_${muestra_id}_${barcode_id}"
  output="$PROJECT_DIR/results/filtered_fastq/${sample}.filtered.fastq.gz"

  echo "Filtrando $sample"
  filtlong \
    --min_length "$MIN_LENGTH" \
    --keep_percent "$KEEP_PERCENT" \
    --target_bases "$TARGET_BASES" \
    "$input" \
    2> "$PROJECT_DIR/logs/${sample}.filtlong.log" \
    | gzip > "$output"
done
