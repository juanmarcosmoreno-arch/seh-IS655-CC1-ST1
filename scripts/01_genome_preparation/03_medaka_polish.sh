#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MEDAKA_ENV="$PROJECT_DIR/.micromamba/envs/gridion_medaka"
METADATA="$PROJECT_DIR/metadata/barcode_muestras.tsv"
THREADS="${THREADS:-4}"
MODEL="${MODEL:-r1041_e82_400bps_sup_v5.2.0}"
ONLY_BARCODE="${ONLY_BARCODE:-}"

export PATH="$MEDAKA_ENV/bin:$PATH"

mkdir -p "$PROJECT_DIR/results/medaka" "$PROJECT_DIR/logs"

tail -n +2 "$METADATA" | while IFS=$'\t' read -r barcode_num barcode_id muestra_id fastq estado notas; do
  if [[ "$estado" != "asignado" ]]; then
    echo "Saltando $barcode_id: $estado"
    continue
  fi

  if [[ -n "$ONLY_BARCODE" && "$barcode_id" != "$ONLY_BARCODE" ]]; then
    continue
  fi

  sample="cepa_${muestra_id}_${barcode_id}"
  reads="$PROJECT_DIR/results/filtered_fastq/${sample}.filtered.fastq.gz"
  outdir="$PROJECT_DIR/results/medaka/$sample"

  # Usar RagTag si existe, si no caer a Flye directo
  ragtag_assembly="$PROJECT_DIR/results/ragtag/$sample/ragtag.scaffold.fasta"
  flye_assembly="$PROJECT_DIR/results/flye/$sample/assembly.fasta"
  if [[ -s "$ragtag_assembly" ]]; then
    assembly="$ragtag_assembly"
    echo "Usando ensamblaje RagTag para $sample"
  else
    assembly="$flye_assembly"
    echo "RagTag no encontrado, usando Flye directo para $sample"
  fi

  if [[ ! -s "$reads" ]]; then
    echo "No existe input filtrado para $sample: $reads" >&2
    exit 1
  fi

  if [[ ! -s "$assembly" ]]; then
    echo "No existe ensamblaje para $sample: $assembly" >&2
    exit 1
  fi

  echo "Puliendo $sample"
  medaka_consensus \
    -i "$reads" \
    -d "$assembly" \
    -o "$outdir" \
    -t "$THREADS" \
    -m "$MODEL" \
    > "$PROJECT_DIR/logs/${sample}.medaka.log" 2>&1
done
