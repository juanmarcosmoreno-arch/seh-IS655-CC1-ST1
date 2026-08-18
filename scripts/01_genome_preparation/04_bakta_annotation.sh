#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAKTA_ENV="$PROJECT_DIR/.micromamba/envs/gridion_bakta"
METADATA="$PROJECT_DIR/metadata/barcode_muestras.tsv"
THREADS="${THREADS:-4}"
ONLY_BARCODE="${ONLY_BARCODE:-}"
BAKTA_DB="${BAKTA_DB:-$PROJECT_DIR/databases/bakta_db}"

export PATH="$BAKTA_ENV/bin:$PATH"
export BLAST_USAGE_REPORT=false
export MPLCONFIGDIR="$PROJECT_DIR/.cache/matplotlib"
export XDG_CACHE_HOME="$PROJECT_DIR/.cache"

mkdir -p "$PROJECT_DIR/results/bakta" "$PROJECT_DIR/logs" "$MPLCONFIGDIR" "$XDG_CACHE_HOME/fontconfig"

if [[ ! -d "$BAKTA_DB" ]]; then
  echo "Base de datos Bakta no encontrada en: $BAKTA_DB" >&2
  echo "Descargar con:" >&2
  echo "  mkdir -p $BAKTA_DB && bakta_db download --output $BAKTA_DB --type light" >&2
  exit 1
fi

tail -n +2 "$METADATA" | while IFS=$'\t' read -r barcode_num barcode_id muestra_id fastq estado notas; do
  if [[ "$estado" != "asignado" ]]; then
    echo "Saltando $barcode_id: $estado"
    continue
  fi

  if [[ -n "$ONLY_BARCODE" && "$barcode_id" != "$ONLY_BARCODE" ]]; then
    continue
  fi

  sample="cepa_${muestra_id}_${barcode_id}"
  assembly="$PROJECT_DIR/results/medaka/$sample/consensus.fasta"
  output_dir="$PROJECT_DIR/results/bakta/$sample"

  if [[ ! -s "$assembly" ]]; then
    echo "No existe ensamblaje Medaka para $sample: $assembly" >&2
    exit 1
  fi

  if [[ -s "$output_dir/${sample}.gff3" && -s "$output_dir/${sample}.tsv" && -s "$output_dir/${sample}.json" ]]; then
    echo "Bakta ya existe, saltando: $sample"
    continue
  fi

  echo "Anotando $sample con Bakta"
  bakta \
    --db "$BAKTA_DB" \
    --output "$output_dir" \
    --prefix "$sample" \
    --threads "$THREADS" \
    --genus Staphylococcus \
    --skip-ori \
    --skip-plot \
    --force \
    "$assembly" \
    > "$PROJECT_DIR/logs/${sample}.bakta.log" 2>&1
done
