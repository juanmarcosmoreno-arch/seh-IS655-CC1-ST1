#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd -- "$SCRIPT_DIR/.." && pwd)}"
cd "$PROJECT_DIR"

export PATH="$PROJECT_DIR/.micromamba/envs/gridion_pangenome_pip/bin:$PATH"

SET_DIR="results/seh_st1_bovine/phylogenomics_set"
PANAROO_OUT="results/seh_st1_bovine/phylogenomics_set/panaroo_core"
TREE_PREFIX="results/seh_st1_bovine/phylogenomics_set/iqtree_core/seh_st1_core"

mkdir -p "$(dirname "$TREE_PREFIX")"

command -v panaroo >/dev/null 2>&1 || {
  echo "ERROR: panaroo is not available in PATH"
  exit 1
}

if command -v iqtree2 >/dev/null 2>&1
then
  IQTREE=iqtree2
elif command -v iqtree3 >/dev/null 2>&1
then
  IQTREE=iqtree3
elif command -v iqtree >/dev/null 2>&1
then
  IQTREE=iqtree
else
  echo "ERROR: iqtree2/iqtree is not available in PATH"
  exit 1
fi

python3 scripts/68_build_seh_phylogenomics_set.py

GFF_COUNT=$(find "$SET_DIR/gff3" -name '*.gff3' | wc -l | tr -d ' ')
echo "GFF3 files available for Panaroo: ${GFF_COUNT}"

panaroo \
  -i "$SET_DIR"/gff3/*.gff3 \
  -o "$PANAROO_OUT" \
  --clean-mode strict \
  --remove-invalid-genes \
  -a core \
  --aligner mafft \
  -t 8

"$IQTREE" \
  -s "$PANAROO_OUT/core_gene_alignment.aln" \
  -m MFP \
  -bb 1000 \
  -alrt 1000 \
  -nt AUTO \
  --prefix "$TREE_PREFIX"
