# Security and sensitive-information audit

Automated scan of all files in this repository for potentially sensitive content.

Scan date: 2026-08-18

---

## Findings

### 1. Local absolute paths — REVIEW REQUIRED

**Files affected:** Multiple scripts in `scripts/01_genome_preparation/`, `scripts/02_typing_and_identity/`, and analysis scripts.

**Pattern:** `/Volumes/JM/BioInf/GridION/` — 3 occurrences confirmed in scripts.

**Assessment:** Not sensitive personal data. These are local filesystem paths pointing to the
analysis directory. They do not reveal security-relevant information. However, before
publication these paths should be replaced with relative paths or a `$PROJECT_ROOT` variable
to allow the scripts to run from any location.

**Action required:** Parameterize all absolute paths before GitHub publication. (See release checklist.)

---

### 2. Local sample identifiers (`cepa_*`, `barcode*`)

**Files affected:** Scripts in `01_genome_preparation/`, `02_typing_and_identity/`, metadata files.

**Pattern:** Sample names of the form `cepa_17_barcode01` etc. (8 local Uruguay isolates).

**Assessment:** These are internal laboratory sample codes, not personally identifying
information (no producer names, farm names, or patient data). The sample IDs are scientifically
relevant and will appear in the deposited sequence records. No action required.

---

### 3. No passwords, API keys, or credentials found

Scan for: `password`, `token`, `api_key`, `secret`, `ssh-rsa`, `BEGIN RSA PRIVATE KEY`,
NCBI API keys, AWS credentials.

**Result:** None found.

---

### 4. No producer/farm/establishment names found

Scan for: Spanish-language identifiers that could reveal the origin farms
(`productor`, `predio`, `tambo`, `establecimiento`).

**Result:** None found in this repository. The local isolates are referred to only by
internal sample IDs.

---

### 5. No personal email addresses found

Scan for: `@gmail`, `@hotmail`, email address patterns.

**Result:** None found in scripts or metadata.

---

### 6. metadata/barcode_muestras.tsv — NOT present in this repo

The file `metadata/barcode_muestras.tsv` referenced by `01_genome_preparation/` scripts
(contains barcode-to-sample mapping for local isolates) is NOT copied into this repository
because it may contain internal lab identifiers. Replace with a redacted version listing
only sample IDs, barcode numbers, and sequence types — no producer or farm identifiers —
before GitHub publication.

---

## Summary

| Finding | Sensitive? | Action |
|---------|-----------|--------|
| `/Volumes/JM/BioInf/GridION/` paths | No | Parameterize before publish |
| `cepa_*` sample IDs | No (lab codes) | None required |
| Passwords / credentials | Not found | None required |
| Farm / producer names | Not found | None required |
| Personal emails | Not found | None required |
| `barcode_muestras.tsv` | Review before adding | Redact if needed |
