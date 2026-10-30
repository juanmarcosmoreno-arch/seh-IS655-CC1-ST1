#!/usr/bin/env python3
"""Figure S2. IS655-like coding capacity schematic (v2 — fixed layout).

Labels inside arrows; stop codons only above arrow; legend below panels.
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUT = Path("/Volumes/JM/BioInf/GridION/results/seh_st1_bovine/final_paper_figures")
OUT.mkdir(exist_ok=True)

ELEM_LEN = 1211
ARROW_H  = 0.36          # arrow height (data units)

# ── Palette ──────────────────────────────────────────────────────────────────
C_ORF1_REF = "#1565C0"   # dark blue  — IS655 ref
C_ORF2_REF = "#B71C1C"   # dark red   — IS655 ref
C_ORF1_CAN = "#90CAF9"   # light blue — IS655-like
C_ORF2_CAN = "#EF9A9A"   # light red  — IS655-like
C_OVL      = "#FF8F00"   # amber      — overlap region
C_STOP     = "#C62828"   # red        — stop codons
C_DEL      = "#6A1B9A"   # purple     — deletion
C_DDE_OK   = "#1B5E20"   # dark green — intact catalytic residue
C_DDE_BAD  = "#D32F2F"   # red        — mutated residue
C_BB       = "#455A64"   # backbone / scale

# ── IS655 reference: ORF coords & DDE positions ───────────────────────────
# ORF1: 1–381 bp; ORF2: 345–1211 bp; overlap 345–381 (37 bp)
# DDE aa in ORF2: D139, D157, E192 → nt = orf2_start + (aa-1)*3
REF_ORF2_START = 345
CAN_ORF2_START = 361

ref_dde = [
    (REF_ORF2_START + (139-1)*3, "D¹³⁹", C_DDE_OK),
    (REF_ORF2_START + (157-1)*3, "D¹⁵⁷", C_DDE_OK),
    (REF_ORF2_START + (192-1)*3, "E¹⁹²", C_DDE_OK),
]
can_dde = [
    (CAN_ORF2_START + (139-1)*3, "D¹³⁹", C_DDE_OK),
    (CAN_ORF2_START + (157-1)*3, "D¹⁵⁷", C_DDE_OK),
    (CAN_ORF2_START + (192-1)*3, "N¹⁹²", C_DDE_BAD),   # E→N mutation
]

PANELS = [
    dict(
        title    = "IS655 (IS3 family) — reference element",
        subtitle = "Intact ORFs | DDE catalytic triad | A₇G frameshift signal → autonomous transposition competent",
        orf1     = (1, 381),    orf1_color = C_ORF1_REF, orf1_tcol = "white",
        orf2     = (345, 1211), orf2_color = C_ORF2_REF, orf2_tcol = "white",
        ovl      = (345, 381),
        stops    = [],
        dde      = ref_dde,
        frameshift = True,
        del_pos  = None,
    ),
    dict(
        title    = "IS655-like (ST1 bovine) — 42/46 loci identical",
        subtitle = "57.4% nt identity to IS655 ORF1  |  49.0% nt identity to IS655 ORF2  |  1,211 bp",
        orf1     = (1, 397),    orf1_color = C_ORF1_CAN, orf1_tcol = "#0D47A1",
        orf2     = (361, 1211), orf2_color = C_ORF2_CAN, orf2_tcol = "#7B1FA2",
        ovl      = (361, 397),
        stops    = [(161, "TAA"), (272, "TGA"), (1069, "TAA")],
        dde      = can_dde,
        frameshift = False,
        del_pos  = 630, del_len = 45,    # 15 aa × 3 nt
    ),
]

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(14, 9.5),
                          gridspec_kw={"height_ratios": [1, 1]})
# bottom=0.26 reserves space below panels for the legend
fig.subplots_adjust(hspace=0.70, left=0.04, right=0.97, top=0.93, bottom=0.26)


def draw_panel(ax, cfg):
    ax.set_xlim(-25, ELEM_LEN + 30)
    ax.set_ylim(-1.05, 2.10)
    ax.axis("off")
    y0 = 0.0

    # ── Title & subtitle ────────────────────────────────────────────────────
    ax.text(ELEM_LEN / 2, 2.00, cfg["title"],
            ha="center", va="top", fontsize=10.5, fontweight="bold")
    ax.text(ELEM_LEN / 2, 1.76, cfg["subtitle"],
            ha="center", va="top", fontsize=8, fontstyle="italic", color="#546E7A")

    # ── Backbone + TIR marks ────────────────────────────────────────────────
    ax.hlines(y0, 0, ELEM_LEN, color=C_BB, lw=1.5, zorder=1)
    for x in (0, ELEM_LEN):
        ax.vlines(x, y0 - 0.22, y0 + 0.22, color=C_BB, lw=2.5)
    ax.text(-22, y0, "TIR", ha="right", va="center", fontsize=6.5, color=C_BB)
    ax.text(ELEM_LEN + 22, y0, "TIR", ha="left",  va="center", fontsize=6.5, color=C_BB)

    # ── Scale bar (top-right, clear of annotations) ─────────────────────────
    sb_x0, sb_x1 = ELEM_LEN - 240, ELEM_LEN - 40
    sb_y = y0 + 1.45
    ax.hlines(sb_y, sb_x0, sb_x1, color=C_BB, lw=1.5)
    ax.vlines([sb_x0, sb_x1], sb_y - 0.04, sb_y + 0.04, color=C_BB, lw=1.5)
    ax.text((sb_x0 + sb_x1) / 2, sb_y + 0.08, "200 bp",
            ha="center", va="bottom", fontsize=7.5, color=C_BB)

    # ── Arrow helper (label INSIDE the arrow body) ──────────────────────────
    def draw_arrow(x0, x1, color, lines, tcol):
        w = x1 - x0
        head_len = 18
        patch = mpatches.FancyArrow(
            x0, y0, w - head_len, 0,
            width=ARROW_H, head_width=ARROW_H, head_length=head_len,
            fc=color, ec="white", lw=0.6,
            length_includes_head=True, zorder=3
        )
        ax.add_patch(patch)
        # Text inside: avoid the arrowhead (place in the first 55% of body)
        tx = x0 + (w - head_len) * 0.42
        n = len(lines)
        for i, line in enumerate(lines):
            ty = y0 + (n - 1) * 0.065 - i * 0.13
            ax.text(tx, ty, line, ha="center", va="center",
                    fontsize=8.0, fontweight="bold", color=tcol, zorder=5)

    # Draw ORF2 first (behind ORF1)
    draw_arrow(*cfg["orf2"], cfg["orf2_color"],
               ["ORF2", "(insB-like transposase)"], cfg["orf2_tcol"])
    draw_arrow(*cfg["orf1"], cfg["orf1_color"],
               ["ORF1", "(insA-like)"], cfg["orf1_tcol"])

    # ── Overlap region ──────────────────────────────────────────────────────
    ov0, ov1 = cfg["ovl"]
    ax.add_patch(mpatches.Rectangle(
        (ov0, y0 - ARROW_H / 2), ov1 - ov0, ARROW_H,
        fc=C_OVL, ec="white", alpha=0.90, zorder=4
    ))
    ax.text((ov0 + ov1) / 2, y0 - ARROW_H / 2 - 0.13,
            f"overlap: {ov1-ov0+1} bp",
            ha="center", va="top", fontsize=6.5, color=C_OVL)

    # ── A7G / No A7G ────────────────────────────────────────────────────────
    motif_x = (ov0 + ov1) / 2
    if cfg["frameshift"]:
        ax.plot(motif_x, y0 - ARROW_H / 2 - 0.42, "^",
                color="#2E7D32", markersize=9, zorder=5)
        ax.text(motif_x, y0 - ARROW_H / 2 - 0.56,
                "A₇G\nframeshift\nmotif",
                ha="center", va="top", fontsize=6.5, color="#2E7D32")
    else:
        ax.plot(motif_x, y0 - ARROW_H / 2 - 0.42, "x",
                color="#B71C1C", markersize=10, markeredgewidth=2.5, zorder=5)
        ax.text(motif_x, y0 - ARROW_H / 2 - 0.56,
                "No A₇G\nframeshift\nmotif",
                ha="center", va="top", fontsize=6.5, color="#B71C1C")

    # ── DDE/DDN catalytic residues (ticks through arrow + label above) ──────
    for (nt_pos, lbl, col) in cfg["dde"]:
        ax.vlines(nt_pos, y0 - ARROW_H / 2, y0 + ARROW_H / 2,
                  color=col, lw=2.2, zorder=6)
        # Stagger DDE labels above the arrow to avoid overlap with stop codons:
        # D139 → low (just above arrow), D157 → mid, E192/N192 → mid
        ax.text(nt_pos, y0 + ARROW_H / 2 + 0.06, lbl,
                ha="center", va="bottom", fontsize=7.0,
                fontweight="bold", color=col, zorder=9)

    # ── Stop codons (above the arrow; dashed tick from arrow top to triangle) ─
    # Positions stagger slightly to avoid label collision for nearby stops
    offsets = [0, 0.20, 0]    # extra y offset per stop (index order)
    for idx, (nt_pos, codon) in enumerate(cfg["stops"]):
        extra = offsets[idx] if idx < len(offsets) else 0
        tri_y  = y0 + ARROW_H / 2 + 0.30 + extra
        text_y = tri_y + 0.14 + extra
        # Short tick from top of arrow to triangle
        ax.vlines(nt_pos, y0 + ARROW_H / 2, tri_y - 0.05,
                  color=C_STOP, lw=1.8, ls="--", zorder=7)
        ax.plot(nt_pos, tri_y, "v",
                color=C_STOP, markersize=9, zorder=8)
        ax.text(nt_pos, text_y, f"stop  {codon}",
                ha="center", va="bottom", fontsize=7.0,
                color=C_STOP, fontweight="bold", zorder=9)

    # ── 15-aa deletion in ORF2 ──────────────────────────────────────────────
    if cfg.get("del_pos"):
        dp = cfg["del_pos"]
        dl = cfg["del_len"]
        d0, d1 = dp - dl / 2, dp + dl / 2
        del_y = y0 - ARROW_H / 2 - 0.55
        ax.annotate("", xy=(d1, del_y), xytext=(d0, del_y),
                    arrowprops=dict(arrowstyle="<->", color=C_DEL, lw=1.8))
        ax.text(dp, del_y - 0.12, "∆15 aa",
                ha="center", va="top",
                fontsize=7.5, color=C_DEL, fontweight="bold")


for ax, cfg in zip(axes, PANELS):
    draw_panel(ax, cfg)

# ── Legend — placed in the white space below the panels (bottom=0.26) ────────
legend_elements = [
    mpatches.Patch(fc=C_ORF1_REF, ec="white",  label="ORF1 (insA-like) — IS655 reference"),
    mpatches.Patch(fc=C_ORF2_REF, ec="white",  label="ORF2 (insB-like transposase) — IS655 reference"),
    mpatches.Patch(fc=C_ORF1_CAN, ec="#90CAF9", label="ORF1 region — IS655-like (ST1 bovine)"),
    mpatches.Patch(fc=C_ORF2_CAN, ec="#EF9A9A", label="ORF2 region — IS655-like (ST1 bovine)"),
    mpatches.Patch(fc=C_OVL,      ec="white",  label="Frameshifting overlap region (37 bp)"),
    plt.Line2D([0], [0], color=C_DDE_OK,  lw=2.5,             label="DDE catalytic residue (intact)"),
    plt.Line2D([0], [0], color=C_DDE_BAD, lw=2.5,             label="DDN — E¹⁹²→N substitution (non-functional)"),
    plt.Line2D([0], [0], color=C_STOP, lw=1.5, ls="--",
               marker="v", markersize=7,                       label="Internal stop codon"),
    mpatches.Patch(fc=C_DEL, ec="white",                       label="15-aa in-frame deletion in ORF2 (∆45 nt)"),
]

fig.legend(
    handles=legend_elements,
    loc="lower center",
    ncol=3,
    fontsize=8.5,
    frameon=True,
    framealpha=0.97,
    edgecolor="#BDBDBD",
    bbox_to_anchor=(0.5, 0.01),   # figure coordinates — in the bottom 0.26 margin
    borderaxespad=0.8,
)

fig.suptitle(
    "Figure S2. IS655-like element coding capacity compared with IS655 (IS3 family) reference",
    fontsize=11, fontweight="bold", y=0.99
)

for fmt in ("pdf", "png"):
    fig.savefig(OUT / f"FigureS2_IS655_coding_capacity.{fmt}",
                dpi=300, bbox_inches="tight", facecolor="white")

print("Saved:", OUT)
