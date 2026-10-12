"""Stage 2 of the framework: border-based minimum covering set + sensitivity analysis.

Reads the cross-impact matrix from _cassava_appendix_matrix.tex (the appendix table),
so every number reported in the manuscript can be reproduced with:  python stage2_selection.py

Overview
--------
Input: the cross-impact matrix X, with one row per candidate action and one column per
essential factor (EF1..EF23, the factors validated as essential by both G1 and G2).
Each cell x_ij is the impact of action i on factor j, scored by Lab's staff on the
scale -2 (very negative) .. 0 (neutral) .. +2 (very positive).
Candidate actions: LA1..LA28 (literature-derived) and CA1..CA16 (proposed by Lab's team).

Selection rule (Section 4.5.2 of the manuscript):
  1. Border: for each factor j, b_j = max_i x_ij, the best score that any candidate
     action reaches on that factor.
  2. An action "strongly treats" factor j if it reaches the border (x_ij == b_j).
  3. The agenda is a set of actions of minimum size that strongly treats every
     essential factor (a minimum set-cover problem). Ties between minimum sets are
     broken outside this script, by preferring literature-derived actions (LA).

Sensitivity analysis: every matrix cell is moved by +1 and by -1 (one cell at a time,
staying within -2..2), and the minimum sets are recomputed for each perturbed matrix.

Implementation note: the set of factors that an action treats strongly is stored as a
bit mask (an integer whose bit j is 1 when factor j is treated strongly). The union of
several actions' coverage is then a bitwise OR, and "covers every factor" means the
union equals FULL (all 23 bits set). This makes the exhaustive search fast.
"""
import itertools
import re
from collections import Counter

# ---------------------------------------------------------------------------
# 1. Read the cross-impact matrix from the LaTeX appendix table.
#    Each data row looks like:  LA1 & 2 & 1 & ~0 & ... \\
#    The regex captures the action code (LA.. or CA..) and the cell values;
#    "~" (a LaTeX spacing character used for alignment) is removed before parsing.
# ---------------------------------------------------------------------------
pat = re.compile(r"\s*((?:LA|CA)\d+)\s*&(.*?)~?\s*\\\\")
rows = {}  # action code -> list of 23 integer scores (EF1..EF23, in column order)
for line in open("_cassava_appendix_matrix.tex", encoding="utf8"):
    m = pat.match(line)
    if m:
        rows[m.group(1)] = [int(v) for v in m.group(2).replace("~", "").split("&")]
keys = list(rows)       # candidate action codes, in the order of the table
n_f = 23                # number of essential factors (matrix columns)
FULL = (1 << n_f) - 1   # bit mask with all 23 factors covered


# ---------------------------------------------------------------------------
# 2. Core functions
# ---------------------------------------------------------------------------
def strong_masks(mat):
    """Compute the border and, for each action, the factors it treats strongly.

    mat: dict action -> list of scores (one per factor).
    Returns (border, masks):
      border[j] = best score reached by any action on factor j (b_j = max_i x_ij);
      masks[k]  = bit mask of the factors on which action k reaches the border.
    The border is recomputed from `mat`, so it also adapts to perturbed matrices.
    """
    border = [max(r[i] for r in mat.values()) for i in range(n_f)]
    return border, {k: sum(1 << i for i, x in enumerate(v) if x == border[i]) for k, v in mat.items()}


def minimum_covers(masks, max_size=6):
    """Find all sets of actions of minimum size that strongly treat every factor.

    Exhaustive search: try every combination of 1 action, then of 2, and so on,
    and stop at the first size r for which at least one combination covers all
    factors. Because the search goes up in size, the sets found are guaranteed to
    be minimum (no smaller set exists), and *all* minimum sets are returned, so
    ties are visible. Returns (size, list of action tuples), or (None, []) if no
    cover with at most `max_size` actions exists.
    """
    for r in range(1, max_size + 1):
        sols = [c for c in itertools.combinations(keys, r)
                if eval_cover(masks, c) == FULL]
        if sols:
            return r, sols
    return None, []


def eval_cover(masks, combo):
    """Return the bit mask of factors strongly treated by at least one action in `combo`."""
    acc = 0
    for k in combo:
        acc |= masks[k]
    return acc


# ---------------------------------------------------------------------------
# 3. Apply the rule to the observed matrix (results reported in Section 4.5.2)
# ---------------------------------------------------------------------------
border, masks = strong_masks(rows)
print("Candidate actions:", len(keys), "| border values:", sorted(set(border)))
size, sols = minimum_covers(masks)
print("Minimum cover size:", size, "| minimum covers:", sols)

# The agenda chosen after the tie-break (literature-derived actions preferred).
SEL = ("LA1", "LA3", "LA13", "LA27")

# Coverage of the selected set: factors reached at the border (strong) and
# factors with at least one positive score (x_ij >= 1).
strong = [i for i in range(n_f) if any(rows[a][i] == border[i] for a in SEL)]
positive = [i for i in range(n_f) if any(rows[a][i] >= 1 for a in SEL)]
print(f"Selected set strong coverage: {len(strong)}/23 | positive coverage: {len(positive)}/23")

# Profile of each selected action: factors with the top score (2; equal to the
# border here, since b_j = 2 for every factor) and factors with a negative score
# (side effects that another selected action must offset).
for a in SEL:
    print(a, "strong on", ["EF%d" % (i + 1) for i, x in enumerate(rows[a]) if x == 2],
          "| negative on", ["EF%d" % (i + 1) for i, x in enumerate(rows[a]) if x < 0])

# Unique strong treatments: factors that only one action of the selected set treats
# strongly. These pairs show what would be lost if that action were dropped.
unique = [(a, "EF%d" % (i + 1)) for a in SEL for i in range(n_f)
          if rows[a][i] == 2 and not any(rows[b][i] == 2 for b in SEL if b != a)]
print("Factor-action pairs with a unique strong treatment:", len(unique), unique)

# Breadth of the two broadest actions: number of factors with a positive impact.
print("Positive-impact count of LA3:", sum(x >= 1 for x in rows["LA3"]), "of 23")
print("Positive-impact count of LA1:", sum(x >= 1 for x in rows["LA1"]), "of 23")

# ---------------------------------------------------------------------------
# 4. Sensitivity analysis: every single-cell +/-1 perturbation that stays on the -2..2 scale
#    For each action k, factor j and direction d, one cell is changed, the border and
#    strong masks are recomputed, and all minimum covers are found again. Counted:
#      - total admissible perturbations;
#      - cases in which the selected set is still one of the minimum covers;
#      - cases in which the core LA1+LA3+LA13 is part of some minimum cover;
#      - the distribution of the minimum cover size.
# ---------------------------------------------------------------------------
stats = Counter()
min_sizes = Counter()
for k in keys:
    for j in range(n_f):
        for d in (-1, 1):
            x = rows[k][j] + d
            if not -2 <= x <= 2:
                continue  # perturbation would leave the scale: not admissible
            mat = dict(rows)       # shallow copy of the matrix...
            mat[k] = rows[k][:]    # ...with a fresh copy of the row being changed
            mat[k][j] = x
            b, m = strong_masks(mat)
            sz, ss = minimum_covers(m)
            min_sizes[sz] += 1
            stats["total"] += 1
            if SEL in ss:
                stats["selected set still a minimum cover"] += 1
            if any({"LA1", "LA3", "LA13"} <= set(s) for s in ss):
                stats["LA1+LA3+LA13 in some minimum cover"] += 1
print("Sensitivity:", dict(stats), "| minimum cover sizes:", dict(min_sizes))
