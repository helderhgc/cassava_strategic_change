# Stage 2 action selection for a cassava-processing MSE

This repository contains the code and data needed to reproduce the Stage 2 results
(action selection and sensitivity analysis) of the manuscript
*Navigating Strategic Change in Agribusiness Micro and Small Enterprises: A
Decision-Support Framework for Identifying Essential Factors and Prioritizing Actions*
(under review at *Strategic Change*).

## Contents

| File | Description |
|---|---|
| `_cassava_appendix_matrix.tex` | Cross-impact matrix (appendix table of the manuscript): impact of each of the 44 candidate actions (LA1–LA28 from the literature, CA1–CA16 proposed by the firm's team) on each of the 23 essential factors (EF1–EF23), on a scale from −2 to +2. |
| `stage2_selection.py` | Reads the matrix, applies the selection rule, and runs the sensitivity analysis. |

## Selection rule

For each essential factor, the *border* is the best score that any candidate action
reaches on it. An action *strongly treats* a factor when it reaches the border. The
agenda is a set of actions of minimum size that strongly treats every essential factor
(a minimum set-cover problem, solved exactly by enumeration). Ties between minimum sets
are broken by preferring literature-derived actions (LA).

The sensitivity analysis changes each cell of the matrix by +1 and by −1, one cell at a
time and within the −2..2 scale (1,830 perturbations), and recomputes the minimum sets.

## How to run

Requires Python 3.8 or later; no external packages are needed.

```
python stage2_selection.py
```

Run the command from this folder, since the script reads `_cassava_appendix_matrix.tex`
from the current directory. The full run (including the sensitivity analysis) takes a
few minutes.

Expected main results:

- minimum set size: 4;
- minimum sets: {LA1, LA3, LA13, LA27}, {LA1, LA3, LA13, CA3}, {LA1, LA3, LA13, CA10};
- sensitivity: the selected set remains a minimum set in 1,814 of 1,830 perturbations,
  and LA1, LA3 and LA13 belong to a minimum set in all 1,830.

## Citation

If you use this code or data, please cite the article (full reference to be added
after publication).

## License

This work is licensed under the Creative Commons Attribution 4.0 International License
(CC BY 4.0). See `LICENSE`.
