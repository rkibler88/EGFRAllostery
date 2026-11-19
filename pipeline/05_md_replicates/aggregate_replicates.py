#!/usr/bin/env python3
"""Between-replicate agreement for CORRECTIONS.md Section 9.1.

Reads egfr_dccm_matrix_rep*.csv (and optionally the published
egfr_dccm_matrix.csv as replicate 1), reports pairwise Pearson r over the
strict upper triangle, and writes the replicate-mean DCCM.

Accepts any of the CSV layouts in this project: the published pandas
to_csv form (header + residue-ID index column), header-only, or bare numeric.
"""
import argparse, glob, itertools, os, re
import numpy as np


def load_dccm(path):
    """Read a DCCM CSV in any of the layouts this project produces.

    Handles: pandas to_csv (header row + residue-ID index column, as in the
    published egfr_dccm_matrix.csv), header-only, index-only, and bare numeric.
    Returns (matrix, labels_or_None).
    """
    with open(path) as fh:
        rows = [ln.rstrip("\n").split(",") for ln in fh if ln.strip()]
    if not rows:
        raise SystemExit(f"{path}: empty")

    def isnum(tok):
        try:
            float(tok)
            return True
        except ValueError:
            return False

    widths = {len(r) for r in rows}
    if len(widths) != 1:
        raise SystemExit(f"{path}: ragged rows, widths {sorted(widths)}")

    # A header is present if its leading field is empty/non-numeric, or if there
    # is exactly one more line than there are columns (a numeric residue-ID
    # header, which `isnum` alone cannot distinguish from data).
    header = None
    if (rows[0][0].strip() == "" or not isnum(rows[0][0])
            or len(rows) == len(rows[0]) + 1):
        header = rows.pop(0)
    n_rows = len(rows)
    n_cols = len(rows[0])

    labels = None
    if n_cols == n_rows + 1:            # leading index column
        labels = [r[0].strip() for r in rows]
        rows = [r[1:] for r in rows]
    elif n_cols != n_rows:
        raise SystemExit(f"{path}: {n_rows} rows x {n_cols} cols is neither "
                         f"square nor square-plus-index")

    if labels is None and header is not None:
        h = header[1:] if len(header) == n_rows + 1 else header
        if len(h) == n_rows:
            labels = [x.strip() for x in h]

    try:
        m = np.array([[float(x) for x in r] for r in rows], dtype=float)
    except ValueError as e:
        raise SystemExit(f"{path}: non-numeric cell ({e})")
    if m.shape[0] != m.shape[1]:
        raise SystemExit(f"{path}: not square after parsing {m.shape}")
    return m, labels


def write_dccm(path, m, labels):
    """Write in the published layout: header row with empty leading field,
    residue ID prefixing each row."""
    lab = labels if labels is not None else [str(i) for i in range(m.shape[0])]
    with open(path, "w") as fh:
        fh.write("," + ",".join(lab) + "\n")
        for name, row in zip(lab, m):
            fh.write(name + "," + ",".join(f"{v:.10f}" for v in row) + "\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--glob", default="egfr_dccm_matrix_rep*.csv")
    p.add_argument("--published", default=None,
                   help="egfr_dccm_matrix.csv from the original run, included as rep1")
    p.add_argument("--out-prefix", default="replicate")
    a = p.parse_args()

    mats, labs, names = [], [], []
    if a.published:
        if not os.path.exists(a.published):
            raise SystemExit(f"--published {a.published} not found")
        m, l = load_dccm(a.published)
        mats.append(m); labs.append(l); names.append("rep1_published")
    for f in sorted(glob.glob(a.glob)):
        m, l = load_dccm(f)
        mats.append(m); labs.append(l)
        names.append(re.sub(r"\.csv$", "", os.path.basename(f)).replace("egfr_dccm_matrix_", ""))
    if len(mats) < 2:
        raise SystemExit(f"need >=2 matrices, found {len(mats)} "
                         f"(glob {a.glob!r} in {os.getcwd()})")

    shapes = {m.shape for m in mats}
    if len(shapes) != 1:
        raise SystemExit("shape mismatch: " +
                         ", ".join(f"{n}{m.shape}" for n, m in zip(names, mats)))

    # residue labels must agree wherever they are present
    known = [(n, l) for n, l in zip(names, labs) if l is not None]
    for n, l in known[1:]:
        if l != known[0][1]:
            raise SystemExit(f"residue labels differ between {known[0][0]} and {n}")
    ref_labels = known[0][1] if known else None

    n = mats[0].shape[0]
    iu = np.triu_indices(n, k=1)
    rows = [(names[ia], names[ib],
             float(np.corrcoef(mats[ia][iu], mats[ib][iu])[0, 1]))
            for ia, ib in itertools.combinations(range(len(mats)), 2)]

    with open(f"{a.out_prefix}_dccm_correlation.csv", "w") as fh:
        fh.write("replicate_a,replicate_b,pearson_r\n")
        for x, y, r in rows:
            fh.write(f"{x},{y},{r:.4f}\n")

    write_dccm(f"{a.out_prefix}_dccm_mean.csv", np.mean(mats, axis=0), ref_labels)

    rs = np.array([r for _, _, r in rows])
    print(f"n_replicates={len(mats)} ({', '.join(names)})  matrix={n}x{n}")
    print(f"pairwise upper-triangle Pearson r: min={rs.min():.4f} "
          f"median={np.median(rs):.4f} max={rs.max():.4f}")
    with_pub = [r for x, y, r in rows if "published" in (x + y)]
    if with_pub:
        wp = np.array(with_pub)
        print(f"  vs published replicate 1:        min={wp.min():.4f} "
              f"median={np.median(wp):.4f} max={wp.max():.4f}")
    print(f"wrote {a.out_prefix}_dccm_correlation.csv and {a.out_prefix}_dccm_mean.csv")


if __name__ == "__main__":
    main()
