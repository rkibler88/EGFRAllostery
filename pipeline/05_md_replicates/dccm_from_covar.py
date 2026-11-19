#!/usr/bin/env python3
"""covar.ascii (3N x 3N Cartesian covariance from `gmx covar -ascii`)
-> N x N dynamic cross-correlation matrix, normalised to unit diagonal,
subset to the residues scored in the co-evolution analysis.

DCCM_ij = Tr(Cov_ij) / sqrt(Tr(Cov_ii) * Tr(Cov_jj))
where Cov_ij is the 3x3 block for the Ca atoms of residues i and j.
"""
import argparse
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--covar", required=True)
    p.add_argument("--n-ca", type=int, required=True, help="Ca atoms in covar group")
    p.add_argument("--traj-first", type=int, required=True, help="first residue in construct")
    p.add_argument("--score-first", type=int, required=True)
    p.add_argument("--score-last", type=int, required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    n3 = 3 * a.n_ca
    vals = np.loadtxt(a.covar).ravel()
    if vals.size != n3 * n3:
        raise SystemExit(f"covar.ascii has {vals.size} values, expected {n3*n3} "
                         f"for {a.n_ca} Ca atoms -- check the group selected in gmx covar")
    cov = vals.reshape(n3, n3)

    # trace of each 3x3 residue-residue block
    blocks = cov.reshape(a.n_ca, 3, a.n_ca, 3)
    c = np.einsum("iajb->ij", blocks * np.eye(3)[None, :, None, :])

    d = np.sqrt(np.diag(c))
    if np.any(d <= 0):
        raise SystemExit("non-positive diagonal variance -- trajectory or group is wrong")
    dccm = c / np.outer(d, d)
    np.fill_diagonal(dccm, 1.0)

    i0 = a.score_first - a.traj_first
    i1 = a.score_last - a.traj_first + 1
    if i0 < 0 or i1 > a.n_ca:
        raise SystemExit("scored range falls outside the construct")
    sub = dccm[i0:i1, i0:i1]
    resids = np.arange(a.score_first, a.score_last + 1)
    assert sub.shape == (resids.size, resids.size)

    # published layout: header row with empty leading field, residue ID per row
    with open(a.out, "w") as fh:
        fh.write("," + ",".join(str(r) for r in resids) + "\n")
        for rid, row in zip(resids, sub):
            fh.write(f"{rid}," + ",".join(f"{v:.10f}" for v in row) + "\n")
    print(f"wrote {a.out}  shape={sub.shape}  "
          f"offdiag[min,max]=[{sub[~np.eye(len(sub),dtype=bool)].min():.3f},"
          f"{sub[~np.eye(len(sub),dtype=bool)].max():.3f}]")


if __name__ == "__main__":
    main()
