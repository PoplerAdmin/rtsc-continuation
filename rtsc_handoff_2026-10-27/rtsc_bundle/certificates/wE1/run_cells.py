"""Adaptive cover of [W_LO, W_HI] by certified cells (even sector).
usage: python3 run_cells.py A0 ETA [N K WLO WHI LOG]
A0, ETA are decimal strings; the certificate is for the exact decimals (arb balls containing them)."""
import sys, time, json
from flint import arb
from cell import certify_cell


def initial_cells(wlo, whi):
    # dyadic breakpoints: width 1/2 on [1,16], 2 on [16,64], 8 on [64,whi]
    pts = [wlo]
    while pts[-1] < whi:
        w = pts[-1]
        step = 0.5 if w < 16 else (2.0 if w < 64 else 8.0)
        pts.append(min(w + step, whi))
    return list(zip(pts[:-1], pts[1:]))


def main():
    A0 = sys.argv[1]; ETA = sys.argv[2]
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 150
    K = int(sys.argv[4]) if len(sys.argv) > 4 else 600
    WLO = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    WHI = float(sys.argv[6]) if len(sys.argv) > 6 else 400.0
    LOG = sys.argv[7] if len(sys.argv) > 7 else "cells_%s.log" % A0
    a = arb(A0); eta = arb(ETA)
    MINW = 2.0 ** -16
    todo = initial_cells(WLO, WHI)
    done, failed = [], []
    t0 = time.time()
    with open(LOG, "w") as f:
        f.write("# even-sector cells a0=%s eta=%s N=%d K=%d range=[%r,%r]\n" % (A0, ETA, N, K, WLO, WHI))
        while todo:
            w1, w2 = todo.pop(0)
            ok, info = certify_cell(a, eta, w1, w2, N=N, K=K)
            if ok:
                done.append((w1, w2))
                f.write("OK   [%.10g, %.10g]  scaled-margin>=%.2e  rem0=%.2e d0=%.2e  t=%.2fs\n"
                        % (w1, w2, info["margin"], info["rem0"], info["d0"], info["time"]))
            else:
                f.write("split [%.10g, %.10g]  %s\n" % (w1, w2, {k: info.get(k) for k in info if k.startswith("ev") or k == "why"}))
                if w2 - w1 <= MINW:
                    failed.append((w1, w2))
                else:
                    m = (w1 + w2) / 2
                    todo[0:0] = [(w1, m), (m, w2)]
            f.flush()
        # check that the certified cells tile [WLO, WHI] exactly
        done.sort()
        tiled = done[0][0] == WLO and done[-1][1] == WHI and all(
            done[i][1] == done[i + 1][0] for i in range(len(done) - 1))
        verdict = (not failed) and tiled
        msg = "VERDICT cells: %s  (%d cells, %d failed, tiled=%s, %.1fs)" % (
            "CERTIFIED" if verdict else "NOT CERTIFIED", len(done), len(failed), tiled, time.time() - t0)
        f.write(msg + "\n")
    print(msg)
    json.dump({"a0": A0, "eta": ETA, "N": N, "K": K, "wlo": WLO, "whi": WHI, "cells": done,
               "failed": failed, "verdict": verdict}, open(LOG.replace(".log", ".json"), "w"))
    return verdict


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
