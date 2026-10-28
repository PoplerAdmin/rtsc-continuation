# Literature and source check: irreducible-wedge q sums for α²F

**Status: NOT ADOPTED.** User decision D3 (2026-10-27; see EPH_CONVERGENCE_PROTOCOL.md §D): fine q grids stay full and uniform. This note is kept for the record.

Date: 2026-10-27. Every item was checked against a paper, a source file or a binary during this session. Short quotes only. "UNCONFIRMED" marks what could not be checked.

## Precedent

| source | what it does | how it was checked |
|---|---|---|
| Savrasov & Savrasov, PRB 54, 16487 (1996), arXiv cond-mat/9604124 | α²F for Al, Cu, Mo, Nb, Pb, Pd, Ta, V from "a set of irreducible q-points at the (8,8,8)-reciprocal lattice grid"; "29 points per 1/48th part of the Brillouin zone" | ar5iv full text |
| Wierzbowska, de Gironcoli & Giannozzi, arXiv cond-mat/0504077 (Nb) | "Symmetry is exploited to reduce the number of k- and q-points"; k sum "restricted to the irreducible BZ calculated with respect to the small group of q"; their Eq. 6 is α²F = (1/2πN(ε_F)) Σ_qν γ_qν/ω_qν δ(ω−ω_qν), with the 1/N_q normalisation given in App. Eq. 23 | ar5iv / arXiv HTML |
| QE ph.x `elphsum` (PHonon/PH/elphon.f90, develop) | e-ph is computed only at irreducible q; the k sum covers the full BZ through the small-group rotation. The star is **unfolded by exact rotation** of the γ matrix (`star_q`, `q2qstar_ph`, elphon.f90:1181, including −q), then Fourier-interpolated (q2r, matdyn `la2F`) to a full grid for a2F.dos | raw source, spot-checked by me |
| QE `lambda.x` (PHonon/PH/lambda.f90) | reads `q(1..3) wk` per irreducible q, normalises the weights, and forms λ = Σ wk·λ_q. The weights are user-supplied and unchecked | raw source, spot-checked by me |
| Lee et al., EPW 2023, arXiv 2302.08085 | "it is useful to also exploit symmetry operations on the fine grids, so as to reduce computational cost", with the caveat that "quantities on the fine grids are not expected to fulfill the required symmetry relations" because Wannierisation does not impose symmetry | ar5iv, read by an agent |
| Poncé et al., EPW 2016, CPC 209, 116 | symmetry is used on the coarse grid only (IBZ→BZ unfolding, time reversal); fine grids are homogeneous, random or Sobol | pdftotext, read by an agent |
| Margine & Giustino, PRB 87, 024505 (2013) | Eliashberg gap computed on the irreducible fine k (`mp_mesh_k`) | EPW slides, docs |

## EPW support (corrects the earlier note)

- **`mp_mesh_q` is live code.** The documentation calls it "Not currently in use", which is stale. In EPW/src/utilities/bzgrid.f90 (`loadqmesh_para` ~L1750), `IF (mp_mesh_q)` calls `kpoint_grid(nsym, time_reversal, ...)`, which gives a **Γ-centred** irreducible-wedge fine q mesh with star weights `wqf`.
- **Our installed binary (conda QE 7.6 / EPW 6.1) contains it.** `strings epw.x` shows `mp_mesh_q` and "Using uniform MP q-mesh".
- **readin.f90 rules:**
  - phonselfen forbids `mp_mesh_k` ("phonon self-energy only works with full uniform k-mesh");
  - elecselfen forbids `mp_mesh_q`;
  - the two cannot both be true.
  - So **phonselfen + mp_mesh_q is the combination EPW permits**: irreducible q, full-BZ k for each q. This matches ph.x and Savrasov.
- selfen.f90 prints `wt:` per q, so the star weights are visible in stdout and can be cross-checked against our own `kpoint_grid` reproduction.

## Theory references

- Maradudin & Vosko, RMP 40, 1 (1968): bibliographic details confirmed; it is the group-theory treatment of normal modes, including time reversal. UNCONFIRMED: no verbatim ω(Sq)=ω(q) quote (full text not reachable).
- Monkhorst & Pack, PRB 13, 5188 (1976): details confirmed. VASP wiki: "Use only Γ-centered meshes for face-centered cubic (fcc), hexagonal, and fcc-orthorhombic" lattices, because even MP meshes are shifted by 1/2.
- Invariance of γ_qν (summed over degenerate subspaces) on the star: **no verbatim primary-source sentence found.** It is standard physics, and QE encodes it, since it rotates γ to the star like the dynamical matrix. A 2006 pw_forum thread reports numerical γ differences within a star (unresolved; page did not load).
- Engel, Marsman, Franchini & Kresse, PRB 101, 184302 (2020), arXiv 2001.08666 (PAW + Wannier, VASP): its matrix elements "lack some symmetries that are usually present" for simple eigenvalue problems. Not directly EPW/NC, but the same caution applies. (Correction: the NSF PAR hit 10517662 is Li et al., CPC 295, 109003, GWPT, which is about gauge-recovering symmetry unfolding.)
- Sakuma, PRB 87, 235109 (2013): MLWFs "do not necessarily reflect the spatial symmetry"; symmetry-adapted WFs.

## Consequences for the amendment

1. Use EPW's native `mp_mesh_q = .true.` with Γ-centred nqf (no shifted grids). Read weights from stdout and verify them against an independent kpoint_grid (spglib) reproduction; any mismatch → `EPH_PIN_VIOLATION`.
2. **Star test (mandatory):** interpolate full stars, including −q, of a random sample of irreducible q via filqf. Compare degenerate-subspace sums of γ and ω. The measured maximum relative deviation enters e; above threshold → `UNCONVERGED_EPH(sym)`.
3. **Full vs IBZ on Al at 16³:** A, λ and ω_log within admission tolerances.
4. G6 (direct DFPT γ at coarse q) and the tetrahedron anchor stay as they are; they are full-grid independent checks.
5. Independent adversarial review per HANDOFF §3.
