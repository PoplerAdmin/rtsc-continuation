/-
DiscoveryLogic.Soundness: `counts_for` (`rtsc/discovery/goals.py`), the single predicate that
decides whether a ledger row makes a material satisfy a goal, and what it means in the model
(via `PipelineLogic.classifyFull_certified_sound`).
-/
import DiscoveryLogic.Ledger
import PipelineLogic.Multiband

namespace DiscoveryLogic

open TcR1 PipelineLogic

noncomputable section

/-- **`counts_for(goal, row, certify_hash=ch, dft_hash=dh)`** (`goals.py:211-233`):
status CERTIFIED_CANDIDATE, stage `pipeline`, origin `rtsc.certify.classify`,
`T_star_K == goal.T_star_K`, material pressure in the CURRENT goal's range, DFT pressure known,
in range and equal to the material pressure up to `tol` (`1e-9` in Python), `certify_hash`
and `dft_hash` equal to the current ones. -/
def countsFor (g : Goal) (ch dh : ℕ) (tol : ℝ) (e : Event) : Prop :=
  e.status = .certifiedCandidate ∧ e.stage = .pipeline ∧ e.origin = .classify ∧
    e.Tstar = g.Tstar ∧ e.pressure ≤ g.Pmax ∧
    (∃ pd, e.dftPressure = some pd ∧ pd ≤ g.Pmax ∧ |pd - e.pressure| ≤ tol) ∧
    e.certifyHash = some ch ∧ e.dftHash = some dh

/-- `certified_ids` (`loop.py:169-175`; reports: `report.py` `goal_rows`): `counts_for` applied to the CURRENT row of
`(m, g)`. -/
def Satisfies (L : List Event) (g : Goal) (ch dh : ℕ) (tol : ℝ) (m : ℕ) : Prop :=
  ∃ e, currentEvent L m g.id = some e ∧ countsFor g ch dh tol e

/-- Run parameters of `rtsc.certify` that a `certify_hash` stands for (`[certify]` table):
the rejection constant `a0`, the envelope `e`, `μ*_hi` and the cutoff rule. -/
structure CertParams where
  a0 : ℝ
  e : ℝ
  μhi : ℝ
  Ncf : ℝ → WithTop ℕ

/-- **Hypothesis on the pipeline (not proved here; it is the PipelineLogic contract plus the
record keeping).** Every `pipeline` / `rtsc.certify.classify` row with hashes `ch`, `dh` and
DFT pressure `pd` was produced by `classifyFull` at the row's `T_star_K`, under the `[certify]`
parameters `params ch`, from alternative spectra that the DFT stage computed for that material
under configuration `dh` at pressure `pd` (`Computed`), each with sound evidence
(`EvidenceSound`: every `true` boolean is backed by its rigorous certificate). -/
def PipelineSound {ι : Type*} (κ : ℝ) (params : ℕ → CertParams)
    (Computed : ℕ → ℕ → ℝ → List (Spectrum ι) → Prop) (L : List Event) : Prop :=
  ∀ ev ∈ L, ev.stage = .pipeline → ev.origin = .classify →
    ∀ ch dh pd, ev.certifyHash = some ch → ev.dftHash = some dh → ev.dftPressure = some pd →
      ∃ (alts : List (Spectrum ι × AltEvidence)) (mb : List Bool) (model : RejectionModel),
        Computed ev.mat dh pd (alts.map Prod.fst) ∧
        (∀ a ∈ alts, EvidenceSound κ (params ch).a0 (params ch).e (params ch).μhi ev.Tstar
          (params ch).Ncf a.1 a.2) ∧
        classifyFull (alts.map Prod.snd) mb model = ev.status

/-- **Goal soundness (post-fix design: `counts_for` is the definition).** If material `m`
satisfies goal `g` under the current configuration hashes `ch`, `dh`, then (given the pipeline
hypothesis) its material pressure and the DFT pressure actually used are in `g`'s range and
agree up to `tol`, and there is a nonempty list of alternative spectra computed for `m` under
the CURRENT DFT configuration at that pressure such that for each of them there is
`T_low ≥ g.T*` at which every spectrum of its envelope, for every `μ* ≤ μ_hi` of the CURRENT
`[certify]` configuration, is unstable in the even sector: `Tc_model ≥ g.T*` in the isotropic
single-band model (`classifyFull_certified_sound`; nothing is claimed about the multiband
model). -/
theorem satisfies_sound {ι : Type*} (κ : ℝ) (params : ℕ → CertParams)
    (Computed : ℕ → ℕ → ℝ → List (Spectrum ι) → Prop) {L : List Event}
    (hP : PipelineSound κ params Computed L) {g : Goal} {ch dh : ℕ} {tol : ℝ} {m : ℕ}
    (hs : Satisfies L g ch dh tol m) :
    ∃ e ∈ L, e.mat = m ∧ e.goal = g.id ∧ e.pressure ≤ g.Pmax ∧
      ∃ pd, pd ≤ g.Pmax ∧ |pd - e.pressure| ≤ tol ∧
      ∃ alts : List (Spectrum ι × AltEvidence), Computed m dh pd (alts.map Prod.fst) ∧
        alts ≠ [] ∧ ∀ a ∈ alts, ∃ T : ℝ, g.Tstar ≤ T ∧
          ∀ sp' : Spectrum ι, InEnv (params ch).e a.1 sp' → ∀ μ : ℝ, μ ≤ (params ch).μhi →
            Unstable (Qp κ sp' μ ((params ch).Ncf T) T) ∧
              ¬ Stable κ sp' μ ((params ch).Ncf T) T := by
  obtain ⟨e, he, hst, hstage, horig, hT, hp, ⟨pd, hpd, hpdr, hpdt⟩, hch, hdh⟩ := hs
  obtain ⟨hmem, hm, hg⟩ := currentEvent_spec he
  obtain ⟨alts, mb, model, hcomp, hev, hcl⟩ := hP e hmem hstage horig ch dh pd hch hdh hpd
  rw [hst] at hcl
  have := classifyFull_certified_sound κ (params ch).a0 (params ch).e (params ch).μhi e.Tstar
    (params ch).Ncf alts hev mb model hcl
  rw [hT] at this
  refine ⟨e, hmem, hm, hg, hp, pd, hpdr, hpdt, alts, ?_, this⟩
  rw [← hm]; exact hcomp

end

end DiscoveryLogic
