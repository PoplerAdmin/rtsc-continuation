/-
DiscoveryLogic.Loop: the write paths of the discovery loop (`rtsc/discovery/loop.py`) as
operations on the ledger state, and the invariants every reachable state satisfies.

The operations are finer than (or as fine as) the Python transactions, so every state the
Python code can commit, including after an exception or a crash at any point, is `Reachable`:
* `genItem`: one generated structure (Python commits the whole generate phase in one
  transaction, `loop.py:411-465`);
* `prescreen`: one material (`loop.py:543-558`, one transaction per material);
* `addJob` / `selectDefer`: one row (`loop.py:591-607`, one transaction for the phase);
* `exec`: one (material, goal) pair (`loop.py:676-678` one transaction per goal group,
  `loop.py:689-690` one per pair);
* `startRound` (`ledger.py:375-379`, `loop.py:149-151`) and `setPhase` (`ledger.py:381-391`).
* `newSession` (`run` / `continue`, `loop.py:113-124, 128-133`): a possibly different config; new goal
  ids get initial statuses for every material (`_add_goal_statuses`); a missing old goal id is
  a `ConfigError` (no change);
* `reinit`, `configDefer`, `exec`: the writes of `_reconcile` (`loop.py:178-247`).
-/
import DiscoveryLogic.Ledger
import Mathlib.Data.List.Nodup

namespace DiscoveryLogic

open PipelineLogic

noncomputable section

/-- `ROUND_PHASES` (`ledger.py:59`). -/
inductive Phase
  | proposed
  | generated
  | prescreened
  | selected
  | executed
  | done
  deriving DecidableEq, Repr

/-! ## The events each write path appends -/

/-- A heuristic-stage row: origin is never `rtsc.certify.classify` / `rtsc.stability*`
(`discovery.generate`, `discovery.scope`, `discovery.config`, `discovery.prescreen:*`,
`discovery.select`, `discovery.executor`), no DFT pressure, no hashes. Python passes
`in_scope=g.in_scope(P)` or the default `True` only in branches where `g.in_scope(P)` holds, or
`False` where it fails; in every case the value is `inScopeB g P`. -/
def hEv (k m : ℕ) (P : ℝ) (g : Goal) (h : HStage) (st : Status) : Event :=
  { mat := m, goal := g.id, round := k, stage := .heuristic h, origin := .other, status := st,
    Tstar := g.Tstar, pressure := P, inScope := inScopeB g P, dftPressure := none,
    dftHash := none, certifyHash := none }

/-- `_initial_status_for` (`loop.py:366-384`): ERROR (stage `generate`) for an invalid
structure, DEFERRED (stage `generate`) if in scope, else DEFERRED (stage `scope`). -/
def initialEvent (k m : ℕ) (P : ℝ) (err : Bool) (g : Goal) : Event :=
  if err then hEv k m P g .generate .error
  else if inScopeB g P then hEv k m P g .generate .deferred
  else hEv k m P g .scope .deferred

/-- ML pre-screen deferral (`loop.py:550-557`), in-scope goals only. -/
def prescreenEvent (k m : ℕ) (P : ℝ) (g : Goal) : Event := hEv k m P g .prescreen .deferred

/-- Not selected within the budget (`loop.py:597-606`), in-scope goals only. -/
def selectEvent (k m : ℕ) (P : ℝ) (g : Goal) : Event := hEv k m P g .select .deferred

/-- Spectra computed under another `dft_hash`: DEFERRED, stage `config`
(`loop.py:205-218`), in-scope goals only. -/
def configEvent (k m : ℕ) (P : ℝ) (g : Goal) : Event := hEv k m P g .config .deferred

/-- No finished pipeline record (`loop.py:762-772`). -/
def execErrorEvent (k m : ℕ) (P : ℝ) (g : Goal) : Event := hEv k m P g .executor .error

/-- A proven verdict whose spectra were computed at an unknown pressure or one different from
the material's (`loop.py:779-789`): ERROR, stage `executor`, the verdict is not recorded. -/
def mismatchEvent (k m : ℕ) (P : ℝ) (g : Goal) (dp : Option ℝ) : Event :=
  { (hEv k m P g .executor .error) with dftPressure := dp }

/-- The pipeline verdict (`loop.py:811-818`): stage `pipeline`, the record's `verdict_origin`
and status, `T_star_K = g.T_star_K` (the pipeline config of `g` certifies at `g.T*`,
`goals.py:170`, and `loop.py:775-778` raises if a proven verdict was made at another T*),
`in_scope = g.in_scope(P)`, the DFT pressure of the record, and the CURRENT `certify_hash(g)`
and `dft_hash`. -/
def pipelineEvent (k m : ℕ) (P : ℝ) (g : Goal) (o : Origin) (st : Status) (dp : Option ℝ)
    (ch dh : ℕ) : Event :=
  { mat := m, goal := g.id, round := k, stage := .pipeline, origin := o, status := st,
    Tstar := g.Tstar, pressure := P, inScope := inScopeB g P, dftPressure := dp,
    dftHash := some dh, certifyHash := some ch }

/-- `1e-9` (GPa), the pressure tolerance of `loop.py:779` and `goals.py:227`. -/
def pTol : ℝ := 1 / 10 ^ 9

/-- `dp is not None and abs(dp - P) <= 1e-9` (`loop.py:779`). -/
def pOk : Option ℝ → ℝ → Bool
  | some x, P => decide (|x - P| ≤ pTol)
  | none, _ => false

/-! ## State -/

/-- The configuration of a session: enabled goals, `DiscoveryConfig.certify_hash(goal)` and
`dft_hash()` (`goals.py:174-186`), strings abstracted to `ℕ`. -/
structure Config where
  goals : List Goal
  certH : Goal → ℕ
  dftH : ℕ

/-- A `materials` row (`ledger.py:76-82`): `material_id` (here `ℕ`, the `seq`), pressure,
round of first generation, `valid`. -/
structure MatRow where
  id : ℕ
  P : ℝ
  round : ℕ
  valid : Bool

/-- One structure returned by a generator: the request's pressure and whether
`validate_structure` accepted it (`loop.py:413-430`). -/
structure GenItem where
  P : ℝ
  valid : Bool

/-- The ledger tables the invariants are about. `generated` counts processed generator outputs
(`counts["generated"]`); `dups` holds the material each `duplicates` row links to;
`screened` the `ml_screen` keys; `jobs` the `dft_jobs` (material, round) pairs. -/
structure Data where
  ledger : List Event
  mats : List MatRow
  dups : List ℕ
  generated : ℕ
  screened : List ℕ
  jobs : List (ℕ × ℕ)

/-- Ledger state, the session's config, and the `rounds` row of the latest round. -/
structure State where
  data : Data
  cfg : Config
  round : ℕ
  phase : Phase

/-- A fresh ledger with config `c`; no round yet (`last_round() is None`, treated like a
finished round 0, `loop.py:149`). -/
def State.empty (c : Config) : State := ⟨⟨[], [], [], 0, [], []⟩, c, 0, .done⟩

/-- `Ledger.material` (`ledger.py`). -/
def Data.mat? (d : Data) (m : ℕ) : Option MatRow := d.mats.find? (fun r => r.id == m)

/-- `seq = MAX(seq) + 1` (`ledger.py:201-203`, `loop.py:412`): a fresh material id. -/
def Data.fresh (d : Data) : ℕ := (d.mats.map MatRow.id).foldr max 0 + 1

/-- The "already recorded (resumed round)" test of `_record` (`loop.py:791-795`): the current
row is a pipeline row of the same record (record path = round directory; digest and config
hash abstracted to the round), same status, current `certify_hash` and `dft_hash`. -/
def alreadyRecorded (L : List Event) (k m g : ℕ) (st : Status) (ch dh : ℕ) : Bool :=
  match currentEvent L m g with
  | some e => e.stage == .pipeline && e.round == k && e.status == st &&
      e.certifyHash == some ch && e.dftHash == some dh
  | none => false

/-- Data-level operations (one per Python write path). -/
inductive DOp
  /-- one generator output; `dd` is the dedupe decision (`_match`, confirmed fingerprint or
  StructureMatcher hit; for an invalid structure the identical-raw-output hit), evaluated on
  the ledger at that moment -/
  | genItem (it : GenItem) (dd : Data → Option ℕ)
  | prescreen (m : ℕ) (defer : Bool)
  | addJob (m : ℕ)
  | selectDefer (m : ℕ)
  /-- `_record` for (m, g); `out = none`: no finished record; `some (origin, status, dft_p)` -/
  | exec (m : ℕ) (g : Goal) (out : Option (Origin × Status × Option ℝ))
  /-- `_reconcile`: `_initial_status_for(g, m, ...)` (scope changed) -/
  | reinit (m : ℕ) (g : Goal)
  /-- `_reconcile`: DEFERRED (stage `config`) for every in-scope goal -/
  | configDefer (m : ℕ)

/-- A new material with its initial status for every goal (`loop.py:443-459`; a hit `existing` is a duplicate, `loop.py:437-442`). -/
def newMat (c : Config) (k : ℕ) (it : GenItem) (d : Data) : Data :=
  { d with mats := d.mats ++ [⟨d.fresh, it.P, k, it.valid⟩], generated := d.generated + 1,
           ledger := commit d.ledger (c.goals.map (initialEvent k d.fresh it.P (!it.valid))) }

/-- Effect of one data operation in round `k` under config `c`. -/
def dataStep (c : Config) (k : ℕ) : DOp → Data → Data
  | .genItem it dd, d =>
    match dd d with
    | some x =>
      if d.mats.any (fun r => r.id == x) then
        { d with dups := d.dups ++ [x], generated := d.generated + 1 }
      else newMat c k it d
    | none => newMat c k it d
  | .prescreen m defer, d =>
    match d.mat? m with
    | none => d
    | some r =>
      { d with screened := d.screened ++ [m],
               ledger := if defer then
                   commit d.ledger ((c.goals.filter (fun g => inScopeB g r.P)).map
                     (prescreenEvent k m r.P))
                 else d.ledger }
  | .addJob m, d => { d with jobs := d.jobs ++ [(m, k)] }
  | .selectDefer m, d =>
    match d.mat? m with
    | none => d
    | some r =>
      { d with ledger := commit d.ledger ((c.goals.filter (fun g => inScopeB g r.P)).map
                 (selectEvent k m r.P)) }
  | .exec m g out, d =>
    match d.mat? m with
    | none => d
    | some r =>
      if inScopeB g r.P then
        match out with
        | none => { d with ledger := commit d.ledger [execErrorEvent k m r.P g] }
        | some (o, st, dp) =>
          if isProven st && !pOk dp r.P then
            { d with ledger := commit d.ledger [mismatchEvent k m r.P g dp] }
          else if alreadyRecorded d.ledger k m g.id st (c.certH g) c.dftH then d
          else
            let e := pipelineEvent k m r.P g o st dp (c.certH g) c.dftH
            { d with ledger := commit d.ledger [e] }
      else d
  | .reinit m g, d =>
    match d.mat? m with
    | none => d
    | some r => { d with ledger := commit d.ledger [initialEvent k m r.P false g] }
  | .configDefer m, d =>
    match d.mat? m with
    | none => d
    | some r =>
      { d with ledger := commit d.ledger ((c.goals.filter (fun g => inScopeB g r.P)).map
                 (configEvent k m r.P)) }

/-- `_add_goal_statuses` (`loop.py:386-391`, called from `loop.py:121-124`): for every material and every goal of the new
config whose id is not yet in the ledger, the initial status (`error=None`). -/
def addGoalEvents (old : List Goal) (c' : Config) (d : Data) : List Event :=
  d.mats.flatMap (fun r => (c'.goals.filter (fun g => !((old.map Goal.id).contains g.id))).map
    (initialEvent r.round r.id r.P false))

/-- All operations on the state. -/
inductive Op
  | startRound
  | setPhase (p : Phase)
  | data (o : DOp)
  | newSession (c' : Config)

def step : Op → State → State
  | .startRound, s => { s with round := s.round + 1, phase := .proposed }
  | .setPhase p, s => { s with phase := p }
  | .data o, s => { s with data := dataStep s.cfg s.round o s.data }
  | .newSession c', s =>
    if (s.cfg.goals.map Goal.id).all (fun x => (c'.goals.map Goal.id).contains x) then
      { s with cfg := c',
               data := { s.data with
                 ledger := (commit s.data.ledger (addGoalEvents s.cfg.goals c' s.data)) } }
    else s

/-- The loop only executes / re-initialises goals of its configuration. -/
def DOpOk (c : Config) : DOp → Prop
  | .exec _ g _ => g ∈ c.goals
  | .reinit _ g => g ∈ c.goals
  | _ => True

def OpOk (s : State) : Op → Prop
  | .data o => DOpOk s.cfg o
  | _ => True

/-- States reachable by ANY sequence of operations from an empty ledger: every `run`, every
crash (a prefix of the operations), every `continue` (with a possibly changed config), with
arbitrary generator, dedupe, pre-screen, policy and pipeline outputs. -/
inductive Reachable : State → Prop
  | init (c : Config) : Reachable (State.empty c)
  | step {s : State} (op : Op) (hop : OpOk s op) : Reachable s → Reachable (step op s)

/-! ## Append-only -/

/-- `L'` extends `L` by events the guard accepted. -/
def Ext (L L' : List Event) : Prop := ∃ es, L' = L ++ es ∧ ∀ e ∈ es, appendAllowed e = true

theorem Ext.refl (L : List Event) : Ext L L := ⟨[], by simp, by simp⟩

theorem Ext.trans {L1 L2 L3 : List Event} (h1 : Ext L1 L2) (h2 : Ext L2 L3) : Ext L1 L3 := by
  obtain ⟨a, rfl, ha⟩ := h1
  obtain ⟨b, rfl, hb⟩ := h2
  refine ⟨a ++ b, by simp, fun e he => ?_⟩
  rcases List.mem_append.1 he with h | h
  · exact ha e h
  · exact hb e h

theorem Ext.prefix {L L' : List Event} (h : Ext L L') : L <+: L' := by
  obtain ⟨a, rfl, -⟩ := h; exact List.prefix_append _ _

theorem ext_commit (L es : List Event) : Ext L (commit L es) := by
  unfold commit; split
  · rename_i hall; exact ⟨es, rfl, List.all_eq_true.1 hall⟩
  · exact Ext.refl L

theorem dataStep_ext (c : Config) (k : ℕ) (o : DOp) (d : Data) :
    Ext d.ledger (dataStep c k o d).ledger := by
  cases o <;> simp only [dataStep, newMat] <;> (repeat' split) <;> (try dsimp only) <;>
    first | exact Ext.refl _ | exact ext_commit _ _

theorem step_ext (op : Op) (s : State) : Ext s.data.ledger (step op s).data.ledger := by
  cases op with
  | startRound => exact Ext.refl _
  | setPhase p => exact Ext.refl _
  | data o => exact dataStep_ext s.cfg s.round o s.data
  | newSession c' =>
    simp only [step]; split
    · exact ext_commit _ _
    · exact Ext.refl _

/-- **Append-only (one operation).** Every operation extends the verdict list (by events the
guard accepted). -/
theorem step_prefix (op : Op) (s : State) : s.data.ledger <+: (step op s).data.ledger :=
  (step_ext op s).prefix

/-! ## Invariants -/

/-- The ledger invariant (what `Ledger.accounting` checks, plus the guard). -/
structure Inv (s : State) : Prop where
  good : ∀ e ∈ s.data.ledger, appendAllowed e = true
  covered : ∀ r ∈ s.data.mats, ∀ g ∈ s.cfg.goals, ∃ e ∈ s.data.ledger, e.mat = r.id ∧ e.goal = g.id
  count : s.data.generated = s.data.mats.length + s.data.dups.length
  dupsLinked : ∀ x ∈ s.data.dups, ∃ r ∈ s.data.mats, r.id = x
  nodup : (s.data.mats.map MatRow.id).Nodup

theorem inv_empty (c : Config) : Inv (State.empty c) := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> simp [State.empty]

theorem initialEvent_mat (k m : ℕ) (P : ℝ) (err : Bool) (g : Goal) :
    (initialEvent k m P err g).mat = m ∧ (initialEvent k m P err g).goal = g.id := by
  unfold initialEvent; split_ifs <;> exact ⟨rfl, rfl⟩

theorem initialEvent_allowed (k m : ℕ) (P : ℝ) (err : Bool) (g : Goal) :
    appendAllowed (initialEvent k m P err g) = true := by
  unfold initialEvent hEv; split_ifs <;> simp [appendAllowed, heuristicAllowed, isProven]

theorem good_ext {L L' : List Event} (hL : ∀ e ∈ L, appendAllowed e = true) (h : Ext L L') :
    ∀ e ∈ L', appendAllowed e = true := by
  obtain ⟨es, rfl, hes⟩ := h
  intro e he
  rcases List.mem_append.1 he with h | h
  · exact hL e h
  · exact hes e h

theorem le_foldr_max (l : List ℕ) : ∀ x ∈ l, x ≤ l.foldr max 0 := by
  induction l with
  | nil => simp
  | cons a l ih =>
    intro x hx
    simp only [List.foldr_cons]
    rcases List.mem_cons.1 hx with rfl | hx
    · exact le_max_left _ _
    · exact (ih x hx).trans (le_max_right _ _)

theorem fresh_not_mem (d : Data) : d.fresh ∉ d.mats.map MatRow.id := by
  intro h
  have := le_foldr_max _ _ h
  unfold Data.fresh at this
  omega

/-- Preservation when only the ledger (by `Ext`), `screened` or `jobs` change. -/
theorem inv_of_ext {s s' : State} (hs : Inv s) (hL : Ext s.data.ledger s'.data.ledger)
    (hm : s'.data.mats = s.data.mats) (hd : s'.data.dups = s.data.dups)
    (hg : s'.data.generated = s.data.generated) (hc : s'.cfg = s.cfg) : Inv s' := by
  refine ⟨good_ext hs.good hL, ?_, ?_, ?_, ?_⟩
  · intro r hr g hg'
    rw [hm] at hr; rw [hc] at hg'
    obtain ⟨e, he, h1⟩ := hs.covered r hr g hg'
    exact ⟨e, hL.prefix.subset he, h1⟩
  · rw [hm, hd, hg]; exact hs.count
  · rw [hm, hd]; exact hs.dupsLinked
  · rw [hm]; exact hs.nodup

theorem newMat_inv (c : Config) (k : ℕ) (it : GenItem) {s : State} (hs : Inv s)
    (hc : s.cfg = c) :
    Inv { s with data := newMat c k it s.data } := by
  subst hc
  have hall : ∀ e ∈ s.cfg.goals.map (initialEvent k s.data.fresh it.P (!it.valid)),
      appendAllowed e = true := by
    intro e he
    obtain ⟨g, _, rfl⟩ := List.mem_map.1 he
    exact initialEvent_allowed _ _ _ _ _
  refine ⟨good_ext hs.good (ext_commit _ _), ?_, ?_, ?_, ?_⟩
  · intro r hr g hg
    simp only [newMat] at hr ⊢
    rw [commit_of_all hall]
    rcases List.mem_append.1 hr with hr | hr
    · obtain ⟨e, he, h1⟩ := hs.covered r hr g hg
      exact ⟨e, List.mem_append_left _ he, h1⟩
    · simp only [List.mem_singleton] at hr
      subst hr
      exact ⟨_, List.mem_append_right _ (List.mem_map_of_mem hg), initialEvent_mat _ _ _ _ _⟩
  · simp [newMat, hs.count]; omega
  · intro x hx
    obtain ⟨r, hr, h1⟩ := hs.dupsLinked x hx
    exact ⟨r, List.mem_append_left _ hr, h1⟩
  · simp only [newMat, List.map_append, List.map_cons, List.map_nil]
    refine List.nodup_append.2 ⟨hs.nodup, List.nodup_singleton _, ?_⟩
    intro a ha b hb
    simp only [List.mem_singleton] at hb
    subst hb
    intro hab; subst hab
    exact fresh_not_mem _ ha

theorem dataStep_inv (o : DOp) (s : State) (hs : Inv s) :
    Inv { s with data := dataStep s.cfg s.round o s.data } := by
  have hext := dataStep_ext s.cfg s.round o s.data
  cases o with
  | genItem it dd =>
    simp only [dataStep] at hext ⊢
    split
    · split
      · rename_i x _ hany
        obtain ⟨r, hr, hrid⟩ := List.any_eq_true.1 hany
        refine ⟨hs.good, hs.covered, ?_, ?_, hs.nodup⟩
        · simp [hs.count]; omega
        · intro y hy
          rcases List.mem_append.1 hy with hy | hy
          · exact hs.dupsLinked y hy
          · simp only [List.mem_singleton] at hy
            exact ⟨r, hr, by simpa [hy] using hrid⟩
      · exact newMat_inv s.cfg s.round it hs rfl
    · exact newMat_inv s.cfg s.round it hs rfl
  | prescreen m defer =>
    refine inv_of_ext hs hext ?_ ?_ ?_ rfl <;> simp only [dataStep] <;> split <;> rfl
  | addJob m => exact inv_of_ext hs hext rfl rfl rfl rfl
  | selectDefer m =>
    refine inv_of_ext hs hext ?_ ?_ ?_ rfl <;> simp only [dataStep] <;> split <;> rfl
  | exec m g out =>
    refine inv_of_ext hs hext ?_ ?_ ?_ rfl <;> simp only [dataStep] <;> (repeat' split) <;> rfl
  | reinit m g =>
    refine inv_of_ext hs hext ?_ ?_ ?_ rfl <;> simp only [dataStep] <;> split <;> rfl
  | configDefer m =>
    refine inv_of_ext hs hext ?_ ?_ ?_ rfl <;> simp only [dataStep] <;> split <;> rfl

theorem newSession_inv (c' : Config) (s : State) (hs : Inv s) : Inv (step (.newSession c') s) := by
  simp only [step]
  split
  · rename_i hsub
    have hall : ∀ e ∈ addGoalEvents s.cfg.goals c' s.data, appendAllowed e = true := by
      intro e he
      simp only [addGoalEvents, List.mem_flatMap, List.mem_map] at he
      obtain ⟨r, -, g, -, rfl⟩ := he
      exact initialEvent_allowed _ _ _ _ _
    refine ⟨good_ext hs.good (ext_commit _ _), ?_, hs.count, hs.dupsLinked, hs.nodup⟩
    intro r hr g hg
    dsimp only at hr hg ⊢
    rw [commit_of_all hall]
    by_cases hold : (s.cfg.goals.map Goal.id).contains g.id = true
    · obtain ⟨g0, hg0, hid⟩ : ∃ g0 ∈ s.cfg.goals, g0.id = g.id := by
        simpa using hold
      obtain ⟨e, he, h1, h2⟩ := hs.covered r hr g0 hg0
      exact ⟨e, List.mem_append_left _ he, h1, h2.trans hid⟩
    · refine ⟨initialEvent r.round r.id r.P false g, List.mem_append_right _ ?_,
        initialEvent_mat _ _ _ _ _⟩
      simp only [addGoalEvents, List.mem_flatMap, List.mem_map, List.mem_filter]
      exact ⟨r, hr, g, ⟨hg, by simpa using hold⟩, rfl⟩
  · exact hs

theorem step_inv (op : Op) (s : State) (hs : Inv s) : Inv (step op s) := by
  cases op with
  | startRound => exact ⟨hs.good, hs.covered, hs.count, hs.dupsLinked, hs.nodup⟩
  | setPhase p => exact ⟨hs.good, hs.covered, hs.count, hs.dupsLinked, hs.nodup⟩
  | data o => exact dataStep_inv o s hs
  | newSession c' => exact newSession_inv c' s hs

theorem reachable_inv {s : State} (h : Reachable s) : Inv s := by
  induction h with
  | init c => exact inv_empty c
  | step op _ _ ih => exact step_inv op _ ih

/-! ## S2 (extended), S1, S3 on reachable ledgers -/

/-- **S2-ext.** In every reachable ledger, every row written by a heuristic stage
(generate / scope / config / prescreen / select / policy / executor) is DEFERRED or ERROR. -/
theorem S2_ext {s : State} (h : Reachable s) {e : Event} (he : e ∈ s.data.ledger) {hs : HStage}
    (hst : e.stage = .heuristic hs) : e.status = .deferred ∨ e.status = .error :=
  appendAllowed_heuristic ((reachable_inv h).good e he) hst

/-- **S2-ext, contrapositive.** No reachable ledger contains REJECTED_PROVEN,
CERTIFIED_CANDIDATE or EXCLUDED_ASSUMPTION from a heuristic stage. -/
theorem no_proven_from_heuristic {s : State} (h : Reachable s) {e : Event}
    (he : e ∈ s.data.ledger) {hs : HStage} (hst : e.stage = .heuristic hs) :
    e.status ≠ .rejectedProven ∧ e.status ≠ .certifiedCandidate ∧
      e.status ≠ .excludedAssumption := by
  rcases S2_ext h he hst with h1 | h1 <;> rw [h1] <;> decide

/-- **S1 on the ledger.** Proven statuses only from `pipeline` / `rtsc.certify.classify`. -/
theorem proven_only_from_classify {s : State} (h : Reachable s) {e : Event}
    (he : e ∈ s.data.ledger)
    (hp : e.status = .rejectedProven ∨ e.status = .certifiedCandidate) :
    e.stage = .pipeline ∧ e.origin = .classify :=
  appendAllowed_proven ((reachable_inv h).good e he) hp

/-- **S3 on the ledger.** EXCLUDED_ASSUMPTION only from `pipeline` / `rtsc.stability*`. -/
theorem excluded_only_from_stability {s : State} (h : Reachable s) {e : Event}
    (he : e ∈ s.data.ledger) (hx : e.status = .excludedAssumption) :
    e.stage = .pipeline ∧ e.origin = .stability :=
  appendAllowed_excluded ((reachable_inv h).good e he) hx

/-- A current CERTIFIED_CANDIDATE row of a reachable ledger comes from `certify.classify` and has
`in_scope = 1` (the guard; whether it COUNTS for a goal is `countsFor`, see
`DiscoveryLogic.Soundness`). -/
theorem current_certified_from_classify {s : State} (h : Reachable s) {m g : ℕ}
    (hc : currentStatus s.data.ledger m g = some .certifiedCandidate) :
    ∃ e, currentEvent s.data.ledger m g = some e ∧ e.stage = .pipeline ∧
      e.origin = .classify ∧ e.inScope = true := by
  obtain ⟨e, he, hst⟩ := currentStatus_eq_some.1 hc
  have hall := (reachable_inv h).good e (currentEvent_spec he).1
  obtain ⟨h1, h2⟩ := appendAllowed_proven hall (Or.inr hst)
  exact ⟨e, he, h1, h2, appendAllowed_certified_inScope hall hst⟩

/-- A proven row recorded by the loop has a DFT pressure equal to the material pressure up to
`1e-9` (`loop.py:779-789`), and the current `certify_hash` / `dft_hash` of the session that
wrote it. -/
theorem dataStep_proven_pressure (c : Config) (k : ℕ) (o : DOp) (d : Data) {e : Event}
    (he : e ∈ (dataStep c k o d).ledger) (hn : e ∉ d.ledger) (hp : isProven e.status = true) :
    ∃ x, e.dftPressure = some x ∧ |x - e.pressure| ≤ pTol := by
  cases o with
  | genItem it dd =>
    simp only [dataStep, newMat] at he
    (repeat' split at he) <;> (try dsimp only at he) <;>
    first
    | exact absurd he hn
    | (rcases mem_commit he with h | ⟨h, -⟩
       · exact absurd h hn
       · obtain ⟨g, -, rfl⟩ := List.mem_map.1 h
         unfold initialEvent hEv at hp; split_ifs at hp <;> simp [isProven] at hp)
  | prescreen m defer =>
    simp only [dataStep] at he
    (repeat' split at he) <;> (try dsimp only at he) <;>
    first
    | exact absurd he hn
    | (rcases mem_commit he with h | ⟨h, -⟩
       · exact absurd h hn
       · obtain ⟨g, -, rfl⟩ := List.mem_map.1 h; simp [prescreenEvent, hEv, isProven] at hp)
  | addJob m => exact absurd he hn
  | selectDefer m =>
    simp only [dataStep] at he
    (repeat' split at he) <;> (try dsimp only at he) <;>
    first
    | exact absurd he hn
    | (rcases mem_commit he with h | ⟨h, -⟩
       · exact absurd h hn
       · obtain ⟨g, -, rfl⟩ := List.mem_map.1 h; simp [selectEvent, hEv, isProven] at hp)
  | configDefer m =>
    simp only [dataStep] at he
    (repeat' split at he) <;> (try dsimp only at he) <;>
    first
    | exact absurd he hn
    | (rcases mem_commit he with h | ⟨h, -⟩
       · exact absurd h hn
       · obtain ⟨g, -, rfl⟩ := List.mem_map.1 h; simp [configEvent, hEv, isProven] at hp)
  | reinit m g =>
    simp only [dataStep] at he
    (repeat' split at he) <;> (try dsimp only at he) <;>
    first
    | exact absurd he hn
    | (rcases mem_commit he with h | ⟨h, -⟩
       · exact absurd h hn
       · simp only [List.mem_singleton] at h; subst h
         unfold initialEvent hEv at hp; split_ifs at hp <;> simp [isProven] at hp)
  | exec m g out =>
    simp only [dataStep] at he
    split at he
    · exact absurd he hn
    · rename_i r _
      split at he
      · split at he
        · dsimp only at he
          rcases mem_commit he with h | ⟨h, -⟩
          · exact absurd h hn
          · simp only [List.mem_singleton] at h; subst h
            simp [execErrorEvent, hEv, isProven] at hp
        · rename_i o st dp
          split at he
          · dsimp only at he
            rcases mem_commit he with h | ⟨h, -⟩
            · exact absurd h hn
            · simp only [List.mem_singleton] at h; subst h
              simp [mismatchEvent, hEv, isProven] at hp
          · rename_i hok
            split at he
            · exact absurd he hn
            · dsimp only at he
              rcases mem_commit he with h | ⟨h, -⟩
              · exact absurd h hn
              · simp only [List.mem_singleton] at h; subst h
                simp only [pipelineEvent] at hp ⊢
                simp only [hp, Bool.true_and, Bool.not_eq_true'] at hok
                cases dp with
                | none => simp [pOk] at hok
                | some x => exact ⟨x, rfl, by simpa [pOk] using hok⟩
      · exact absurd he hn

/-- On every reachable ledger, every REJECTED_PROVEN / CERTIFIED_CANDIDATE row has a known DFT
pressure within `1e-9` GPa of the material pressure (`loop.py:779-789`, audit F3). -/
theorem reachable_proven_pressure {s : State} (h : Reachable s) {e : Event}
    (he : e ∈ s.data.ledger) (hp : isProven e.status = true) :
    ∃ x, e.dftPressure = some x ∧ |x - e.pressure| ≤ pTol := by
  induction h with
  | init c => simp [State.empty] at he
  | step op _ _ ih =>
    rename_i s _ _
    by_cases hold : e ∈ s.data.ledger
    · exact ih hold
    · cases op with
      | startRound => exact absurd he hold
      | setPhase p => exact absurd he hold
      | data o => exact dataStep_proven_pressure s.cfg s.round o s.data he hold hp
      | newSession c' =>
        simp only [step] at he
        split at he
        · dsimp only at he
          rcases mem_commit he with h' | ⟨h', -⟩
          · exact absurd h' hold
          · simp only [addGoalEvents, List.mem_flatMap, List.mem_map] at h'
            obtain ⟨r, -, g, -, rfl⟩ := h'
            unfold initialEvent hEv at hp; split_ifs at hp <;> simp [isProven] at hp
        · exact absurd he hold

/-! ## Accounting (S4 for the ledger) -/

/-- **Accounting.** In every reachable state: every registered material has exactly one current
status for every goal of the current config (a value of the function `currentStatus`, the
latest row; `none` never occurs); every generated structure is either one material or one
duplicate (`generated = materials + duplicates`); every duplicate links to an existing
material; material ids are unique. -/
theorem accounting {s : State} (h : Reachable s) :
    (∀ r ∈ s.data.mats, ∀ g ∈ s.cfg.goals,
      ∃! st, currentStatus s.data.ledger r.id g.id = some st) ∧
    s.data.generated = s.data.mats.length + s.data.dups.length ∧
    (∀ x ∈ s.data.dups, ∃ r ∈ s.data.mats, r.id = x) ∧
    (s.data.mats.map MatRow.id).Nodup := by
  have hi := reachable_inv h
  refine ⟨?_, hi.count, hi.dupsLinked, hi.nodup⟩
  intro r hr g hg
  obtain ⟨e, he, h1, h2⟩ := hi.covered r hr g hg
  have hs : (currentEvent s.data.ledger r.id g.id).isSome :=
    currentEvent_isSome.2 ⟨e, he, h1, h2⟩
  obtain ⟨e', he'⟩ := Option.isSome_iff_exists.1 hs
  refine ⟨e'.status, ?_, ?_⟩
  · simp [currentStatus, he']
  · intro st hst
    simp [currentStatus, he'] at hst
    exact hst.symm

end

end DiscoveryLogic
