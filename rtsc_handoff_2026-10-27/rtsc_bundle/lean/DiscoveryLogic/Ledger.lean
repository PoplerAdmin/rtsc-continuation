/-
DiscoveryLogic.Ledger: the verdict ledger of the discovery loop
(`/home/claude/rtsc/rtsc/discovery/ledger.py`) as a list of verdict events, the soundness guard
`append_verdict` and the current-status view.

Every definition carries the Python file:line it mirrors (`ledger.py` = `rtsc/discovery/ledger.py`,
`goals.py` = `rtsc/discovery/goals.py`, `status.py` = `rtsc/status.py`).

Abstraction of strings:
* `stage` (`ledger.py:51-58`, `ALL_STAGES`) is a closed type `Stage`; Python raises `ValueError`
  for any other string (`ledger.py:308-309`; the SQL trigger `ledger.py:117` likewise), so no other value can reach the table.
* `origin` is abstracted to `Origin`: `.classify` iff `origin == CLASSIFY_ORIGIN`
  (`= "rtsc.certify.classify"`, `pipeline.py:55`), `.stability` iff
  `origin.startswith("rtsc.stability")`, `.other` otherwise. The two predicates are disjoint
  (`"rtsc.certify.classify"` does not start with `"rtsc.stability"`), so this is a function of
  the string.
* Real numbers stand for the Python floats `T_star_K`, `pressure_GPa`, `P_max_GPa`.
-/
import PipelineLogic.Classify

namespace DiscoveryLogic

open PipelineLogic

noncomputable section

/-- Heuristic ledger stages, the keys of `HEURISTIC_STAGES` (`ledger.py:51-57`). The task's
names: generator = `generate`, budget = `select`, error = `executor`. -/
inductive HStage
  | generate
  | scope
  | config
  | prescreen
  | select
  | policy
  | executor
  deriving DecidableEq, Repr

/-- `ALL_STAGES = set(HEURISTIC_STAGES) | {"pipeline"}` (`ledger.py:58`). -/
inductive Stage
  | pipeline
  | heuristic (h : HStage)
  deriving DecidableEq, Repr

/-- Abstraction of the `origin` string (see the header). -/
inductive Origin
  | classify
  | stability
  | other
  deriving DecidableEq, Repr

/-- `HEURISTIC_STAGES[h]` (`ledger.py:51-57`): the statuses stage `h` may write. -/
def heuristicAllowed : HStage → Status → Bool
  | .generate, s => s == .deferred || s == .error
  | .scope, s => s == .deferred
  | .config, s => s == .deferred
  | .prescreen, s => s == .deferred
  | .select, s => s == .deferred
  | .policy, s => s == .deferred
  | .executor, s => s == .error || s == .deferred

/-- `PROVEN_STATUSES` (`status.py:27`). -/
def isProven (s : Status) : Bool := s == .rejectedProven || s == .certifiedCandidate

/-- A discovery goal `(id, T*, P_max)` (`goals.py:67-76`); ids are strings in Python, here `ℕ`. -/
structure Goal where
  id : ℕ
  Tstar : ℝ
  Pmax : ℝ

/-- `Goal.in_scope` (`goals.py:78-79`): `P ≤ P_max`. -/
def inScopeB (g : Goal) (P : ℝ) : Bool := decide (P ≤ g.Pmax)

theorem inScopeB_iff (g : Goal) (P : ℝ) : inScopeB g P = true ↔ P ≤ g.Pmax := by
  simp [inScopeB]

/-- One row of the `verdicts` table (`ledger.py:99-106`), restricted to the columns the
invariants are about (bracket, rule ids, flags, paths, digests are omitted). `round` is the
`round` column; `inScope` the `in_scope` column; `pressure` the material pressure
`pressure_GPa`; `dftPressure` the DFT pressure actually used, `dft_pressure_GPa` (NULL =
`none`); `dftHash`, `certifyHash` the config hashes (strings in Python, abstracted to `ℕ`;
NULL = `none`). -/
structure Event where
  mat : ℕ
  goal : ℕ
  round : ℕ
  stage : Stage
  origin : Origin
  status : Status
  Tstar : ℝ
  pressure : ℝ
  inScope : Bool
  dftPressure : Option ℝ
  dftHash : Option ℕ
  certifyHash : Option ℕ

/-- **The guard of `Ledger.append_verdict`, exactly as implemented** (`ledger.py:307-318`):
* `ledger.py:310`: a non-pipeline stage may only write `HEURISTIC_STAGES[stage]`;
* `ledger.py:312`: a proven status needs `stage == "pipeline"` and `origin == CLASSIFY_ORIGIN`;
* `ledger.py:314-316`: EXCLUDED_ASSUMPTION needs `stage == "pipeline"` and an
  `rtsc.stability*` origin;
* `ledger.py:317`: CERTIFIED_CANDIDATE needs `in_scope`.
The guard does NOT check `T_star_K` or `in_scope` against any goal: both are supplied by the
caller; whether a row counts for a goal is decided by `counts_for` (`DiscoveryLogic.Soundness`). -/
def appendAllowed (e : Event) : Bool :=
  (match e.stage with
    | .pipeline => true
    | .heuristic h => heuristicAllowed h e.status) &&
  (!isProven e.status || (e.stage == .pipeline && e.origin == .classify)) &&
  (e.status != .excludedAssumption || (e.stage == .pipeline && e.origin == .stability)) &&
  (e.status != .certifiedCandidate || e.inScope)

/-- The SQL trigger `verdicts_guard` (`ledger.py:115-126`), clause by clause:
`stage = pipeline OR status = DEFERRED OR (status = ERROR AND stage IN (generate, executor))`,
proven statuses need `origin = rtsc.certify.classify`, EXCLUDED_ASSUMPTION needs
`substr(origin,1,14) = 'rtsc.stability'`, CERTIFIED_CANDIDATE needs `in_scope = 1`.
(The `stage IN (...)` / `status IN (...)` clauses are the closed types.) -/
def sqlGuard (e : Event) : Bool :=
  (e.stage == .pipeline || e.status == .deferred ||
    (e.status == .error &&
      (e.stage == .heuristic .generate || e.stage == .heuristic .executor))) &&
  (!isProven e.status || e.origin == .classify) &&
  (e.status != .excludedAssumption || e.origin == .stability) &&
  (e.status != .certifiedCandidate || e.inScope)

/-- **The SQL trigger and the Python check agree on every row**, so a raw connection that
bypasses `append_verdict` cannot insert anything `append_verdict` would refuse. -/
theorem appendAllowed_eq_sqlGuard (e : Event) : appendAllowed e = sqlGuard e := by
  obtain ⟨m, g, r, st, o, s, T, P, b, dp, dh, ch⟩ := e
  cases st with
  | pipeline => cases s <;> cases o <;> cases b <;> rfl
  | heuristic h => cases h <;> cases s <;> cases o <;> cases b <;> rfl

/-- One SQLite transaction appending `es` (`ledger.py:319-331` inside `Ledger.tx`,
`ledger.py:160-175`): if any `append_verdict` raises, the outermost transaction rolls back, so
either all of `es` is appended or nothing. The only `INSERT INTO verdicts` in the code base is
`ledger.py:320-321`, and triggers `ledger.py:108-114` (with `recursive_triggers=ON`,
`ledger.py:149`) reject UPDATE, DELETE, INSERT OR REPLACE and re-used ids: the table only
grows. -/
def commit (L es : List Event) : List Event :=
  if es.all appendAllowed then L ++ es else L

/-- The `current_status` VIEW (`ledger.py:127-128`, `MAX(id)` per (material, goal)): the latest
event of `(m, g)`. `none` = no row. -/
def currentEvent (L : List Event) (m g : ℕ) : Option Event :=
  L.reverse.find? (fun e => e.mat == m && e.goal == g)

/-- Status of the current verdict (`Ledger.current`, `ledger.py:333-345`). -/
def currentStatus (L : List Event) (m g : ℕ) : Option Status :=
  (currentEvent L m g).map Event.status

/-! ## The guard -/

/-- **S2 for one event.** An event accepted from a heuristic stage is DEFERRED or ERROR. -/
theorem appendAllowed_heuristic {e : Event} (he : appendAllowed e = true) {h : HStage}
    (hs : e.stage = .heuristic h) : e.status = .deferred ∨ e.status = .error := by
  obtain ⟨m, g, r, st, o, s, T, P, b, dp, dh, ch⟩ := e
  simp only at hs
  subst hs
  cases h <;> cases s <;> simp_all [appendAllowed, heuristicAllowed, isProven]

/-- **S1 for one event.** A proven status is accepted only from `pipeline` /
`rtsc.certify.classify`. -/
theorem appendAllowed_proven {e : Event} (he : appendAllowed e = true)
    (hp : e.status = .rejectedProven ∨ e.status = .certifiedCandidate) :
    e.stage = .pipeline ∧ e.origin = .classify := by
  obtain ⟨m, g, r, st, o, s, T, P, b, dp, dh, ch⟩ := e
  simp only at hp ⊢
  rcases hp with rfl | rfl <;>
    simp_all [appendAllowed, isProven]

/-- **S3 for one event.** EXCLUDED_ASSUMPTION only from `pipeline` / `rtsc.stability*`. -/
theorem appendAllowed_excluded {e : Event} (he : appendAllowed e = true)
    (hx : e.status = .excludedAssumption) :
    e.stage = .pipeline ∧ e.origin = .stability := by
  obtain ⟨m, g, r, st, o, s, T, P, b, dp, dh, ch⟩ := e
  simp only at hx ⊢
  subst hx
  simp_all [appendAllowed, isProven]

/-- CERTIFIED_CANDIDATE is accepted only with `in_scope = 1`. -/
theorem appendAllowed_certified_inScope {e : Event} (he : appendAllowed e = true)
    (hc : e.status = .certifiedCandidate) : e.inScope = true := by
  obtain ⟨m, g, r, st, o, s, T, P, b, dp, dh, ch⟩ := e
  simp only at hc ⊢
  subst hc
  simp_all [appendAllowed, isProven]

/-! ## `commit` and `currentEvent` -/

theorem prefix_commit (L es : List Event) : L <+: commit L es := by
  unfold commit; split
  · exact List.prefix_append L es
  · exact List.prefix_refl L

theorem mem_commit {L es : List Event} {e : Event} (h : e ∈ commit L es) :
    e ∈ L ∨ (e ∈ es ∧ appendAllowed e = true) := by
  unfold commit at h
  split at h
  · rename_i hall
    rcases List.mem_append.1 h with h | h
    · exact Or.inl h
    · exact Or.inr ⟨h, List.all_eq_true.1 hall e h⟩
  · exact Or.inl h

theorem commit_of_all {L es : List Event} (h : ∀ e ∈ es, appendAllowed e = true) :
    commit L es = L ++ es := by
  unfold commit
  split
  · rfl
  · rename_i hn; exact absurd (List.all_eq_true.2 h) hn

theorem currentEvent_spec {L : List Event} {m g : ℕ} {e : Event}
    (h : currentEvent L m g = some e) : e ∈ L ∧ e.mat = m ∧ e.goal = g := by
  unfold currentEvent at h
  have h1 := List.mem_of_find?_eq_some h
  have h2 := List.find?_some h
  simp only [Bool.and_eq_true, beq_iff_eq] at h2
  exact ⟨List.mem_reverse.1 h1, h2.1, h2.2⟩

theorem currentEvent_isSome {L : List Event} {m g : ℕ} :
    (currentEvent L m g).isSome ↔ ∃ e ∈ L, e.mat = m ∧ e.goal = g := by
  unfold currentEvent
  simp [List.find?_isSome]

/-- The current event is the latest one: appending an event of `(m, g)` makes it current;
appending events of other pairs leaves the current event unchanged. -/
theorem currentEvent_append (L es : List Event) (m g : ℕ) :
    currentEvent (L ++ es) m g =
      (currentEvent es m g).or (currentEvent L m g) := by
  unfold currentEvent
  rw [List.reverse_append, List.find?_append]

theorem currentStatus_eq_some {L : List Event} {m g : ℕ} {s : Status} :
    currentStatus L m g = some s ↔ ∃ e, currentEvent L m g = some e ∧ e.status = s := by
  unfold currentStatus
  simp [Option.map_eq_some_iff]

end

end DiscoveryLogic
