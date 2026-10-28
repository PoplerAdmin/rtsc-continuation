/-
DiscoveryLogic.Continue: round phases (`loop.py` `_run_round`) and `continue` after a crash.

Determinism assumptions are explicit (`EnvAgree`): the original run and the continued run see
the same generator outputs for the stored requests of round `k`, the same pre-screen decisions,
the same selection (policy deterministic given seed, round and ledger; same DFT budget) and the
same pipeline outcomes (pipeline resume reuses finished records). Nothing about the RNG is
proved.
-/
import DiscoveryLogic.Loop

namespace DiscoveryLogic

open PipelineLogic

noncomputable section

/-- The outside world of one run: generators, ML pre-screen, search policy, pipeline. -/
structure Env where
  /-- Generator outputs for the requests stored with round `k` (`loop.py:393-410`; the
  requests are read back from the `rounds` row, so a resumed round does not re-propose). -/
  gen : ℕ → List GenItem
  /-- Dedupe decision for a generator output in round `k`, given the ledger at that moment
  (`_match`, `loop.py:347-360`). -/
  dedupe : ℕ → GenItem → Data → Option ℕ
  /-- `_ml_decision` / relaxed-duplicate deferral for material `m` in round `k`
  (`loop.py:466-481, 521-542`). -/
  pre : ℕ → ℕ → Bool
  /-- `_select`'s chosen materials (`loop.py:572-590`): policy ranking within the budget. -/
  select : ℕ → State → List ℕ
  /-- Pipeline outcome for (round, material, goal id) (`loop.py:650-690, 758-818`); `none` = no
  finished record. -/
  outcome : ℕ → ℕ → ℕ → Option (Origin × Status × Option ℝ)

/-- **Determinism hypotheses** for round `k` between the original run `e1` and the
continued run `e2`. -/
structure EnvAgree (k : ℕ) (e1 e2 : Env) : Prop where
  gen : e2.gen k = e1.gen k
  dedupe : ∀ it, e2.dedupe k it = e1.dedupe k it
  pre : ∀ m, e2.pre k m = e1.pre k m
  select : ∀ s, e2.select k s = e1.select k s
  outcome : ∀ m g, e2.outcome k m g = e1.outcome k m g


/-- A sequence of data operations in round `k`. -/
def runData (c : Config) (k : ℕ) (ds : List DOp) (d : Data) : Data :=
  ds.foldl (fun d o => dataStep c k o d) d

/-- Valid materials first generated in round `k` (`ledger.materials(k)` with `valid`,
`loop.py:494`, `loop.py:595`). -/
def roundMats (d : Data) (k : ℕ) : List ℕ :=
  (d.mats.filter (fun r => r.round == k && r.valid)).map MatRow.id

/-- Materials of round `k` without an `ml_screen` row (`loop.py:494-495`). -/
def toScreen (d : Data) (k : ℕ) : List ℕ :=
  (roundMats d k).filter (fun m => !(d.screened.contains m))

def genOps (env : Env) (k : ℕ) : List DOp :=
  (env.gen k).map (fun it => DOp.genItem it (env.dedupe k it))

def preOps (env : Env) (k : ℕ) (d : Data) : List DOp :=
  (toScreen d k).map (fun m => DOp.prescreen m (env.pre k m))

/-- `_select` (`loop.py:589-607`): jobs for the chosen, DEFERRED for the other new ones. -/
def selOps (env : Env) (s : State) : List DOp :=
  (env.select s.round s).map DOp.addJob ++
    ((roundMats s.data s.round).filter (fun m => !((env.select s.round s).contains m))).map
      DOp.selectDefer

/-- `_execute` (`loop.py:650-697`): every job of round `k`, every in-scope goal. -/
def execOps (env : Env) (c : Config) (k : ℕ) (d : Data) : List DOp :=
  ((d.jobs.filter (fun j => j.2 == k)).map Prod.fst).flatMap
    (fun m => c.goals.map (fun g => DOp.exec m g (env.outcome k m g.id)))

/-- `start_round` (`loop.py:149-151`, `ledger.py:375-379`). -/
def startRound (s : State) : State := ⟨s.data, s.cfg, s.round + 1, .proposed⟩

/-- Phase `generated` (`loop.py:320-322`). -/
def runGen (env : Env) (s : State) : State :=
  ⟨runData s.cfg s.round (genOps env s.round) s.data, s.cfg, s.round, .generated⟩

/-- Phase `prescreened` (`loop.py:323-325`). -/
def runPre (env : Env) (s : State) : State :=
  ⟨runData s.cfg s.round (preOps env s.round s.data) s.data, s.cfg, s.round, .prescreened⟩

/-- Phase `selected` (`loop.py:326-328`). -/
def runSel (env : Env) (s : State) : State :=
  ⟨runData s.cfg s.round (selOps env s) s.data, s.cfg, s.round, .selected⟩

/-- Phase `executed` (`loop.py:329-331`). -/
def runExec (env : Env) (s : State) : State :=
  ⟨runData s.cfg s.round (execOps env s.cfg s.round s.data) s.data, s.cfg, s.round, .executed⟩

/-- Phase `done` (`loop.py:332-336`; policy snapshot and reports write no verdicts). -/
def runDone (s : State) : State := ⟨s.data, s.cfg, s.round, .done⟩

/-- Complete a round from the recorded phase (`_run_round`, `loop.py:318-336`). -/
def finish (env : Env) : Phase → State → State
  | .proposed, s => runDone (runExec env (runSel env (runPre env
      (runGen env s))))
  | .generated, s => runDone (runExec env (runSel env (runPre env s)))
  | .prescreened, s => runDone (runExec env (runSel env s))
  | .selected, s => runDone (runExec env s)
  | .executed, s => runDone s
  | .done, s => s

/-- One iteration of `run` (`loop.py:144-157`): resume the unfinished round at its phase, or
start the next one. `run` and `continue` both execute this. -/
def resume (env : Env) (s : State) : State :=
  if s.phase = .done then finish env .proposed (startRound s) else finish env s.phase s

/-- The committed states of one uninterrupted round from `s0` before it is done: after
`start_round` (0) and after phases generated..executed (1-4; `j ≥ 4` all mean 4). A crash during phase `j+1` leaves exactly snapshot `j` (each phase is
committed atomically in this model; see `continue_prescreen_partial` for the finer commits). -/
def snap (env : Env) (s0 : State) : ℕ → State
  | 0 => startRound s0
  | 1 => runGen env (startRound s0)
  | 2 => runPre env (runGen env (startRound s0))
  | 3 => runSel env (runPre env (runGen env (startRound s0)))
  | _ => runExec env (runSel env (runPre env (runGen env (startRound s0))))

/-! ## Phases are sequences of operations -/

theorem runData_reachable {ds : List DOp} {s : State}
    (hok : ∀ o ∈ ds, DOpOk s.cfg o) (h : Reachable s) :
    Reachable ⟨runData s.cfg s.round ds s.data, s.cfg, s.round, s.phase⟩ := by
  induction ds generalizing s with
  | nil => exact h
  | cons o ds ih =>
    have h1 := Reachable.step (.data o) (hok o (List.mem_cons_self ..)) h
    exact ih (s := step (.data o) s) (fun o' ho' => hok o' (List.mem_cons_of_mem _ ho')) h1

theorem runData_prefix (c : Config) (k : ℕ) (ds : List DOp) (d : Data) :
    d.ledger <+: (runData c k ds d).ledger := by
  induction ds generalizing d with
  | nil => exact List.prefix_refl _
  | cons o ds ih => exact (dataStep_ext c k o d).prefix.trans (ih _)

theorem execOps_ok (env : Env) (c : Config) (k : ℕ) (d : Data) :
    ∀ o ∈ execOps env c k d, DOpOk c o := by
  intro o ho
  simp only [execOps, List.mem_flatMap, List.mem_map] at ho
  obtain ⟨m, -, g, hg, rfl⟩ := ho
  exact hg

theorem setPhase_reachable {s : State} (p : Phase) (h : Reachable s) :
    Reachable ⟨s.data, s.cfg, s.round, p⟩ :=
  Reachable.step (.setPhase p) trivial h

theorem runGen_reachable (env : Env) {s : State} (h : Reachable s) : Reachable (runGen env s) :=
  setPhase_reachable .generated (runData_reachable (by simp [genOps, DOpOk]) h)

theorem runPre_reachable (env : Env) {s : State} (h : Reachable s) : Reachable (runPre env s) :=
  setPhase_reachable .prescreened (runData_reachable (by simp [preOps, DOpOk]) h)

theorem selOps_ok (env : Env) (s : State) : ∀ o ∈ selOps env s, DOpOk s.cfg o := by
  intro o ho
  simp only [selOps, List.mem_append, List.mem_map] at ho
  rcases ho with (⟨_, _, rfl⟩ | ⟨_, _, rfl⟩) <;> trivial

theorem runSel_reachable (env : Env) {s : State} (h : Reachable s) : Reachable (runSel env s) :=
  setPhase_reachable .selected (runData_reachable (selOps_ok env s) h)

theorem runExec_reachable (env : Env) {s : State} (h : Reachable s) : Reachable (runExec env s) :=
  setPhase_reachable .executed (runData_reachable (execOps_ok _ _ _ _) h)

theorem finish_reachable (env : Env) (p : Phase) {s : State} (h : Reachable s) :
    Reachable (finish env p s) := by
  have hD : ∀ {s}, Reachable s → Reachable (runDone s) :=
    fun h => setPhase_reachable .done h
  cases p <;> simp only [finish]
  · exact hD (runExec_reachable env (runSel_reachable env
      (runPre_reachable env (runGen_reachable env h))))
  · exact hD (runExec_reachable env (runSel_reachable env
      (runPre_reachable env h)))
  · exact hD (runExec_reachable env (runSel_reachable env h))
  · exact hD (runExec_reachable env h)
  · exact hD h
  · exact h

theorem resume_reachable (env : Env) {s : State} (h : Reachable s) :
    Reachable (resume env s) := by
  unfold resume
  split
  · exact finish_reachable env _ (Reachable.step .startRound trivial h)
  · exact finish_reachable env _ h

theorem finish_prefix (env : Env) (p : Phase) (s : State) :
    s.data.ledger <+: (finish env p s).data.ledger := by
  have hG : ∀ s, s.data.ledger <+: (runGen env s).data.ledger :=
    fun s => runData_prefix _ _ _ _
  have hP : ∀ s, s.data.ledger <+: (runPre env s).data.ledger :=
    fun s => runData_prefix _ _ _ _
  have hS : ∀ s, s.data.ledger <+: (runSel env s).data.ledger :=
    fun s => runData_prefix _ _ _ _
  have hE : ∀ s, s.data.ledger <+: (runExec env s).data.ledger :=
    fun s => runData_prefix _ _ _ _
  cases p <;> simp only [finish, runDone]
  · exact (hG _).trans ((hP _).trans ((hS _).trans (hE _)))
  · exact (hP _).trans ((hS _).trans (hE _))
  · exact (hS _).trans (hE _)
  · exact hE _
  · exact List.prefix_refl _
  · exact List.prefix_refl _

/-- **Append-only (rounds).** A resumed / completed round extends the ledger. -/
theorem resume_prefix (env : Env) (s : State) :
    s.data.ledger <+: (resume env s).data.ledger := by
  unfold resume
  split
  · exact finish_prefix env _ (startRound s)
  · exact finish_prefix env _ s

/-! ## Continuation after a crash at a phase boundary -/

theorem runGen_congr {k : ℕ} {e1 e2 : Env} (h : EnvAgree k e1 e2) (s : State)
    (hk : s.round = k) : runGen e2 s = runGen e1 s := by
  subst hk; simp only [runGen, genOps, h.gen, h.dedupe]

theorem runPre_congr {k : ℕ} {e1 e2 : Env} (h : EnvAgree k e1 e2) (s : State)
    (hk : s.round = k) : runPre e2 s = runPre e1 s := by
  subst hk; simp only [runPre, preOps, h.pre]

theorem runSel_congr {k : ℕ} {e1 e2 : Env} (h : EnvAgree k e1 e2) (s : State)
    (hk : s.round = k) : runSel e2 s = runSel e1 s := by
  subst hk; simp only [runSel, selOps, h.select]

theorem runExec_congr {k : ℕ} {e1 e2 : Env} (h : EnvAgree k e1 e2) (s : State)
    (hk : s.round = k) : runExec e2 s = runExec e1 s := by
  subst hk; simp only [runExec, execOps, h.outcome]

/-- **Continuation.** Let `s0` be a state between rounds (`phase = done`) and `e1` the original
run. If the process crashes after committing any phase `j` of the next round (`snap e1 s0 j`,
`j = 0..4`), or before `start_round` committed (`s0` itself), then `continue` under
deterministic outside behaviour (`EnvAgree (s0.round + 1) e1 e2`) produces exactly the state
of the uninterrupted round. -/
theorem continue_after_crash (e1 e2 : Env) (s0 : State) (h0 : s0.phase = .done)
    (h : EnvAgree (s0.round + 1) e1 e2) (j : ℕ) :
    resume e2 (snap e1 s0 j) = resume e1 s0 ∧
      resume e2 s0 = resume e1 s0 := by
  have hr : ∀ e, resume e s0 = runDone (runExec e (runSel e (runPre e
      (runGen e (startRound s0))))) := by
    intro e; simp only [resume, h0, ↓reduceIte]; rfl
  have k0 : (startRound s0).round = s0.round + 1 := rfl
  have hG := runGen_congr h (startRound s0) k0
  have hP := runPre_congr h (runGen e1 (startRound s0)) k0
  have hS := runSel_congr h (runPre e1 (runGen e1 (startRound s0))) k0
  have hE := runExec_congr h
    (runSel e1 (runPre e1 (runGen e1 (startRound s0)))) k0
  have hfull : resume e2 s0 = resume e1 s0 := by
    rw [hr e2, hr e1, hG, hP, hS, hE]
  refine ⟨?_, hfull⟩
  match j with
  | 0 =>
    show finish e2 .proposed (startRound s0) = _
    rw [hr e1]; simp only [finish]; rw [hG, hP, hS, hE]
  | 1 =>
    show finish e2 .generated (runGen e1 (startRound s0)) = _
    rw [hr e1]; simp only [finish]; rw [hP, hS, hE]
  | 2 =>
    show finish e2 .prescreened (runPre e1 (runGen e1 (startRound s0))) = _
    rw [hr e1]; simp only [finish]; rw [hS, hE]
  | 3 =>
    show finish e2 .selected
      (runSel e1 (runPre e1 (runGen e1 (startRound s0)))) = _
    rw [hr e1]; simp only [finish]; rw [hE]
  | _ + 4 =>
    show finish e2 .executed
      (runExec e1 (runSel e1 (runPre e1 (runGen e1 (startRound s0))))) = _
    rw [hr e1]; rfl

/-- Crash snapshots are truncations of the final ledger. -/
theorem snap_prefix (env : Env) (s0 : State) (h0 : s0.phase = .done) (j : ℕ) :
    (snap env s0 j).data.ledger <+: (resume env s0).data.ledger := by
  have hr : resume env s0 = runDone (runExec env (runSel env (runPre env
      (runGen env (startRound s0))))) := by
    simp only [resume, h0, ↓reduceIte]; rfl
  rw [hr]
  have hG : ∀ s, s.data.ledger <+: (runGen env s).data.ledger :=
    fun s => runData_prefix _ _ _ _
  have hP : ∀ s, s.data.ledger <+: (runPre env s).data.ledger :=
    fun s => runData_prefix _ _ _ _
  have hS : ∀ s, s.data.ledger <+: (runSel env s).data.ledger :=
    fun s => runData_prefix _ _ _ _
  have hE : ∀ s, s.data.ledger <+: (runExec env s).data.ledger :=
    fun s => runData_prefix _ _ _ _
  match j with
  | 0 => exact (hG _).trans ((hP _).trans ((hS _).trans (hE _)))
  | 1 => exact (hP _).trans ((hS _).trans (hE _))
  | 2 => exact (hS _).trans (hE _)
  | 3 => exact hE _
  | _ + 4 => exact List.prefix_refl _

/-- **Continuation, current statuses.** After a crash at any phase boundary and `continue`,
every material/goal pair has the same current status as in the uninterrupted run. -/
theorem continue_currentStatus (e1 e2 : Env) (s0 : State) (h0 : s0.phase = .done)
    (h : EnvAgree (s0.round + 1) e1 e2) (j m g : ℕ) :
    currentStatus (resume e2 (snap e1 s0 j)).data.ledger m g =
      currentStatus (resume e1 s0).data.ledger m g := by
  rw [(continue_after_crash e1 e2 s0 h0 h j).1]

/-- `continue` with the SAME config: the session start (`_add_goal_statuses`) adds nothing. -/
theorem newSession_same (s : State) : step (.newSession s.cfg) s = s := by
  have h1 : ((s.cfg.goals.map Goal.id).all fun x => (s.cfg.goals.map Goal.id).contains x) =
      true := by
    simp only [List.all_eq_true, List.contains_iff_mem]; exact fun x hx => hx
  have h2 : addGoalEvents s.cfg.goals s.cfg s.data = [] := by
    simp only [addGoalEvents]
    simp only [List.flatMap_eq_nil_iff, List.map_eq_nil_iff, List.filter_eq_nil_iff]
    intro r _ g hg
    simpa using List.mem_map_of_mem (f := Goal.id) hg
  simp only [step, h1, h2, ↓reduceIte, commit, List.all_nil, List.append_nil]

/-- **Continuation including the session start.** `continue` first starts a session with the
same config (no new goal statuses, `newSession_same`) and runs `_reconcile`, modelled as an
arbitrary state transformation `R`. Under the explicit hypothesis that `R` writes nothing at the
crash point (true for `_reconcile` with an unchanged config, except for the re-evaluation of
DFT-pressure-mismatch rows, see the correspondence notes), `continue` reproduces the
uninterrupted round. -/
theorem continue_with_session (e1 e2 : Env) (R : State → State) (s0 : State)
    (h0 : s0.phase = .done) (h : EnvAgree (s0.round + 1) e1 e2) (j : ℕ)
    (hR : R (snap e1 s0 j) = snap e1 s0 j) :
    resume e2 (R (step (.newSession (snap e1 s0 j).cfg) (snap e1 s0 j))) = resume e1 s0 := by
  rw [newSession_same, hR]
  exact (continue_after_crash e1 e2 s0 h0 h j).1

/-! ## Continuation after a crash inside the pre-screen phase

`_prescreen` commits one transaction per material (`loop.py:543-558`) and, when resumed, only
screens materials without an `ml_screen` row (`loop.py:494-495`). A crash after `i` materials
leaves the phase marker at `generated` with `i` materials screened. -/

theorem prescreen_ops_effect (c : Config) (k : ℕ) (f : ℕ → Bool) (ms : List ℕ) (d : Data)
    (hms : ∀ m ∈ ms, ∃ r ∈ d.mats, r.id = m) :
    (runData c k (ms.map fun m => DOp.prescreen m (f m)) d).mats = d.mats ∧
      (runData c k (ms.map fun m => DOp.prescreen m (f m)) d).screened =
        d.screened ++ ms := by
  induction ms generalizing d with
  | nil => simp [runData]
  | cons m ms ih =>
    obtain ⟨r0, hr0, hid⟩ := hms m (List.mem_cons_self ..)
    have hsome : (d.mat? m).isSome := by
      unfold Data.mat?; rw [List.find?_isSome]; exact ⟨r0, hr0, by simp [hid]⟩
    obtain ⟨r, hr⟩ := Option.isSome_iff_exists.1 hsome
    have hstep : (dataStep c k (DOp.prescreen m (f m)) d).mats = d.mats ∧
        (dataStep c k (DOp.prescreen m (f m)) d).screened = d.screened ++ [m] := by
      simp [dataStep, hr]
    have := ih (dataStep c k (DOp.prescreen m (f m)) d)
      (fun m' hm' => by rw [hstep.1]; exact hms m' (List.mem_cons_of_mem _ hm'))
    simp only [List.map_cons, runData, List.foldl_cons] at this ⊢
    refine ⟨this.1.trans hstep.1, ?_⟩
    rw [this.2, hstep.2]; simp

theorem filter_not_take (X : List ℕ) (hX : X.Nodup) (i : ℕ) :
    X.filter (fun m => !((X.take i).contains m)) = X.drop i := by
  have hd : (X.take i ++ X.drop i).Nodup := by rwa [List.take_append_drop]
  obtain ⟨-, -, hdisj⟩ := List.nodup_append.1 hd
  have e : X.filter (fun m => !((X.take i).contains m)) =
      (X.take i ++ X.drop i).filter (fun m => !((X.take i).contains m)) := by
    rw [List.take_append_drop]
  rw [e, List.filter_append]
  have h1 : (X.take i).filter (fun m => !((X.take i).contains m)) = [] :=
    List.filter_eq_nil_iff.2 (fun a ha => by simp [ha])
  have h2 : (X.drop i).filter (fun m => !((X.take i).contains m)) = X.drop i :=
    List.filter_eq_self.2 (fun a ha => by
      have : a ∉ X.take i := fun h' => hdisj a h' a ha rfl
      simpa using this)
  rw [h1, h2, List.nil_append]

theorem toScreen_nodup {s : State} (hd : Inv s) (k : ℕ) :
    (toScreen s.data k).Nodup := by
  unfold toScreen roundMats
  exact ((hd.nodup.sublist (List.filter_sublist.map _))).filter _

theorem toScreen_mem {d : Data} {k m : ℕ} (h : m ∈ toScreen d k) : ∃ r ∈ d.mats, r.id = m := by
  unfold toScreen roundMats at h
  obtain ⟨h1, -⟩ := List.mem_filter.1 h
  obtain ⟨r, hr, rfl⟩ := List.mem_map.1 h1
  exact ⟨r, (List.mem_filter.1 hr).1, rfl⟩

/-- The state after a crash that committed the first `i` pre-screen transactions. -/
def preCrash (env : Env) (s : State) (i : ℕ) : State :=
  ⟨runData s.cfg s.round ((preOps env s.round s.data).take i) s.data, s.cfg, s.round, s.phase⟩

theorem preCrash_reachable (env : Env) {s : State} (hs : Reachable s) (i : ℕ) :
    Reachable (preCrash env s i) :=
  runData_reachable (fun o ho => by
    obtain ⟨m, -, rfl⟩ := List.mem_map.1 (List.mem_of_mem_take ho); trivial) hs

/-- **Continuation inside the pre-screen phase.** If the process crashes after any number `i`
of per-material pre-screen commits, re-running the phase (which skips screened materials)
under deterministic pre-screen decisions yields exactly the uninterrupted phase result. -/
theorem continue_prescreen_partial (e1 e2 : Env) {s : State} (hs : Reachable s)
    (h : EnvAgree s.round e1 e2) (i : ℕ) :
    runPre e2 (preCrash e1 s i) = runPre e1 s := by
  have hinv := reachable_inv hs
  set c := s.cfg with hcc
  set k := s.round with hk
  set d := s.data with hdd
  set X := toScreen d k with hX
  have hXm : ∀ m ∈ X.take i, ∃ r ∈ d.mats, r.id = m :=
    fun m hm => toScreen_mem (List.mem_of_mem_take hm)
  have hpre1 : preOps e1 k d = X.map (fun m => DOp.prescreen m (e1.pre k m)) := rfl
  have heff := prescreen_ops_effect c k (e1.pre k) (X.take i) d hXm
  have hc : (preCrash e1 s i).data =
      runData c k ((X.take i).map fun m => DOp.prescreen m (e1.pre k m)) d := by
    simp only [preCrash, ← hcc, ← hk, ← hdd, hpre1, List.map_take]
  have hscreen : toScreen (preCrash e1 s i).data k = X.drop i := by
    rw [hc]
    unfold toScreen roundMats
    rw [heff.1, heff.2]
    rw [← filter_not_take X (toScreen_nodup hinv k) i]
    simp only [hX, toScreen, roundMats, List.filter_filter]
    apply List.filter_congr
    intro x _
    simp [Bool.and_comm]
  simp only [runPre, preCrash]
  congr 1
  show runData c k (preOps e2 k (preCrash e1 s i).data) (preCrash e1 s i).data =
    runData c k (preOps e1 k d) d
  rw [show preOps e2 k (preCrash e1 s i).data =
      (X.drop i).map (fun m => DOp.prescreen m (e1.pre k m)) by
    simp only [preOps, hscreen, h.pre]]
  rw [hc, hpre1]
  unfold runData
  rw [← List.foldl_append, ← List.map_append, List.take_append_drop]

end

end DiscoveryLogic
