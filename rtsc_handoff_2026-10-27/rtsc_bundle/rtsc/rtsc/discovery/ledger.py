"""SQLite ledger of a discovery run (DISCOVERY.md section 3).

Tables
------
``meta``            key/value (config, seed, schema)
``materials``       one row per unique structure x pressure (fingerprint), with the generation
                    request, generator provenance and the round it first appeared in
``duplicates``      every later generated structure that matched an existing material (linked,
                    never re-run)
``ml_screen``       heuristic ML pre-screen results (S2: ordering / DEFERRED only)
``dft_jobs``        one row per material sent to the DFT stage (state PENDING -> SUBMITTED -> DONE)
``verdicts``        APPEND-ONLY history: status, bracket, rule ids, assumptions, flags,
                    certificate path, goal, T*, material pressure, DFT pressure actually used,
                    stage + origin, input digest, pipeline config hash, dft_hash, certify_hash
``current_status``  a VIEW: for each (material, goal) the verdict with the largest id
``rounds``          requests, phase (proposed -> generated -> prescreened -> selected -> executed
                    -> done), counts
``policy_state``    per-round snapshot of the bandit statistics
``sessions``        one row per ``run`` / ``continue`` invocation (fresh stop limits; the config
                    used by session N is in meta ``config_session_N``)

Invariants (enforced IN SQL, so a raw connection cannot bypass them either):

* append-only: triggers reject UPDATE, DELETE and an INSERT whose id already exists
  (``recursive_triggers`` makes INSERT OR REPLACE's implicit delete hit the DELETE trigger too);
* guard (ARCH S1/S2/S3), trigger ``verdicts_guard``: a row is accepted only if
  - stage = pipeline, or status = DEFERRED, or (status = ERROR and stage in {generate, executor});
  - REJECTED_PROVEN / CERTIFIED_CANDIDATE have origin = rtsc.certify.classify (hence stage pipeline);
  - EXCLUDED_ASSUMPTION has an rtsc.stability* origin (hence stage pipeline);
  - CERTIFIED_CANDIDATE has in_scope = 1.
  ``append_verdict`` checks the same matrix in Python first (clearer errors);
* current status = the latest verdict per (material, goal) (a view, never stored separately).

Every write is one SQLite transaction (ACID; WAL journal).  Whether a CERTIFIED row COUNTS for a
goal is decided by ``rtsc.discovery.goals.counts_for`` (T*, pressures, config hashes).
"""
from __future__ import annotations

import contextlib
import datetime as _dt
import json
import sqlite3
from pathlib import Path
from typing import Any, Iterator

from ..pipeline import CLASSIFY_ORIGIN
from ..status import PROVEN_STATUSES, Status

SCHEMA_VERSION = "rtsc.discovery.ledger/2"

HEURISTIC_STAGES = {"generate": {Status.DEFERRED, Status.ERROR},
                    "scope": {Status.DEFERRED},
                    "config": {Status.DEFERRED},
                    "prescreen": {Status.DEFERRED},
                    "select": {Status.DEFERRED},
                    "policy": {Status.DEFERRED},
                    "executor": {Status.ERROR, Status.DEFERRED}}
ALL_STAGES = set(HEURISTIC_STAGES) | {"pipeline"}
ROUND_PHASES = ["proposed", "generated", "prescreened", "selected", "executed", "done"]


class LedgerSoundnessError(AssertionError):
    """A heuristic stage tried to write a proven or excluded status (S1/S2/S3)."""


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def _j(x: Any) -> str:
    return json.dumps(x, sort_keys=True, default=str)


_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS materials (
  material_id TEXT PRIMARY KEY, seq INTEGER NOT NULL UNIQUE,
  fingerprint TEXT NOT NULL, fingerprint_relaxed TEXT, coarse_key TEXT NOT NULL,
  formula TEXT NOT NULL, chemsys TEXT NOT NULL, pressure_GPa REAL NOT NULL,
  structure TEXT NOT NULL, generator TEXT NOT NULL, request_id TEXT NOT NULL,
  request TEXT NOT NULL, arm TEXT NOT NULL, round INTEGER NOT NULL, provenance TEXT NOT NULL,
  dup_of TEXT, valid INTEGER NOT NULL DEFAULT 1, created_utc TEXT NOT NULL);
CREATE INDEX IF NOT EXISTS materials_coarse ON materials(coarse_key);
CREATE INDEX IF NOT EXISTS materials_fp ON materials(fingerprint);
CREATE INDEX IF NOT EXISTS materials_fpr ON materials(fingerprint_relaxed);
CREATE TABLE IF NOT EXISTS duplicates (
  id INTEGER PRIMARY KEY AUTOINCREMENT, material_id TEXT NOT NULL REFERENCES materials(material_id),
  fingerprint TEXT NOT NULL, generator TEXT NOT NULL, request_id TEXT NOT NULL, arm TEXT NOT NULL,
  round INTEGER NOT NULL, how TEXT NOT NULL, provenance TEXT NOT NULL, created_utc TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS ml_screen (
  material_id TEXT PRIMARY KEY REFERENCES materials(material_id), round INTEGER NOT NULL,
  prescreener TEXT NOT NULL, ok INTEGER NOT NULL, result TEXT NOT NULL, ml_deferred INTEGER NOT NULL,
  created_utc TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS dft_jobs (
  id INTEGER PRIMARY KEY AUTOINCREMENT, material_id TEXT NOT NULL REFERENCES materials(material_id),
  round INTEGER NOT NULL, state TEXT NOT NULL, executor TEXT NOT NULL, external_id TEXT,
  attempts INTEGER NOT NULL DEFAULT 0, note TEXT, created_utc TEXT NOT NULL, updated_utc TEXT NOT NULL,
  UNIQUE(material_id, round));
CREATE TABLE IF NOT EXISTS verdicts (
  id INTEGER PRIMARY KEY AUTOINCREMENT, material_id TEXT NOT NULL REFERENCES materials(material_id),
  goal_id TEXT NOT NULL, round INTEGER NOT NULL, stage TEXT NOT NULL, origin TEXT NOT NULL,
  status TEXT NOT NULL, T_low_K REAL, T_up_K REAL, T_star_K REAL NOT NULL, pressure_GPa REAL NOT NULL,
  rule_ids TEXT NOT NULL, assumptions TEXT NOT NULL, flags TEXT NOT NULL, reasons TEXT NOT NULL,
  certificate_path TEXT, record_path TEXT, input_digest TEXT, config_hash TEXT, dft_hash TEXT,
  certify_hash TEXT, dft_pressure_GPa REAL,
  retryable INTEGER NOT NULL DEFAULT 0, in_scope INTEGER NOT NULL DEFAULT 1, created_utc TEXT NOT NULL);
CREATE INDEX IF NOT EXISTS verdicts_mg ON verdicts(material_id, goal_id);
CREATE TRIGGER IF NOT EXISTS verdicts_no_update BEFORE UPDATE ON verdicts
  BEGIN SELECT RAISE(ABORT, 'verdicts are append-only'); END;
CREATE TRIGGER IF NOT EXISTS verdicts_no_delete BEFORE DELETE ON verdicts
  BEGIN SELECT RAISE(ABORT, 'verdicts are append-only'); END;
CREATE TRIGGER IF NOT EXISTS verdicts_no_overwrite BEFORE INSERT ON verdicts
  WHEN NEW.id IS NOT NULL AND EXISTS (SELECT 1 FROM verdicts WHERE id = NEW.id)
  BEGIN SELECT RAISE(ABORT, 'verdicts are append-only: id exists'); END;
CREATE TRIGGER IF NOT EXISTS verdicts_guard BEFORE INSERT ON verdicts
  WHEN NOT (
    NEW.stage IN ('generate','scope','prescreen','select','policy','config','executor','pipeline')
    AND NEW.status IN ('REJECTED_PROVEN','CERTIFIED_CANDIDATE','UNDECIDED','EXCLUDED_ASSUMPTION',
                       'DEFERRED','ERROR')
    AND (NEW.stage = 'pipeline' OR NEW.status = 'DEFERRED'
         OR (NEW.status = 'ERROR' AND NEW.stage IN ('generate','executor')))
    AND (NEW.status NOT IN ('REJECTED_PROVEN','CERTIFIED_CANDIDATE')
         OR NEW.origin = 'rtsc.certify.classify')
    AND (NEW.status <> 'EXCLUDED_ASSUMPTION' OR substr(NEW.origin, 1, 14) = 'rtsc.stability')
    AND (NEW.status <> 'CERTIFIED_CANDIDATE' OR NEW.in_scope = 1))
  BEGIN SELECT RAISE(ABORT, 'ledger guard (S1/S2/S3): stage/status/origin not allowed'); END;
CREATE VIEW IF NOT EXISTS current_status AS
  SELECT material_id, goal_id, MAX(id) AS verdict_id FROM verdicts GROUP BY material_id, goal_id;
CREATE TABLE IF NOT EXISTS rounds (
  round INTEGER PRIMARY KEY, session INTEGER NOT NULL, phase TEXT NOT NULL, requests TEXT NOT NULL,
  counts TEXT NOT NULL, notes TEXT NOT NULL, started_utc TEXT NOT NULL, finished_utc TEXT);
CREATE TABLE IF NOT EXISTS policy_state (
  round INTEGER PRIMARY KEY, state TEXT NOT NULL, created_utc TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS sessions (
  id INTEGER PRIMARY KEY AUTOINCREMENT, command TEXT NOT NULL, limits TEXT NOT NULL,
  started_utc TEXT NOT NULL, finished_utc TEXT, stop_reason TEXT);
"""


class Ledger:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.con = sqlite3.connect(str(self.path), isolation_level=None, timeout=60.0)
        self.con.row_factory = sqlite3.Row
        self.con.execute("PRAGMA journal_mode=WAL")
        self.con.execute("PRAGMA foreign_keys=ON")
        self.con.execute("PRAGMA synchronous=FULL")
        self.con.execute("PRAGMA recursive_triggers=ON")   # REPLACE's implicit delete fires triggers
        self._depth = 0
        self.con.executescript(_SCHEMA)
        if self.get_meta("schema") is None:
            self.set_meta("schema", SCHEMA_VERSION)

    def close(self) -> None:
        self.con.close()

    # ------------------------------------------------------------------ transactions
    @contextlib.contextmanager
    def tx(self) -> Iterator[sqlite3.Connection]:
        """Nested-safe transaction: only the outermost level commits / rolls back."""
        if self._depth == 0:
            self.con.execute("BEGIN IMMEDIATE")
        self._depth += 1
        try:
            yield self.con
        except BaseException:
            self._depth -= 1
            if self._depth == 0:
                self.con.execute("ROLLBACK")
            raise
        else:
            self._depth -= 1
            if self._depth == 0:
                self.con.execute("COMMIT")

    def q(self, sql: str, args: tuple = ()) -> list[sqlite3.Row]:
        return list(self.con.execute(sql, args))

    # ------------------------------------------------------------------ meta / sessions
    def get_meta(self, key: str, default: Any = None) -> Any:
        r = self.q("SELECT value FROM meta WHERE key=?", (key,))
        return json.loads(r[0]["value"]) if r else default

    def set_meta(self, key: str, value: Any) -> None:
        with self.tx():
            self.con.execute("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)", (key, _j(value)))

    def start_session(self, command: str, limits: dict[str, Any]) -> int:
        with self.tx():
            cur = self.con.execute("INSERT INTO sessions(command, limits, started_utc) VALUES (?,?,?)",
                                   (command, _j(limits), _now()))
            return int(cur.lastrowid)

    def finish_session(self, sid: int, reason: str) -> None:
        with self.tx():
            self.con.execute("UPDATE sessions SET finished_utc=?, stop_reason=? WHERE id=?",
                             (_now(), reason, sid))

    # ------------------------------------------------------------------ materials
    def next_seq(self) -> int:
        r = self.q("SELECT COALESCE(MAX(seq), 0) AS m FROM materials")
        return int(r[0]["m"]) + 1

    def add_material(self, *, material_id: str, seq: int, fingerprint: str, coarse_key: str,
                     formula: str, chemsys: tuple[str, ...], pressure_GPa: float,
                     structure: dict, generator: str, request_id: str, request: dict, arm: str,
                     round_: int, provenance: dict, valid: bool = True) -> None:
        with self.tx():
            self.con.execute(
                "INSERT INTO materials(material_id, seq, fingerprint, coarse_key, formula, chemsys, "
                "pressure_GPa, structure, generator, request_id, request, arm, round, provenance, "
                "valid, created_utc) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (material_id, seq, fingerprint, coarse_key, formula, "-".join(chemsys),
                 float(pressure_GPa), _j(structure), generator, request_id, _j(request), arm,
                 round_, _j(provenance), int(valid), _now()))

    def add_duplicate(self, *, material_id: str, fingerprint: str, generator: str,
                      request_id: str, arm: str, round_: int, how: str, provenance: dict) -> None:
        with self.tx():
            self.con.execute(
                "INSERT INTO duplicates(material_id, fingerprint, generator, request_id, arm, round, "
                "how, provenance, created_utc) VALUES (?,?,?,?,?,?,?,?,?)",
                (material_id, fingerprint, generator, request_id, arm, round_, how, _j(provenance),
                 _now()))

    def materials_by_fingerprint(self, fp: str) -> list[sqlite3.Row]:
        """Fingerprint hits are only CANDIDATE matches (confirmed by the caller)."""
        return self.q("SELECT * FROM materials WHERE fingerprint=? OR fingerprint_relaxed=? "
                      "ORDER BY seq", (fp, fp))

    def material_by_fingerprint(self, fp: str) -> sqlite3.Row | None:
        r = self.q("SELECT * FROM materials WHERE fingerprint=? OR fingerprint_relaxed=? "
                   "ORDER BY seq LIMIT 1", (fp, fp))
        return r[0] if r else None

    def materials_with_coarse_key(self, key: str) -> list[sqlite3.Row]:
        return self.q("SELECT * FROM materials WHERE coarse_key=? AND valid=1 ORDER BY seq", (key,))

    def material(self, mid: str) -> sqlite3.Row:
        r = self.q("SELECT * FROM materials WHERE material_id=?", (mid,))
        if not r:
            raise KeyError(mid)
        return r[0]

    def materials(self, round_: int | None = None) -> list[sqlite3.Row]:
        if round_ is None:
            return self.q("SELECT * FROM materials ORDER BY seq")
        return self.q("SELECT * FROM materials WHERE round=? ORDER BY seq", (round_,))

    def set_relaxed_fingerprint(self, mid: str, fp: str | None, dup_of: str | None) -> None:
        with self.tx():
            self.con.execute("UPDATE materials SET fingerprint_relaxed=?, dup_of=? WHERE material_id=?",
                             (fp, dup_of, mid))

    # ------------------------------------------------------------------ ML pre-screen
    def add_ml_screen(self, mid: str, round_: int, prescreener: str, ok: bool, result: dict,
                      ml_deferred: bool) -> None:
        with self.tx():
            self.con.execute("INSERT OR REPLACE INTO ml_screen(material_id, round, prescreener, ok, "
                             "result, ml_deferred, created_utc) VALUES (?,?,?,?,?,?,?)",
                             (mid, round_, prescreener, int(ok), _j(result), int(ml_deferred), _now()))

    def ml_screen(self, mid: str) -> dict | None:
        r = self.q("SELECT * FROM ml_screen WHERE material_id=?", (mid,))
        if not r:
            return None
        d = dict(r[0])
        d["result"] = json.loads(d["result"])
        return d

    # ------------------------------------------------------------------ DFT jobs
    def add_dft_job(self, mid: str, round_: int, executor: str) -> None:
        with self.tx():
            self.con.execute("INSERT OR IGNORE INTO dft_jobs(material_id, round, state, executor, "
                             "created_utc, updated_utc) VALUES (?,?,?,?,?,?)",
                             (mid, round_, "PENDING", executor, _now(), _now()))

    def set_job_state(self, mid: str, round_: int, state: str, *, external_id: str | None = None,
                      note: str | None = None, attempt: bool = False) -> None:
        with self.tx():
            self.con.execute(
                "UPDATE dft_jobs SET state=?, external_id=COALESCE(?, external_id), "
                "note=COALESCE(?, note), attempts=attempts+?, updated_utc=? WHERE material_id=? AND round=?",
                (state, external_id, note, int(attempt), _now(), mid, round_))

    def jobs(self, round_: int | None = None, states: tuple[str, ...] | None = None) -> list[sqlite3.Row]:
        sql, args = "SELECT * FROM dft_jobs WHERE 1=1", []
        if round_ is not None:
            sql += " AND round=?"
            args.append(round_)
        if states:
            sql += f" AND state IN ({','.join('?' * len(states))})"
            args.extend(states)
        return self.q(sql + " ORDER BY id", tuple(args))

    # ------------------------------------------------------------------ verdicts
    def append_verdict(self, *, material_id: str, goal_id: str, round_: int, stage: str, origin: str,
                       status: Status, T_star_K: float, pressure_GPa: float,
                       T_low_K: float | None = None, T_up_K: float | None = None,
                       rule_ids=(), assumptions=(), flags=(), reasons=(),
                       certificate_path: str | None = None, record_path: str | None = None,
                       input_digest: str | None = None, config_hash: str | None = None,
                       dft_hash: str | None = None, certify_hash: str | None = None,
                       dft_pressure_GPa: float | None = None, retryable: bool = False,
                       in_scope: bool = True) -> int:
        status = Status(status)
        if stage not in ALL_STAGES:
            raise ValueError(f"unknown ledger stage {stage!r}")
        if stage != "pipeline" and status not in HEURISTIC_STAGES[stage]:
            raise LedgerSoundnessError(f"S2 violated: stage {stage!r} may not produce {status.value}")
        if status in PROVEN_STATUSES and (stage != "pipeline" or origin != CLASSIFY_ORIGIN):
            raise LedgerSoundnessError(f"S1 violated: {status.value} from {stage}/{origin}")
        if status is Status.EXCLUDED_ASSUMPTION and (stage != "pipeline"
                                                     or not origin.startswith("rtsc.stability")):
            raise LedgerSoundnessError(f"S3 violated: EXCLUDED_ASSUMPTION from {stage}/{origin}")
        if status is Status.CERTIFIED_CANDIDATE and not in_scope:
            raise LedgerSoundnessError("CERTIFIED for a goal whose pressure range excludes it")
        with self.tx():
            cur = self.con.execute(
                "INSERT INTO verdicts(material_id, goal_id, round, stage, origin, status, T_low_K, "
                "T_up_K, T_star_K, pressure_GPa, rule_ids, assumptions, flags, reasons, "
                "certificate_path, record_path, input_digest, config_hash, dft_hash, certify_hash, "
                "dft_pressure_GPa, retryable, in_scope, created_utc) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (material_id, goal_id, round_, stage, origin, status.value, T_low_K, T_up_K,
                 float(T_star_K), float(pressure_GPa), _j(list(rule_ids)), _j(list(assumptions)),
                 _j(list(flags)), _j(list(reasons)), certificate_path, record_path, input_digest,
                 config_hash, dft_hash, certify_hash, dft_pressure_GPa, int(retryable),
                 int(in_scope), _now()))
            return int(cur.lastrowid)

    def current(self, goal_id: str | None = None, material_id: str | None = None) -> list[dict]:
        sql = ("SELECT v.*, m.formula, m.chemsys, m.arm, m.generator, m.round AS gen_round, "
               "m.fingerprint, m.dup_of, m.seq FROM current_status c "
               "JOIN verdicts v ON v.id=c.verdict_id JOIN materials m ON m.material_id=c.material_id "
               "WHERE 1=1")
        args: list[Any] = []
        if goal_id is not None:
            sql += " AND c.goal_id=?"
            args.append(goal_id)
        if material_id is not None:
            sql += " AND c.material_id=?"
            args.append(material_id)
        return [_vrow(r) for r in self.q(sql + " ORDER BY m.seq, c.goal_id", tuple(args))]

    def current_one(self, material_id: str, goal_id: str) -> dict | None:
        r = self.current(goal_id, material_id)
        return r[0] if r else None

    def history(self, material_id: str | None = None) -> list[dict]:
        if material_id is None:
            return [_vrow(r) for r in self.q("SELECT * FROM verdicts ORDER BY id")]
        return [_vrow(r) for r in self.q("SELECT * FROM verdicts WHERE material_id=? ORDER BY id",
                                         (material_id,))]

    # ------------------------------------------------------------------ rounds / policy
    def last_round(self) -> sqlite3.Row | None:
        r = self.q("SELECT * FROM rounds ORDER BY round DESC LIMIT 1")
        return r[0] if r else None

    def round_row(self, k: int) -> sqlite3.Row | None:
        r = self.q("SELECT * FROM rounds WHERE round=?", (k,))
        return r[0] if r else None

    def rounds(self) -> list[dict]:
        out = []
        for r in self.q("SELECT * FROM rounds ORDER BY round"):
            d = dict(r)
            for k in ("requests", "counts", "notes"):
                d[k] = json.loads(d[k])
            out.append(d)
        return out

    def start_round(self, k: int, session: int, requests: list[dict]) -> None:
        with self.tx():
            self.con.execute("INSERT INTO rounds(round, session, phase, requests, counts, notes, "
                             "started_utc) VALUES (?,?,?,?,?,?,?)",
                             (k, session, "proposed", _j(requests), _j({}), _j([]), _now()))

    def set_phase(self, k: int, phase: str, counts: dict | None = None, note: str | None = None) -> None:
        assert phase in ROUND_PHASES
        with self.tx():
            r = self.round_row(k)
            c = json.loads(r["counts"])
            c.update(counts or {})
            notes = json.loads(r["notes"])
            if note:
                notes.append(note)
            self.con.execute("UPDATE rounds SET phase=?, counts=?, notes=?, finished_utc=? WHERE round=?",
                             (phase, _j(c), _j(notes), _now() if phase == "done" else None, k))

    def add_round_note(self, k: int, note: str) -> None:
        r = self.round_row(k)
        self.set_phase(k, r["phase"], note=note)

    def save_policy_state(self, k: int, state: dict) -> None:
        with self.tx():
            self.con.execute("INSERT OR REPLACE INTO policy_state(round, state, created_utc) VALUES (?,?,?)",
                             (k, _j(state), _now()))

    # ------------------------------------------------------------------ accounting
    def accounting(self, goal_ids: list[str]) -> dict[str, Any]:
        """S4 for the ledger: every material has exactly one current status per goal, each
        pointing at its latest verdict; every duplicate links to an existing material; per round,
        the recorded number of generator outputs = material rows + duplicate rows of that round.
        Raises AssertionError."""
        mats = [r["material_id"] for r in self.q("SELECT material_id FROM materials")]
        cur = self.q("SELECT material_id, goal_id, verdict_id FROM current_status")
        pairs: dict[tuple[str, str], int] = {}
        for r in cur:
            k = (r["material_id"], r["goal_id"])
            if k in pairs:
                raise AssertionError(f"two current statuses for {k}")
            pairs[k] = r["verdict_id"]
        problems = []
        for m in mats:
            for g in goal_ids:
                if (m, g) not in pairs:
                    problems.append(f"{m}/{g}: no current status")
        latest = {(r["material_id"], r["goal_id"]): r["m"] for r in
                  self.q("SELECT material_id, goal_id, MAX(id) AS m FROM verdicts GROUP BY material_id, goal_id")}
        for k, vid in pairs.items():
            if latest.get(k) != vid:
                problems.append(f"{k}: current status is not the latest verdict")
        orphans = self.q("SELECT d.id FROM duplicates d LEFT JOIN materials m ON m.material_id=d.material_id "
                         "WHERE m.material_id IS NULL")
        if orphans:
            problems.append(f"{len(orphans)} duplicates without material")
        # every generator output of a round is exactly one material row or one duplicate row
        # created in that round (counts.generated is written in the same transaction)
        generated = 0
        for r in self.q("SELECT round, phase, counts FROM rounds ORDER BY round"):
            if r["phase"] == "proposed":
                continue
            c = json.loads(r["counts"])
            if "generated" not in c:
                problems.append(f"round {r['round']}: generator output count missing")
                continue
            n_m = self.q("SELECT COUNT(*) AS n FROM materials WHERE round=?", (r["round"],))[0]["n"]
            n_d = self.q("SELECT COUNT(*) AS n FROM duplicates WHERE round=?", (r["round"],))[0]["n"]
            if int(c["generated"]) != n_m + n_d:
                problems.append(f"round {r['round']}: {c['generated']} generator outputs but "
                                f"{n_m} materials + {n_d} duplicates")
            generated += int(c["generated"])
        if problems:
            raise AssertionError("ledger accounting failed: " + "; ".join(problems[:20]))
        n_dup = self.q("SELECT COUNT(*) AS n FROM duplicates")[0]["n"]
        return {"materials": len(mats), "duplicates": int(n_dup),
                "generated": generated, "goals": list(goal_ids),
                "current_pairs": len(pairs)}


def _vrow(r: sqlite3.Row) -> dict:
    d = dict(r)
    for k in ("rule_ids", "assumptions", "flags", "reasons"):
        if k in d and isinstance(d[k], str):
            d[k] = json.loads(d[k])
    return d
