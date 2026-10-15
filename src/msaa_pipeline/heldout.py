from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

from .decode import CONDITIONS, decode_subjects, get_subjects


def _load_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _split_outer_subject_indices(n_subjects: int, rng: np.random.Generator):
    idx = np.arange(n_subjects)
    rng.shuffle(idx)
    midpoint = n_subjects // 2
    return idx[:midpoint], idx[midpoint:]


def _subset_subjects(subjects, indices):
    return [subjects[int(i)] for i in indices]


def _n_archetypes(subjects) -> int:
    return int(np.asarray(subjects[0]["S"]).shape[0])


def rank_archetypes_on_subjects(
    ranking_subjects,
    analysis_type: str,
    *,
    nfolds: int = 2,
    nreps: int = 20,
    seed: int = 42,
):
    """Rank single archetypes using ranking participants only.

    This mirrors the final held-out cross-fitting notebook: accuracy is the
    primary ranking criterion, then mean rank, then archetype index for a
    deterministic tie break.
    """
    K = _n_archetypes(ranking_subjects)
    rows = []
    for archetype in range(K):
        dec = decode_subjects(
            ranking_subjects,
            analysis_type,
            [archetype],
            nfolds=nfolds,
            nreps=nreps,
            seed=seed,
        )
        rows.append({
            "archetype": int(archetype),
            "mean_accuracy": float(dec["accuracy"].mean()),
            "sem_accuracy": float(dec["accuracy"].std() / np.sqrt(max(len(dec), 1))),
            "rank_mean": float(dec["rank"].mean()),
            "error_mean": float(dec["error"].mean()),
        })

    per_arch = pd.DataFrame(rows)
    ranked = (
        per_arch
        .sort_values(
            ["mean_accuracy", "rank_mean", "archetype"],
            ascending=[False, False, True],
        )["archetype"]
        .astype(int)
        .tolist()
    )
    return ranked, per_arch


def evaluate_topm_on_subjects(
    evaluation_subjects,
    analysis_type: str,
    ranked_archetypes,
    *,
    top_m_values=(1, 2, 3, 5, 10),
    nfolds: int = 2,
    nreps: int = 20,
    seed: int = 42,
):
    """Evaluate a ranking exclusively on held-out evaluation participants."""
    K = _n_archetypes(evaluation_subjects)

    full = decode_subjects(
        evaluation_subjects,
        analysis_type,
        None,
        nfolds=nfolds,
        nreps=nreps,
        seed=seed,
    )
    full_mean = float(full["accuracy"].mean())
    full_err = float(full["accuracy"].std() / np.sqrt(max(len(full), 1)))

    # All empirical Pieman reconstructions used here have 300 timepoints after
    # the canonical rest trim. Read T from the reconstructed decoder output
    # indirectly via the first subject's sXC/S orientation would be brittle;
    # 1/300 is the historical analysis definition and is explicit here.
    chance = 1.0 / 300.0

    rows = []
    for requested_m in top_m_values:
        m = min(int(requested_m), K)
        selected = ranked_archetypes[:m]
        dec = decode_subjects(
            evaluation_subjects,
            analysis_type,
            selected,
            nfolds=nfolds,
            nreps=nreps,
            seed=seed,
        )
        mean = float(dec["accuracy"].mean())
        denom = full_mean - chance
        retention = (mean - chance) / denom if abs(denom) > 1e-12 else np.nan
        rows.append({
            "top_m": m,
            "selected_archetypes": str(selected),
            "mean": mean,
            "err": float(dec["accuracy"].std() / np.sqrt(max(len(dec), 1))),
            "std": float(dec["accuracy"].std()),
            "count": int(len(dec)),
            "rank_mean": float(dec["rank"].mean()),
            "error_mean": float(dec["error"].mean()),
            "full_mean_same_eval": full_mean,
            "full_err_same_eval": full_err,
            "chance": chance,
            "retention_above_chance": float(retention) if np.isfinite(retention) else np.nan,
        })
    return pd.DataFrame(rows)


def _completed_outer_reps(df: pd.DataFrame, analysis_type: str, condition: str, K: int):
    if df.empty:
        return set()
    required = {"analysis_type", "fit_scope", "condition", "K", "outer_rep", "direction"}
    if not required.issubset(df.columns):
        return set()
    sub = df[
        (df["analysis_type"].astype(str) == str(analysis_type))
        & (df["fit_scope"].astype(str) == "across")
        & (df["condition"].astype(str) == str(condition))
        & (df["K"].astype(int) == int(K))
    ]
    complete = set()
    for outer_rep, grp in sub.groupby("outer_rep"):
        directions = set(grp["direction"].dropna().astype(str))
        if {"A_to_B", "B_to_A"}.issubset(directions):
            complete.add(int(outer_rep))
    return complete


def _dedupe_and_save(run_df: pd.DataFrame, rank_df: pd.DataFrame, run_path: Path, rank_path: Path):
    if not run_df.empty:
        run_df = run_df.drop_duplicates(
            subset=["analysis_type", "fit_scope", "condition", "K", "outer_rep", "direction", "top_m"],
            keep="last",
        ).sort_values(["condition", "K", "outer_rep", "direction", "top_m"]).reset_index(drop=True)
        _write_csv(run_df, run_path)
    if not rank_df.empty:
        rank_df = rank_df.drop_duplicates(
            subset=["analysis_type", "fit_scope", "condition", "K", "outer_rep", "direction", "archetype"],
            keep="last",
        ).sort_values(["condition", "K", "outer_rep", "direction", "archetype"]).reset_index(drop=True)
        _write_csv(rank_df, rank_path)
    return run_df, rank_df


def _sem(x):
    x = pd.Series(x).dropna()
    return np.nan if len(x) <= 1 else float(x.std(ddof=1) / np.sqrt(len(x)))


def write_heldout_summaries(run_df: pd.DataFrame, rank_df: pd.DataFrame, outdir: Path, n_outer_reps: int):
    """Write corrected outer-split summaries and selection stability tables.

    A->B and B->A are two directions of the same outer participant split, so
    they are averaged within outer_rep before SEM is computed across outer
    repetitions. This is the uncertainty calculation used by the current
    manuscript Figure 4.
    """
    if run_df.empty:
        return

    outer = (
        run_df
        .groupby(["analysis_type", "fit_scope", "condition", "K", "top_m", "outer_rep"], as_index=False)
        .agg(
            mean=("mean", "mean"),
            full_mean_same_eval=("full_mean_same_eval", "mean"),
            retention_above_chance=("retention_above_chance", "mean"),
        )
    )
    summary = (
        outer
        .groupby(["analysis_type", "fit_scope", "condition", "K", "top_m"], as_index=False)
        .agg(
            mean=("mean", "mean"),
            err=("mean", _sem),
            std=("mean", "std"),
            count=("outer_rep", "size"),
            mean_full_same_eval=("full_mean_same_eval", "mean"),
            mean_retention=("retention_above_chance", "mean"),
            err_retention=("retention_above_chance", _sem),
        )
    )
    _write_csv(summary, outdir / "heldout_topm_summary.csv")

    if rank_df.empty:
        return
    stability_rows = []
    group_cols = ["analysis_type", "fit_scope", "condition", "K", "outer_rep", "direction"]
    top_ms = sorted(set(run_df["top_m"].astype(int)))
    for keys, sub in rank_df.groupby(group_cols):
        ranked = (
            sub.sort_values(
                ["mean_accuracy", "rank_mean", "archetype"],
                ascending=[False, False, True],
            )["archetype"].astype(int).tolist()
        )
        metadata = dict(zip(group_cols, keys))
        for top_m in top_ms:
            m = min(int(top_m), len(ranked))
            for rank_position, archetype in enumerate(ranked[:m], start=1):
                stability_rows.append({
                    **metadata,
                    "top_m": m,
                    "rank_position": int(rank_position),
                    "archetype": int(archetype),
                })
    stability_long = pd.DataFrame(stability_rows)
    stability = (
        stability_long
        .groupby(["analysis_type", "fit_scope", "condition", "K", "top_m", "archetype"], as_index=False)
        .agg(
            selected_count=("archetype", "size"),
            mean_rank_position=("rank_position", "mean"),
        )
    )
    # Each completed outer split contributes two independent directions.
    denom = max(int(n_outer_reps) * 2, 1)
    stability["selection_frequency"] = stability["selected_count"] / denom
    stability = stability.sort_values(
        ["condition", "K", "top_m", "selection_frequency", "archetype"],
        ascending=[True, True, True, False, True],
    ).reset_index(drop=True)
    _write_csv(stability, outdir / "heldout_archetype_selection_stability.csv")


def run_heldout_topm_job(
    output_root: Path,
    analysis_type: str,
    K_values,
    *,
    top_m=(1, 2, 3, 5, 10),
    n_outer_reps: int = 20,
    nfolds: int = 2,
    nreps: int = 20,
    seed: int = 42,
):
    """Run/resume the canonical held-out top-m cross-fit analysis.

    Only across-condition fits are supported. For each condition/K/outer split,
    one half ranks individual archetypes and the other half evaluates top-m;
    the direction is then reversed. No same-data top-m result is produced.
    """
    outdir = output_root / f"msaa_heldout_topm_outputs_{analysis_type}_across"
    outdir.mkdir(parents=True, exist_ok=True)
    run_path = outdir / "heldout_topm_crossfit_runs.csv"
    rank_path = outdir / "heldout_per_archetype_ranking_runs.csv"

    run_df = _load_csv(run_path)
    rank_df = _load_csv(rank_path)

    for K in map(int, K_values):
        for condition_index, condition in enumerate(CONDITIONS):
            run_done = _completed_outer_reps(run_df, analysis_type, condition, K)
            rank_done = _completed_outer_reps(rank_df, analysis_type, condition, K)
            completed = run_done.intersection(rank_done)
            if len(completed) >= int(n_outer_reps):
                print(f"Held-out {analysis_type}/across/{condition}/K={K}: already complete; skipping")
                continue

            subjects, source = get_subjects(output_root, analysis_type, "across", K, condition)
            if len(subjects) < 4:
                raise ValueError(
                    f"Held-out top-m needs at least 4 subjects; got {len(subjects)} for {condition}/K={K}"
                )
            print(f"Held-out {analysis_type}/across/{condition}/K={K} | N={len(subjects)}")

            for outer_rep in range(int(n_outer_reps)):
                if outer_rep in completed:
                    continue

                seed_sequence = np.random.SeedSequence([
                    int(seed), int(K), int(condition_index), int(outer_rep)
                ])
                outer_seed = int(seed_sequence.generate_state(1, dtype=np.uint32)[0])
                outer_rng = np.random.default_rng(outer_seed)
                A_idx, B_idx = _split_outer_subject_indices(len(subjects), outer_rng)

                new_runs = []
                new_ranks = []
                for direction, ranking_idx, evaluation_idx in [
                    ("A_to_B", A_idx, B_idx),
                    ("B_to_A", B_idx, A_idx),
                ]:
                    ranking_subjects = _subset_subjects(subjects, ranking_idx)
                    evaluation_subjects = _subset_subjects(subjects, evaluation_idx)
                    direction_offset = 0 if direction == "A_to_B" else 1
                    ranking_seed = seed + 100000 * outer_rep + 1000 * K + direction_offset
                    evaluation_seed = seed + 200000 * outer_rep + 1000 * K + direction_offset

                    ranked, per_arch = rank_archetypes_on_subjects(
                        ranking_subjects,
                        analysis_type,
                        nfolds=nfolds,
                        nreps=nreps,
                        seed=ranking_seed,
                    )
                    per_arch = per_arch.copy()
                    per_arch["analysis_type"] = analysis_type
                    per_arch["fit_scope"] = "across"
                    per_arch["condition"] = condition
                    per_arch["K"] = K
                    per_arch["outer_rep"] = outer_rep
                    per_arch["direction"] = direction
                    per_arch["outer_seed"] = outer_seed
                    per_arch["n_ranking_subjects"] = len(ranking_subjects)
                    new_ranks.append(per_arch)

                    heldout = evaluate_topm_on_subjects(
                        evaluation_subjects,
                        analysis_type,
                        ranked,
                        top_m_values=top_m,
                        nfolds=nfolds,
                        nreps=nreps,
                        seed=evaluation_seed,
                    )
                    heldout["analysis_type"] = analysis_type
                    heldout["fit_scope"] = "across"
                    heldout["condition"] = condition
                    heldout["K"] = K
                    heldout["outer_rep"] = outer_rep
                    heldout["direction"] = direction
                    heldout["outer_seed"] = outer_seed
                    heldout["n_ranking_subjects"] = len(ranking_subjects)
                    heldout["n_evaluation_subjects"] = len(evaluation_subjects)
                    heldout["fit_source_file"] = str(source)
                    heldout["fit_source_basename"] = source.name
                    new_runs.append(heldout)

                run_df = pd.concat([run_df, *new_runs], ignore_index=True) if not run_df.empty else pd.concat(new_runs, ignore_index=True)
                rank_df = pd.concat([rank_df, *new_ranks], ignore_index=True) if not rank_df.empty else pd.concat(new_ranks, ignore_index=True)
                run_df, rank_df = _dedupe_and_save(run_df, rank_df, run_path, rank_path)

                if (outer_rep + 1) % 5 == 0 or (outer_rep + 1) == int(n_outer_reps):
                    print(f"  completed outer rep {outer_rep + 1}/{n_outer_reps}")

            write_heldout_summaries(run_df, rank_df, outdir, n_outer_reps)

    write_heldout_summaries(run_df, rank_df, outdir, n_outer_reps)
    return run_df, rank_df
