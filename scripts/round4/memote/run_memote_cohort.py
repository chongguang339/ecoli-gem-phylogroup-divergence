#!/usr/bin/env python3
"""Resume-safe MEMOTE snapshot run for the frozen final-model cohort."""
import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def models_in(root):
    models = sorted(p for p in list(root.glob("*.xml")) + list(root.glob("*/*.xml"))
                    if not p.name.endswith("-draft.xml"))
    ids = [p.stem for p in models]
    if len(ids) != len(set(ids)):
        raise SystemExit("Duplicate final-model basenames; resolve before running.")
    return models


def score_from_html(path):
    data = path.read_bytes()
    match = re.search(rb'"score":(\{"sections":.*?\}),"cards"', data)
    if not match:
        raise ValueError("MEMOTE score object absent from HTML")
    item = json.loads(match.group(1))
    result = {"total_score": item.get("total_score")}
    for section in item.get("sections", []):
        result["section_" + str(section.get("section"))] = section.get("score")
    return result


def run_one(model, output, executable, timeout_seconds):
    tag = model.stem
    html = output / "reports" / (tag + ".html")
    log = output / "logs" / (tag + ".log")
    status = output / "status" / (tag + ".json")
    started = time.time()
    if status.exists() and html.exists() and html.stat().st_size > 1000:
        try:
            previous = json.loads(status.read_text(encoding="utf-8"))
            if previous.get("exit_code") == 0 and previous.get("model_xml") == str(model):
                score_from_html(html)
                return previous
        except (ValueError, KeyError, OSError, json.JSONDecodeError):
            pass
    command = [executable, "report", "snapshot", "--filename", str(html), str(model)]
    with log.open("wb") as stream:
        try:
            completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                       timeout=timeout_seconds, check=False)
            code = completed.returncode
            error = ""
        except subprocess.TimeoutExpired:
            code, error = 124, "timeout"
    result = {"model_id": tag, "model_xml": str(model), "report_html": str(html),
              "log_file": str(log), "exit_code": code, "error": error,
              "wall_seconds": round(time.time() - started, 1)}
    if code == 0:
        try:
            result.update(score_from_html(html))
        except (ValueError, OSError, json.JSONDecodeError) as exc:
            result["exit_code"], result["error"] = 125, "score parse: " + str(exc)
    temp = status.with_suffix(".json.tmp")
    temp.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(temp, status)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", type=Path, required=True, help="Folder containing one subfolder per final SBML model")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout-hours", type=float, default=3)
    parser.add_argument("--memote", default="memote")
    parser.add_argument("--cohort-manifest", type=Path, help="Frozen annotated-model TSV with Archive_Member and Annotated_SHA256")
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    models = models_in(args.models.resolve())
    print(f"Discovered {len(models)} final XML files in {args.models}", flush=True)
    if len(models) != args.expected:
        raise SystemExit(f"Expected {args.expected} models; check cohort path and archive before running.")
    if args.cohort_manifest:
        with args.cohort_manifest.open(newline="", encoding="utf-8-sig") as handle:
            expected_rows = list(csv.DictReader(handle, delimiter="\t"))
        expected_by_id = {row["Genome_ID"]: row for row in expected_rows}
        if len(expected_by_id) != args.expected or {p.stem for p in models} != set(expected_by_id):
            raise SystemExit("Model IDs do not match the frozen cohort manifest.")
        for path in models:
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest != expected_by_id[path.stem]["Annotated_SHA256"]:
                raise SystemExit(f"Annotated model hash mismatch: {path.name}")
        print("All model IDs and annotated SHA256 values match the frozen cohort manifest.", flush=True)
    executable = shutil.which(args.memote)
    if not executable:
        raise SystemExit("MEMOTE executable not found in active environment.")
    args.output.mkdir(parents=True, exist_ok=True)
    for name in ("reports", "logs", "status"):
        (args.output / name).mkdir(exist_ok=True)
    manifest = args.output / "model_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["model_id", "model_xml", "bytes"])
        writer.writerows((p.stem, str(p), p.stat().st_size) for p in models)
    version = subprocess.run([executable, "--version"], capture_output=True, text=True, check=False)
    env = {"memote_executable": executable, "memote_version": (version.stdout + version.stderr).strip(),
           "model_root": str(args.models.resolve()), "expected": args.expected, "workers": args.workers}
    (args.output / "environment.json").write_text(json.dumps(env, indent=2), encoding="utf-8")
    if args.preflight:
        print(f"Preflight passed. Manifest: {manifest}")
        return
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = [pool.submit(run_one, p, args.output, executable, int(args.timeout_hours * 3600)) for p in models]
        for index, future in enumerate(as_completed(jobs), 1):
            row = future.result()
            if index % 25 == 0 or row["exit_code"] != 0:
                print(f"{index}/{len(models)} {row['model_id']} exit={row['exit_code']}", flush=True)
    rows = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((args.output / "status").glob("*.json"))]
    fields = sorted({key for row in rows for key in row})
    with (args.output / "memote_all_models.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    passed = sum(row.get("exit_code") == 0 for row in rows)
    print(f"Completed: {passed}/{len(models)} successful. Summary: {args.output / 'memote_all_models.csv'}")
    if passed != len(models):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
