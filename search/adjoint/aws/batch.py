"""batch.py -- send fine adjoint jobs to AWS Batch and bring their gradients back (owner 2026-11-28: AWS Batch).

    python batch.py package <job.json> <in.tgz>                      # the job + both repos' code (git HEAD)
    python batch.py submit  <job.json> [--config aws_config.json]    # package, upload, submit; prints the job id
    python batch.py status  [<job id>] [--config ...]                # every job of the ledger, or one
    python batch.py wait    <job id> [--poll 120] [--config ...]     # until it ends
    python batch.py fetch   <job id> <dest dir> [--config ...]       # out.tgz -> dest (result.json, gradient.json, ...)

aws_config.json (aws_config.example.json; the stack's outputs, cloudformation.yaml):
  {"region": "eu-west-2", "bucket": "...", "prefix": "fine-adjoint", "job_queue": "...", "job_definition": "...",
   "vcpus": 64, "memory_mib": 480000, "np": 64, "keep_state": false, "ledger": "aws_jobs.jsonl",
   "adjoint_repo": "../stem-racing-adjoint"}
AWS credentials come from the usual places (environment, ~/.aws, an instance role); nothing is stored here.

The code each job runs is a snapshot of both repositories at their committed HEAD (`git archive`), so a result
says exactly which code made it (the commits are in result.json's `code`). Uncommitted changes are NOT sent:
package refuses a dirty tree unless --allow-dirty.
"""
from __future__ import annotations

import argparse
import io
import json
import os
import subprocess
import sys
import tarfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVO = HERE.parents[2]                                       # the evolution repo root
EVO_PATHS = ["car", "search", "cfd", "repo_paths.py", "requirements.txt"]
ADJ_PATHS = ["part2-simulation", "part4-simulation", "research/sidepod_opt/dafoam"]
ADJ_TOP = "research/sidepod_opt"                            # its top-level files only (the pipeline's scripts)
# The shared fine-adjoint core (owner 2026-11-28: "one core, two users"): the whole-car hybrid's DAFoam patch (lift,
# per-solve GMRES record, mesh-quality limits by env) and its gradient map (drag, lift, mass, race time), taken from
# that branch at a pinned commit -- no copy lives here. link.py reproduces grad.py's drag map (cosine 1.0000 on the
# p0 team-car adjoint, search/adjoint/xcheck) and stays as the per-job cross-check.
CORE_REF = "3101248d08aeabaa9f1d8a2c7aca7cea994a08d9"     # claude/wizardly-brahmagupta-1edalc, 2026-11-28 14:08
CORE_PATHS = ["carhybrid/patch_daf.py", "carhybrid/grad.py"]


def cfg_load(p):
    c = json.loads(Path(p).read_text())
    c.setdefault("prefix", "fine-adjoint"); c.setdefault("ledger", str(HERE / "aws_jobs.jsonl"))
    return c


def _git(repo, *a):
    return subprocess.run(["git", "-C", str(repo), *a], check=True, capture_output=True, text=True).stdout


def _archive(repo: Path, paths, prefix: str, tf: tarfile.TarFile, allow_dirty: bool, ref="HEAD"):
    if ref == "HEAD" and not allow_dirty and _git(repo, "status", "--porcelain", "--", *paths).strip():
        raise SystemExit(f"{repo}: uncommitted changes under {paths}; commit them (or --allow-dirty)")
    if ref != "HEAD" and subprocess.run(["git", "-C", str(repo), "cat-file", "-e", f"{ref}^{{commit}}"], stderr=subprocess.DEVNULL).returncode != 0:
        subprocess.run(["git", "-C", str(repo), "fetch", "-q", "--depth", "1", "origin", ref], check=True)
    data = subprocess.run(["git", "-C", str(repo), "archive", "--format=tar", f"--prefix={prefix}/", ref, "--", *paths],
                          check=True, capture_output=True).stdout
    with tarfile.open(fileobj=io.BytesIO(data)) as src:
        for m in src.getmembers():
            tf.addfile(m, src.extractfile(m) if m.isfile() else None)
    return _git(repo, "rev-parse", ref).strip()


def package(job_path: Path, out: Path, adjoint_repo: Path, allow_dirty=False) -> dict:
    job = json.loads(Path(job_path).read_text())
    car = job.get("car", {})
    with tarfile.open(out, "w:gz") as tf:
        code = {"evolution": _archive(EVO, EVO_PATHS, "evolution", tf, allow_dirty)}
        top = [p for p in _git(adjoint_repo, "ls-files", ADJ_TOP).split() if p.count("/") == 2]
        code["adjoint"] = _archive(adjoint_repo, ADJ_PATHS + top, "adjoint", tf, allow_dirty)
        if not job.get("no_core"):
            code["core"] = _archive(EVO, CORE_PATHS, "core", tf, allow_dirty, ref=job.get("core_ref", CORE_REF))
        if car.get("cfd_car_dir"):                          # a car built elsewhere travels with the job
            tf.add(car["cfd_car_dir"], arcname="car/cfd_car")
            car["cfd_car_dir"] = "car/cfd_car"
            if car.get("layout_car_json"):
                tf.add(car["layout_car_json"], arcname="car/layout_car.json")
                car["layout_car_json"] = "car/layout_car.json"
        job["code"] = code
        b = json.dumps(job, indent=1).encode()
        ti = tarfile.TarInfo("job.json"); ti.size = len(b); ti.mtime = int(time.time())
        tf.addfile(ti, io.BytesIO(b))
    return job


def _clients(c):
    import boto3
    s = boto3.session.Session(region_name=c.get("region"))
    return s.client("s3"), s.client("batch")


def submit(job_path, c, allow_dirty=False, s3=None, batch=None):
    job = json.loads(Path(job_path).read_text())
    jid = job.get("id") or Path(job_path).stem
    stamp = time.strftime("%Y%m%d-%H%M%S", time.gmtime())
    key = f"{c['prefix'].strip('/')}/{jid}-{stamp}"
    tgz = Path(job_path).with_suffix(".in.tgz")
    job = package(Path(job_path), tgz, Path(c.get("adjoint_repo", EVO.parent / "stem-racing-adjoint")), allow_dirty)
    if s3 is None:
        s3, batch = _clients(c)
    with open(tgz, "rb") as f:
        s3.put_object(Bucket=c["bucket"], Key=f"{key}/in.tgz", Body=f.read())
    env = [{"name": "S3_IN", "value": f"s3://{c['bucket']}/{key}/in.tgz"},
           {"name": "S3_OUT", "value": f"s3://{c['bucket']}/{key}/"},
           {"name": "NP", "value": str(c.get("np", c.get("vcpus", 64)))},
           {"name": "KEEP_STATE", "value": "1" if c.get("keep_state") else "0"}]
    r = batch.submit_job(jobName=f"fine-adjoint-{jid}"[:128].replace(".", "-"), jobQueue=c["job_queue"],
                         jobDefinition=c["job_definition"],
                         containerOverrides={"environment": env, "resourceRequirements": [
                             {"type": "VCPU", "value": str(c.get("vcpus", 64))},
                             {"type": "MEMORY", "value": str(c.get("memory_mib", 480000))}]})
    row = {"job_id": r["jobId"], "id": jid, "key": key, "bucket": c["bucket"], "submitted": stamp, "code": job["code"]}
    with open(c["ledger"], "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row))
    return row


def ledger(c):
    p = Path(c["ledger"])
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()] if p.exists() else []


def _row(c, job_id):
    for r in ledger(c):
        if job_id in (r["job_id"], r["id"]):
            return r
    raise SystemExit(f"{job_id}: not in the ledger {c['ledger']}")


def status(c, job_id=None, batch=None):
    rows = [_row(c, job_id)] if job_id else ledger(c)
    if not rows:
        return []
    if batch is None:
        _, batch = _clients(c)
    out = []
    for i in range(0, len(rows), 100):
        part = rows[i:i + 100]
        d = {j["jobId"]: j for j in batch.describe_jobs(jobs=[r["job_id"] for r in part])["jobs"]}
        for r in part:
            j = d.get(r["job_id"], {})
            out.append({"id": r["id"], "job_id": r["job_id"], "status": j.get("status"), "reason": j.get("statusReason"),
                        "started": j.get("startedAt"), "stopped": j.get("stoppedAt")})
    for o in out:
        print(json.dumps(o))
    return out


def wait(c, job_id, poll=120, batch=None):
    while True:
        s = status(c, job_id, batch=batch)[0]
        if s["status"] in ("SUCCEEDED", "FAILED"):
            return s
        time.sleep(poll)


def fetch(c, job_id, dest, s3=None):
    r = _row(c, job_id)
    if s3 is None:
        s3, _ = _clients(c)
    dest = Path(dest); dest.mkdir(parents=True, exist_ok=True)
    got = []
    for name in ("out.tgz", "result.json", "fine_adjoint.log", "state.tgz"):
        try:
            body = s3.get_object(Bucket=r["bucket"], Key=f"{r['key']}/{name}")["Body"].read()
        except Exception as e:                                 # noqa: BLE001  (a job that died early has no out.tgz)
            if name in ("out.tgz", "result.json"):
                print(f"{name}: not there ({type(e).__name__})")
            continue
        (dest / name).write_bytes(body); got.append(name)
    if "out.tgz" in got:
        with tarfile.open(dest / "out.tgz") as tf:
            tf.extractall(dest, filter="data")
    if (dest / "result.json").exists():
        res = json.loads((dest / "result.json").read_text())
        print(json.dumps({k: res.get(k) for k in ("id", "usable", "error", "production", "steady", "wall_s")}, default=str))
    return got


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["package", "submit", "status", "wait", "fetch"])
    ap.add_argument("args", nargs="*")
    ap.add_argument("--config", default=str(HERE / "aws_config.json"))
    ap.add_argument("--poll", type=int, default=120)
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--adjoint-repo", default=None)
    a = ap.parse_args(argv)
    if a.cmd == "package":
        package(Path(a.args[0]), Path(a.args[1]), Path(a.adjoint_repo or EVO.parent / "stem-racing-adjoint"), a.allow_dirty)
        print(a.args[1])
        return
    c = cfg_load(a.config)
    if a.cmd == "submit":
        submit(a.args[0], c, a.allow_dirty)
    elif a.cmd == "status":
        status(c, a.args[0] if a.args else None)
    elif a.cmd == "wait":
        print(json.dumps(wait(c, a.args[0], a.poll), default=str))
    elif a.cmd == "fetch":
        fetch(c, a.args[0], a.args[1])


if __name__ == "__main__":
    main()
