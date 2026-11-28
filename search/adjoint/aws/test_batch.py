"""batch.py without AWS: packaging for real (git archive of both repositories and the shared core), and the S3 / Batch
calls against botocore's Stubber (the exact requests are checked, nothing is sent).

    python -m pytest search/adjoint/aws/test_batch.py -q      (needs ../stem-racing-adjoint; boto3)
"""
import io
import json
import tarfile
from pathlib import Path

import pytest

import batch as B

ADJ = B.EVO.parent / "stem-racing-adjoint"
pytestmark = pytest.mark.skipif(not (ADJ / ".git").exists(), reason="needs the adjoint repository next to this one")


def _job(tmp):
    j = tmp / "job.json"
    j.write_text(json.dumps({"id": "t1", "car": {"genes": [0.5] * 12}, "mesh": {"res": "coarse"}}))
    return j


def test_package_contents(tmp_path):
    out = tmp_path / "in.tgz"
    job = B.package(_job(tmp_path), out, ADJ, allow_dirty=True)
    names = tarfile.open(out).getnames()
    assert "job.json" in names
    assert any(n.startswith("evolution/search/adjoint/fine_adjoint.py") for n in names)
    assert any(n.startswith("adjoint/research/sidepod_opt/dafoam/runScript.py") for n in names)
    assert any(n.startswith("adjoint/research/sidepod_opt/prodcheck.py") for n in names)
    assert "core/carhybrid/grad.py" in names and "core/carhybrid/patch_daf.py" in names
    assert set(job["code"]) == {"evolution", "adjoint", "core"} and job["code"]["core"] == B.CORE_REF
    packed = json.load(tarfile.open(out).extractfile("job.json"))
    assert packed["code"] == job["code"]


def test_submit_status_fetch(tmp_path):
    boto3 = pytest.importorskip("boto3")
    from botocore.stub import Stubber, ANY
    c = {"region": "us-east-2", "bucket": "bkt", "prefix": "fa", "job_queue": "q", "job_definition": "jd",
         "vcpus": 64, "memory_mib": 480000, "np": 64, "ledger": str(tmp_path / "ledger.jsonl"), "adjoint_repo": str(ADJ)}
    s3 = boto3.client("s3", region_name="us-east-2", aws_access_key_id="x", aws_secret_access_key="x")
    bt = boto3.client("batch", region_name="us-east-2", aws_access_key_id="x", aws_secret_access_key="x")
    with Stubber(s3) as ss, Stubber(bt) as sb:
        ss.add_response("put_object", {}, {"Bucket": "bkt", "Key": ANY, "Body": ANY})
        sb.add_response("submit_job", {"jobId": "J1", "jobName": "fine-adjoint-t1", "jobArn": "arn:j"},
                        {"jobName": "fine-adjoint-t1", "jobQueue": "q", "jobDefinition": "jd", "containerOverrides": ANY})
        row = B.submit(_job(tmp_path), c, allow_dirty=True, s3=s3, batch=bt)
        assert row["job_id"] == "J1" and row["key"].startswith("fa/t1-")
        sb.add_response("describe_jobs", {"jobs": [{"jobId": "J1", "jobName": "n", "jobQueue": "q", "status": "SUCCEEDED",
                                                    "startedAt": 1, "jobDefinition": "jd"}]}, {"jobs": ["J1"]})
        st = B.status(c, "J1", batch=bt)
        assert st[0]["status"] == "SUCCEEDED"
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tf:
            b = json.dumps({"usable": True}).encode(); ti = tarfile.TarInfo("gradient.json"); ti.size = len(b)
            tf.addfile(ti, io.BytesIO(b))
        from botocore.response import StreamingBody
        for name, body in (("out.tgz", buf.getvalue()), ("result.json", json.dumps({"id": "t1", "usable": True}).encode())):
            ss.add_response("get_object", {"Body": StreamingBody(io.BytesIO(body), len(body))}, {"Bucket": "bkt", "Key": f"{row['key']}/{name}"})
        ss.add_client_error("get_object", "NoSuchKey", expected_params={"Bucket": "bkt", "Key": f"{row['key']}/fine_adjoint.log"})
        ss.add_client_error("get_object", "NoSuchKey", expected_params={"Bucket": "bkt", "Key": f"{row['key']}/state.tgz"})
        got = B.fetch(c, "J1", tmp_path / "res", s3=s3)
        assert got == ["out.tgz", "result.json"] and (tmp_path / "res" / "gradient.json").exists()
