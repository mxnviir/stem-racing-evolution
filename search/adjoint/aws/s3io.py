"""Tiny S3 helper for the container (boto3; the job's IAM role gives the access).
    s3io.py get <s3://bucket/key> <local file>
    s3io.py put <local file> <s3://bucket/key>
"""
import sys

import boto3


def split(uri):
    assert uri.startswith("s3://"), uri
    b, _, k = uri[5:].partition("/")
    return b, k


def main(a):
    s3 = boto3.client("s3")
    if a[1] == "get":
        b, k = split(a[2]); s3.download_file(b, k, a[3])
    elif a[1] == "put":
        b, k = split(a[3]); s3.upload_file(a[2], b, k)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv)
