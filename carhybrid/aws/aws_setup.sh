#!/bin/bash
# One-time AWS setup for the fine-mesh DAFoam runs (carhybrid/aws/README.md). Run on YOUR computer with the AWS CLI
# v2 configured (`aws configure` or AWS_PROFILE). Idempotent: re-running it finds what already exists.
#
#   REGION=us-east-1 BUCKET=<globally unique name> bash carhybrid/aws/aws_setup.sh
#
# Creates (all named carhybrid-dafoam unless overridden):
#   - an S3 bucket (BUCKET; public access blocked, default encryption) that carries the code in and the results out;
#   - an IAM role + instance profile whose only permission is read/write on that bucket (plus, with SSM=1, the AWS
#     managed policy for Session Manager, to open a shell on the instance without SSH);
#   - a security group in the default VPC (or VPC_ID) with NO inbound rules (the instance only makes outbound HTTPS).
# Writes the names to ~/.config/carhybrid/aws.env, which launch.sh / status.sh / fetch.sh read. Nothing secret is
# written anywhere: the instance gets its permissions from the role, not from keys.
set -euo pipefail
REGION=${REGION:-$(aws configure get region || true)}
: "${REGION:?set REGION (e.g. us-east-1)}"
: "${BUCKET:?set BUCKET to a globally unique S3 bucket name}"
NAME=${NAME:-carhybrid-dafoam}
SSM=${SSM:-1}
CONF=${CONF:-$HOME/.config/carhybrid/aws.env}
export AWS_DEFAULT_REGION=$REGION
say() { echo "[aws_setup] $*"; }

aws sts get-caller-identity --query Account --output text > /dev/null || { echo "AWS CLI is not configured"; exit 1; }

# ---------------------------------------------------------------- bucket
if aws s3api head-bucket --bucket "$BUCKET" 2>/dev/null; then
  say "bucket $BUCKET exists"
else
  if [ "$REGION" = us-east-1 ]; then aws s3api create-bucket --bucket "$BUCKET" > /dev/null
  else aws s3api create-bucket --bucket "$BUCKET" --create-bucket-configuration LocationConstraint="$REGION" > /dev/null; fi
  say "bucket $BUCKET created"
fi
aws s3api put-public-access-block --bucket "$BUCKET" --public-access-block-configuration \
  BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=true,RestrictPublicBuckets=true
aws s3api put-bucket-encryption --bucket "$BUCKET" --server-side-encryption-configuration \
  '{"Rules":[{"ApplyServerSideEncryptionByDefault":{"SSEAlgorithm":"AES256"}}]}'

# ---------------------------------------------------------------- role and instance profile
TRUST='{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"ec2.amazonaws.com"},"Action":"sts:AssumeRole"}]}'
if ! aws iam get-role --role-name "$NAME" > /dev/null 2>&1; then
  aws iam create-role --role-name "$NAME" --assume-role-policy-document "$TRUST" \
    --description "carhybrid fine-mesh DAFoam instances: read/write one S3 bucket" > /dev/null
  say "role $NAME created"
fi
POLICY=$(cat <<EOF
{"Version":"2012-10-17","Statement":[
 {"Effect":"Allow","Action":["s3:ListBucket"],"Resource":"arn:aws:s3:::$BUCKET"},
 {"Effect":"Allow","Action":["s3:GetObject","s3:PutObject"],"Resource":"arn:aws:s3:::$BUCKET/*"}]}
EOF
)
aws iam put-role-policy --role-name "$NAME" --policy-name "$NAME-s3" --policy-document "$POLICY"
if [ "$SSM" = 1 ]; then
  aws iam attach-role-policy --role-name "$NAME" --policy-arn arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore
fi
if ! aws iam get-instance-profile --instance-profile-name "$NAME" > /dev/null 2>&1; then
  aws iam create-instance-profile --instance-profile-name "$NAME" > /dev/null
  aws iam add-role-to-instance-profile --instance-profile-name "$NAME" --role-name "$NAME"
  say "instance profile $NAME created (IAM can take ~10 s to make it usable)"
fi

# ---------------------------------------------------------------- security group (no inbound)
VPC_ID=${VPC_ID:-$(aws ec2 describe-vpcs --filters Name=isDefault,Values=true --query 'Vpcs[0].VpcId' --output text)}
[ "$VPC_ID" != None ] || { echo "no default VPC in $REGION: set VPC_ID and SUBNET_ID (a public subnet)"; exit 1; }
SG_ID=$(aws ec2 describe-security-groups --filters Name=group-name,Values="$NAME" Name=vpc-id,Values="$VPC_ID" \
        --query 'SecurityGroups[0].GroupId' --output text)
if [ "$SG_ID" = None ]; then
  SG_ID=$(aws ec2 create-security-group --group-name "$NAME" --vpc-id "$VPC_ID" \
          --description "carhybrid DAFoam: no inbound, outbound only" --query GroupId --output text)
  say "security group $SG_ID created (no inbound rules)"
fi

mkdir -p "$(dirname "$CONF")"
cat > "$CONF" <<EOF
# written by carhybrid/aws/aws_setup.sh $(date -u +%FT%TZ)
REGION=$REGION
BUCKET=$BUCKET
PROFILE_NAME=$NAME
SG_ID=$SG_ID
SUBNET_ID=${SUBNET_ID:-}
KEY_NAME=${KEY_NAME:-}
EOF
say "config written to $CONF:"; cat "$CONF"
say "check your EC2 quota for the instance family before the first launch: Service Quotas -> EC2 ->"
say "'Running On-Demand Standard (A, C, D, H, I, M, R, T, Z) instances' must be >= 64 vCPU (r7i.16xlarge)."
