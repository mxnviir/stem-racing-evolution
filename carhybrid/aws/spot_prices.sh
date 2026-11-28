#!/bin/bash
# Live spot prices (USD per hour, Linux) of the candidate instance types, per Availability Zone, from your AWS CLI
# (the authoritative number; the tables in README.md are a snapshot). On your computer:
#   bash carhybrid/aws/spot_prices.sh                     # us-east-2 and ap-south-1, the README's shortlist
#   REGIONS="us-east-2" TYPES="m7a.32xlarge r7a.16xlarge" bash carhybrid/aws/spot_prices.sh
set -uo pipefail
REGIONS=${REGIONS:-"us-east-2 ap-south-1"}
TYPES=${TYPES:-"m7a.32xlarge m7a.48xlarge r7a.16xlarge r7a.24xlarge m8a.48xlarge r8a.16xlarge m6a.32xlarge m6a.48xlarge r6a.16xlarge m8i.32xlarge m8i.48xlarge r8i.16xlarge r8i.24xlarge r7i.16xlarge c7a.48xlarge m7a.24xlarge"}
for R in $REGIONS; do
  echo "== $R (latest price per zone; a type missing here has no spot offer in $R)"
  # shellcheck disable=SC2086
  aws ec2 describe-spot-price-history --region "$R" --product-descriptions "Linux/UNIX" --instance-types $TYPES \
      --start-time "$(date -u +%FT%TZ)" --query 'SpotPriceHistory[].[InstanceType,AvailabilityZone,SpotPrice]' --output text \
    | sort -k1,1 -k3,3n | awk '{printf "  %-16s %-12s $%.3f\n", $1, $2, $3}'
done
echo "Interruption frequency per type and region: https://aws.amazon.com/ec2/spot/instance-advisor/"
