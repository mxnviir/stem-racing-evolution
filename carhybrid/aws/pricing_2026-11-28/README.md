Snapshot behind carhybrid/aws/README.md section 5 (2026-11-28):
- `od_<region>.json`: on-demand Linux prices of every instance type with >= 256 GiB, from the AWS Price List API
  (`https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonEC2/current/<region>/index.csv`, publication
  2026-11-24), with vCPU, memory and processor.
- `model.py`: the time / cost model per 2.1 M-cell gradient, run in a folder that also holds AWS's public spot feed
  (`https://spot-price.s3.amazonaws.com/spot.js`) and the Spot Instance Advisor data
  (`https://spot-bid-advisor.s3.amazonaws.com/spot-advisor-data.json`); `model.json` is its output.
Live spot prices: `carhybrid/aws/spot_prices.sh`.
