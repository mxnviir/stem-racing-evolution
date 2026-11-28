import json
adv=json.load(open("spot-advisor-data.json")); A=adv["spot_advisor"]; rng={x["index"]:x["label"] for x in adv["ranges"]}
s=open("9a826e60.dat").read(); j=json.loads(s[s.index("(")+1:s.rindex(")")])
feed={}
for r in j["config"]["regions"]:
    for t in r["instanceTypes"]:
        for sz in t["sizes"]:
            v=[x for x in sz["valueColumns"] if x["name"]=="linux"][0]["prices"]["USD"]
            try: feed[(r["region"],sz["size"])]=float(v)
            except: pass
# per physical core speed vs the hosted runner's Milan core (sparse CFD: memory bandwidth + IPC), rough
SPEED={"7R13":1.0,"8375C":1.0,"Icelake":1.0,"Ice Lake":1.0,"9R14":1.35,"9R45":1.55,"Sapphire":1.3,"Granite":1.5,"9575F":1.6}
def speed(cpu):
    for k,v in SPEED.items():
        if k in cpu: return v
    return None
def eff(p):  # parallel + bandwidth-saturation efficiency at 2.1 M cells
    return 0.9 if p<=64 else 0.8 if p<=96 else 0.75 if p<=128 else 0.65
# work for ONE 2.1 M-cell gradient, Milan-core-hours (from the hosted measurements), + serial hours
WORK={"low":48,"mid":81,"high":161}; SERIAL=0.8
FAMS=("r6a","r6i","r7a","r7i","r8a","r8i","r8i-flex","m6a","m6i","m7a","m7i","m8a","m8i","c7a","c8a","c7i","c8i","hpc6a","hpc7a","x2idn","x8i","r6id","r6in","m8azn","r7iz")
rows=[]
for reg in ("ap-south-1","us-east-2"):
    od=json.load(open(f"od_{reg}.json"))
    for t,v in od.items():
        fam=t.split(".")[0]
        if fam not in FAMS or v["mem"]<380 or v["mem"]>1100: continue
        smt=not(fam.endswith("7a") or fam.endswith("8a") or fam.startswith("hpc") or fam=="m8azn")
        cores=v["vcpu"]//2 if smt else v["vcpu"]
        sp=speed(v["cpu"]); 
        if sp is None: continue
        eq=cores*sp*(1.15 if smt else 1.0)*eff(cores)
        a=A.get(reg,{}).get("Linux",{}).get(t)
        spot=feed.get((reg,t)); src="feed"
        if spot is None and a: spot=v["usd"]*(1-a["s"]/100); src="adv"
        h={k:SERIAL+w/eq for k,w in WORK.items()}
        rows.append(dict(reg=reg,t=t,cores=cores,vcpu=v["vcpu"],mem=v["mem"],cpu=v["cpu"],od=v["usd"],spot=spot,src=src,
                         intr=(rng.get(a["r"]) if a else "no spot"),h=h,
                         cod=v["usd"]*h["mid"],csp=(spot*h["mid"] if spot else None)))
json.dump(rows,open("model.json","w"),indent=0)
for reg in ("us-east-2","ap-south-1"):
  for memmin,label in ((480,"512+ GiB (2.1 M cells safe)"),(380,"384 GiB (2.1 M only if F0 measures <~350 GB)")):
    R=[r for r in rows if r["reg"]==reg and (r["mem"]>=480 if memmin==480 else r["mem"]<480)]
    R.sort(key=lambda r:(r["csp"] if r["csp"] else 9e9))
    print(f"==== {reg}  {label}")
    for r in R[:14]:
        print(f"{r['t']:17s} {r['cores']:4d}c {r['mem']:5.0f}G {r['cpu'][:22]:22s} h {r['h']['low']:.1f}-{r['h']['mid']:.1f}-{r['h']['high']:.1f}  od ${r['od']:6.2f}/h -> ${r['cod']:6.1f}/grad  spot ${(r['spot'] or 0):5.2f}/h({r['src']}) -> ${(r['csp'] or 0):5.1f}/grad  int {r['intr']}")
