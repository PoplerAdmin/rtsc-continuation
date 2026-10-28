import json,sys,glob
for f in sorted(glob.glob(sys.argv[1])):
    d=json.load(open(f)); print("==",f.split('/')[-1],"wall",round(d['wall_s'],1))
    for r in d['rows']:
        if len(sys.argv)>2 and not r['label'].startswith(sys.argv[2]): continue
        print(r['label']); print("  "+"\n  ".join(r['notes'][3:]))
