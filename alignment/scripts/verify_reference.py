"""Compare two FASTAs sequence-by-sequence (md5 of uppercased sequence); print effective genome size."""
import hashlib, sys

def read(p):
    d, n = {}, None
    for l in open(p):
        if l[0] == ">":
            n = l[1:].split()[0]; d[n] = []
        else:
            d[n].append(l.strip().upper())
    return {k: "".join(v) for k, v in d.items()}

a, b = read(sys.argv[1]), read(sys.argv[2])
h = lambda s: hashlib.md5(s.encode()).hexdigest()
hb = {h(v): k for k, v in b.items()}
ok = True
print(f"sequences: repo={len(a)} ncbi={len(b)}")
for k, v in a.items():
    m = hb.get(h(v))
    ok &= m is not None
    print(f"{k}\t{len(v)}\t{'identical to ' + m if m else 'NO MATCH'}")
tot = sum(len(v) for v in a.values())
nn = sum(v.count("N") for v in a.values())
print(f"total_bp={tot}\tN_bp={nn}\teffective_genome_size_nonN={tot-nn}")
print("RESULT:", "ALL IDENTICAL" if ok and len(a) == len(b) else "MISMATCH")
sys.exit(0 if ok else 1)
