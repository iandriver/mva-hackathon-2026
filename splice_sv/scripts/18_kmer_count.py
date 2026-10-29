#!/usr/bin/env python3
"""Count BUB1B target + control k-mers by streaming FASTQs, no alignment, no storage.

Each FASTQ is streamed from its URL (or read from a local path), gunzipped on the
fly, and its reads matched against the panel from script 17. Reads are never
written to disk, so the 85 GB release never lands locally.

  --probe URL [BYTES]   grab the first BYTES (default 150MB) of one file and
                        report throughput, for runtime estimation
  --files LIST          newline list of "name<TAB>url_or_path"; counts all and
                        writes per-window and per-control-region counts as JSON

Panel k-mers are stored in BOTH orientations so each read k-mer is a single
substring lookup with no per-k-mer canonicalization.
"""
import gzip, json, os, sys, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = os.path.join(HERE, "..", "alphagenome_control", "data", "kmer_panel.json")
OUT = os.environ.get("MVA_RESULTS", os.path.join(HERE, "..", "..", "results"))
COMP = str.maketrans("ACGT", "TGCA")

def load_panel():
    p = json.load(open(PANEL))
    k = p["k"]
    kmer_win = {}          # kmer(both strands) -> window_start (int)
    for w, kms in p["target"].items():
        for km in kms:
            rc = km.translate(COMP)[::-1]
            kmer_win[km] = int(w); kmer_win[rc] = int(w)
    ctrl = set()
    for km in p["control"]:
        ctrl.add(km); ctrl.add(km.translate(COMP)[::-1])
    return p, k, kmer_win, ctrl

def stream_reads(src, byte_limit=None):
    """Yield read sequences from a gzipped FASTQ URL or local path."""
    if src.startswith("http"):
        req = urllib.request.Request(src)
        req.add_header("Authorization", "Bearer " + open(
            os.path.expanduser("~/.cache/huggingface/token")).read().strip())
        if byte_limit:
            req.add_header("Range", f"bytes=0-{byte_limit-1}")
        raw = urllib.request.urlopen(req, timeout=120)
    else:
        raw = open(src, "rb")
    gz = gzip.GzipFile(fileobj=raw)
    line_no = 0
    try:
        for line in gz:
            line_no += 1
            if line_no % 4 == 2:
                yield line.rstrip()
    except (EOFError, OSError):
        return                     # truncated stream (probe): stop cleanly

def count(src, k, kmer_win, ctrl, win_counts, ctrl_counts, byte_limit=None):
    n = 0
    for seq in stream_reads(src, byte_limit):
        n += 1
        s = seq.decode("ascii", "ignore")
        for i in range(len(s) - k + 1):
            km = s[i:i+k]
            w = kmer_win.get(km)
            if w is not None:
                win_counts[w] = win_counts.get(w, 0) + 1
            elif km in ctrl:
                ctrl_counts[0] += 1
    return n

def main():
    p, k, kmer_win, ctrl = load_panel()
    print(f"panel: {len(kmer_win)//2:,} target + {len(ctrl)//2:,} control k-mers "
          f"(both strands loaded)", file=sys.stderr)

    if sys.argv[1] == "--probe":
        url = sys.argv[2]
        limit = int(sys.argv[3]) if len(sys.argv) > 3 else 150_000_000
        wc, cc = {}, [0]
        t0 = time.time()
        n = count(url, k, kmer_win, ctrl, wc, cc, byte_limit=limit)
        dt = time.time() - t0
        target_hits = sum(wc.values())
        print(f"probe: {n:,} reads in {dt:.1f}s = {n/dt:,.0f} reads/s "
              f"from {limit/1e6:.0f}MB compressed")
        print(f"  target hits {target_hits:,}  control hits {cc[0]:,}  "
              f"windows touched {len(wc)}/{len(p['windows'])}")
        # ~10.5GB/file, 8 files. reads scale with compressed bytes.
        full_reads = n * (84.67e9 / limit)
        print(f"  extrapolated full run: ~{full_reads/1e6:.0f}M reads, "
              f"~{full_reads/(n/dt)/3600:.1f} h single-threaded at this rate")
        return

    if sys.argv[1] == "--one":
        name, src, outjson = sys.argv[2], sys.argv[3], sys.argv[4]
        wc, cc = {}, [0]
        t0 = time.time()
        n = count(src, k, kmer_win, ctrl, wc, cc)
        json.dump(dict(name=name, n_reads=n,
                       window_counts={str(w): c for w, c in wc.items()},
                       control_total=cc[0], n_control_kmers=len(ctrl)//2),
                  open(outjson, "w"))
        print(f"{name}: {n:,} reads in {time.time()-t0:.0f}s -> {outjson}", flush=True)
        return

    if sys.argv[1] == "--files":
        files = [l.split("\t") for l in open(sys.argv[2]) if l.strip()]
        wc, cc = {}, [0]
        for name, src in files:
            t0 = time.time()
            n = count(src.strip(), k, kmer_win, ctrl, wc, cc)
            print(f"{name}: {n:,} reads in {time.time()-t0:.0f}s", flush=True)
        os.makedirs(OUT, exist_ok=True)
        out = os.path.join(OUT, "bub1b_kmer_counts.json")
        json.dump(dict(window_counts={str(w): c for w, c in wc.items()},
                       control_total=cc[0], n_control_kmers=len(ctrl)//2,
                       panel=os.path.basename(PANEL)), open(out, "w"))
        print(f"wrote {out}")

if __name__ == "__main__":
    main()
