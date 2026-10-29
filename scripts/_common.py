"""Shared helpers.

The gnomAD wrapper exists because "absent from gnomAD" is evidence -- it is what supports
PM2 for a private variant. A failed query must therefore never be reported as absence.
Every result carries an explicit status, and callers are expected to branch on it.
"""
import json, time, urllib.error, urllib.request

GNOMAD = "https://gnomad.broadinstitute.org/api"

def gnomad_variant(variant_id, dataset="gnomad_r4", retries=4, backoff=3.0, timeout=90):
    """Returns (status, payload).

    status is one of:
      "ok"       payload is the variant record (present in gnomAD)
      "absent"   the query succeeded and gnomAD holds no such variant; payload is None
      "failed"   the query did not complete; payload is the reason string. NOT absence.

    payload is a dict only when status == "ok". Callers must branch on status before
    treating it as one.
    """
    q = ('{variant(variantId:"%s",dataset:%s){variant_id rsids '
         'exome{ac an af} genome{ac an af} in_silico_predictors{id value}}}'
         % (variant_id, dataset))
    body = json.dumps({"query": q}).encode()
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                GNOMAD, data=body,
                headers={"Content-Type": "application/json",
                         "User-Agent": "mva-hackathon-2026/1.0"})
            d = json.load(urllib.request.urlopen(req, timeout=timeout))
            errs = d.get("errors") or []
            if errs:
                # A server-side timeout is a failure, not an answer about the variant.
                last = "; ".join(e.get("message", "?") for e in errs)
                if any("not found" in (e.get("message") or "").lower() for e in errs):
                    return "absent", None
                time.sleep(backoff * (attempt + 1)); continue
            v = (d.get("data") or {}).get("variant")
            return ("ok", v) if v else ("absent", None)
        except (urllib.error.URLError, TimeoutError, OSError, ValueError) as e:
            last = f"{type(e).__name__}: {e}"
            time.sleep(backoff * (attempt + 1))
    return "failed", last

def af_of(status, payload):
    """Highest AF across exome/genome, or None when unknown. Callers must check status."""
    if status != "ok" or not payload: return None
    return max([(payload.get(k) or {}).get("af") or 0.0 for k in ("exome", "genome")])

def fmt_af(status, payload):
    if status == "failed": return "QUERY FAILED"
    if status == "absent": return "absent"
    af = af_of(status, payload)
    return f"{af:.2e}" if af else "absent"
