"""
Canada's River Pulse 2025 - data download step (run inside QGIS Python).

Downloads, with the laptop's internet:
  1. GloFAS v4 daily river discharge (Open-Meteo Flood API) for ~190 sample reaches,
     after screening a 3x3 GloFAS neighbourhood to find the cell that is really the river.
  2. Daily mean 2 m air temperature for 2025 (NASA POWER, MERRA-2) on a 2-degree grid.

Everything runs in a background thread and is cached in raw/cache, so it can be
stopped and restarted. Progress is in builtins.RP_FETCH.
Usage inside QGIS:  exec(open(r"D:\\River Pulse\\fetch_river_pulse_data.py").read()); start()
"""
import os, json, time, math, threading, builtins, traceback, urllib.request, urllib.error
import numpy as np

ROOT = r"D:\River Pulse"
DATA = os.path.join(ROOT, "data")
CACHE = os.path.join(ROOT, "raw", "cache")
os.makedirs(CACHE, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/128.0 Safari/537.36"}
FLOOD = "https://flood-api.open-meteo.com/v1/flood"
YEAR0, YEAR1 = "2025-01-01", "2025-12-31"
SCREEN0, SCREEN1 = "2025-09-01", "2025-09-14"   # open-water, near-average flow window
BATCH = 15

if not hasattr(builtins, "RP_FETCH"):
    builtins.RP_FETCH = {"log": [], "phase": "idle", "done": False, "err": None}
ST = builtins.RP_FETCH


def log(msg):
    ST["log"].append(time.strftime("%H:%M:%S ") + msg)
    ST["log"] = ST["log"][-40:]


# ---------- polite rate limiting for Open-Meteo (600/min, 5000/h) ----------
_calls = getattr(builtins, 'RP_CALLS', None)
if _calls is None:
    _calls = builtins.RP_CALLS = []  # (time, weight), shared across re-runs


def _throttle(weight):
    while True:
        now = time.time()
        _calls[:] = [c for c in _calls if now - c[0] < 3600]
        m = sum(w for t, w in _calls if now - t < 60)
        h = sum(w for t, w in _calls)
        if m + weight <= 540 and h + weight <= 4800:
            _calls.append((now, weight))
            return
        time.sleep(5)


def get_json(url, weight=0):
    if weight:
        _throttle(weight)
    tries = 0
    while True:
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=180) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            body = e.read().decode(errors="ignore")
            if e.code == 429:
                if "aily" in body:
                    raise RuntimeError("Open-Meteo daily limit reached: " + body)
                wait = 600 if "ourly" in body else 65
                log(f"429 ({body[:60]}), waiting {wait}s")
                time.sleep(wait)
                continue
            tries += 1
            if tries > 4:
                raise
            log(f"HTTP {e.code}, retry"); time.sleep(10 * tries)
        except Exception as e:
            tries += 1
            if tries > 4:
                raise
            log(f"net error {e}, retry"); time.sleep(10 * tries)


def flood(lats, lons, d0, d1, tag):
    fn = os.path.join(CACHE, tag + ".json")
    if os.path.exists(fn):
        return json.load(open(fn))
    days = (np.datetime64(d1) - np.datetime64(d0)).astype(int) + 1
    w = len(lats) * math.ceil(days / 14)
    url = (f"{FLOOD}?latitude={','.join(f'{v:.3f}' for v in lats)}&longitude={','.join(f'{v:.3f}' for v in lons)}"
           f"&daily=river_discharge&start_date={d0}&end_date={d1}")
    d = get_json(url, w)
    if isinstance(d, dict):
        d = [d]
    res = [[(x if x is not None else np.nan) for x in r["daily"]["river_discharge"]] for r in d]
    json.dump(res, open(fn, "w"))
    return res


# ---------- sample selection ----------
def km(lo1, la1, lo2, la2):
    p = np.radians
    a = np.sin(p(la2 - la1) / 2) ** 2 + np.cos(p(la1)) * np.cos(p(la2)) * np.sin(p(lo2 - lo1) / 2) ** 2
    return 12742 * np.arcsin(np.sqrt(a))


def pick_samples(z):
    """Big sea outlets first (>=200 m3/s, 30 km apart), then every reach >=30 m3/s,
    largest first, thinned to 200 km apart."""
    dis, lon, lat, out = z["dis"], z["lon"], z["lat"], z["outlet"]
    sel = []
    big = np.where(out & (dis >= 200))[0]
    for i in big[np.argsort(-dis[big])]:
        if sel and km(lon[i], lat[i], lon[sel], lat[sel]).min() < 30:
            continue
        sel.append(i)
    cand = np.where(dis >= 30)[0]
    for i in cand[np.argsort(-dis[cand])]:
        if km(lon[i], lat[i], lon[sel], lat[sel]).min() < 200:
            continue
        sel.append(i)
    return np.array(sel)


def snap(v):
    return math.floor(v / 0.05) * 0.05 + 0.025


def run_flood():
    z = np.load(os.path.join(DATA, "reaches_attr.npz"))
    sel = pick_samples(z)
    ST["n_samples"] = len(sel)
    # 3x3 candidate cells per sample
    cells = []
    for i in sel:
        cy, cx = snap(z["lat"][i]), snap(z["lon"][i])
        cells += [(i, cy + dy * 0.05, cx + dx * 0.05) for dy in (-1, 0, 1) for dx in (-1, 0, 1)]
    ST["phase"] = "screening GloFAS cells"
    means = []
    for b in range(0, len(cells), BATCH):
        chunk = cells[b:b + BATCH]
        r = flood([c[1] for c in chunk], [c[2] for c in chunk], SCREEN0, SCREEN1, f"screen_{b:05d}")
        means += [np.nanmean(v) if np.isfinite(v).any() else np.nan for v in np.array(r, float)]
        ST["screen"] = f"{min(b + BATCH, len(cells))}/{len(cells)}"
    means = np.array(means).reshape(len(sel), 9)
    best = []
    for k, i in enumerate(sel):
        m = means[k]
        score = np.abs(np.log((m + 0.01) / z["dis"][i]))
        score[~np.isfinite(score)] = 99
        j = int(np.argmin(score))
        best.append((int(i), cells[k * 9 + j][1], cells[k * 9 + j][2], float(m[j]), float(score[j])))
    with open(os.path.join(DATA, "glofas_samples.csv"), "w") as f:
        f.write("reach_index,HYRIV_ID,DIS_AV_CMS,cell_lat,cell_lon,sep1_14_mean,abs_log_mismatch\n")
        for i, la, lo, m, s in best:
            f.write(f"{i},{z['ids'][i]},{z['dis'][i]:.2f},{la:.3f},{lo:.3f},{m:.2f},{s:.3f}\n")
    ST["phase"] = "downloading daily GloFAS 2025"
    Q = []
    for b in range(0, len(best), BATCH):
        chunk = best[b:b + BATCH]
        Q += flood([c[1] for c in chunk], [c[2] for c in chunk], YEAR0, YEAR1, f"year_{b:04d}")
        ST["year"] = f"{min(b + BATCH, len(best))}/{len(best)}"
    Q = np.array(Q, float)
    np.savez(os.path.join(DATA, "glofas_2025.npz"), reach_index=np.array([b[0] for b in best]),
             cell_lat=np.array([b[1] for b in best]), cell_lon=np.array([b[2] for b in best]), Q=Q)
    log(f"GloFAS done: {Q.shape}")


def run_temp():
    z = np.load(os.path.join(DATA, "reaches_attr.npz"))
    lon, lat = z["lon"], z["lat"]
    pts = []
    for la in np.arange(42, 84, 2.0):
        for lo in np.arange(-141, -51, 2.0):
            near = (np.abs(lat - la) < 1.5) & (np.abs(lon - lo) < 2.5)
            if near.any():
                pts.append((la, lo))
    ST["t_points"] = len(pts)
    T = []
    for k, (la, lo) in enumerate(pts):
        fn = os.path.join(CACHE, f"t2m_{la:.0f}_{lo:.0f}.json")
        if os.path.exists(fn):
            v = json.load(open(fn))
        else:
            u = ("https://power.larc.nasa.gov/api/temporal/daily/point?parameters=T2M&community=RE"
                 f"&longitude={lo}&latitude={la}&start=20241225&end=20251231&format=JSON")
            d = get_json(u)["properties"]["parameter"]["T2M"]
            v = [d[kk] if d[kk] > -900 else None for kk in sorted(d)]
            json.dump(v, open(fn, "w"))
        T.append([np.nan if x is None else x for x in v])
        ST["temp"] = f"{k + 1}/{len(pts)}"
    np.savez(os.path.join(DATA, "t2m_2025.npz"), lat=np.array([p[0] for p in pts]), lon=np.array([p[1] for p in pts]),
             T=np.array(T, float))
    log("temperature done")


def _wrap(fn, key):
    try:
        fn()
    except Exception:
        ST["err"] = traceback.format_exc()
        log("ERROR " + ST["err"][-300:])
    ST[key] = True


def start():
    ST.update({"done": False, "err": None, "flood_done": False, "temp_done": False})
    threading.Thread(target=_wrap, args=(run_flood, "flood_done"), daemon=True).start()
    threading.Thread(target=_wrap, args=(run_temp, "temp_done"), daemon=True).start()


def run_rescreen(th=math.log(3)):
    """Second pass for samples whose best 3x3 cell still disagrees with HydroRIVERS by more than 3x:
    search the next ring of GloFAS cells (5x5). If still no match, the sample is dropped so its reaches
    fall back to the next sampled reach downstream (or the nearest one)."""
    ST["phase"] = "re-screening poorly matched samples"
    z = np.load(os.path.join(DATA, "reaches_attr.npz"))
    rows = [l.strip().split(",") for l in open(os.path.join(DATA, "glofas_samples.csv"))][1:]
    g = dict(np.load(os.path.join(DATA, "glofas_2025.npz")))
    bad = [k for k, r in enumerate(rows) if float(r[6]) > th]
    cells = []
    for k in bad:
        i = int(rows[k][0])
        cy, cx = snap(z["lat"][i]), snap(z["lon"][i])
        cells += [(k, cy + dy * 0.05, cx + dx * 0.05) for dy in range(-2, 3) for dx in range(-2, 3)
                  if max(abs(dy), abs(dx)) == 2]
    means = []
    for b in range(0, len(cells), BATCH):
        chunk = cells[b:b + BATCH]
        r = flood([c[1] for c in chunk], [c[2] for c in chunk], SCREEN0, SCREEN1, f"rescreen_{b:05d}")
        means += [np.nanmean(v) if np.isfinite(v).any() else np.nan for v in np.array(r, float)]
        ST["rescreen"] = f"{min(b + BATCH, len(cells))}/{len(cells)}"
    means = np.array(means)
    fixed, dropped = [], []
    for k in bad:
        idx = [j for j, c in enumerate(cells) if c[0] == k]
        dis = float(rows[k][2])
        sc = np.abs(np.log((means[idx] + 0.01) / dis))
        sc[~np.isfinite(sc)] = 99
        j = idx[int(np.argmin(sc))]
        if sc.min() < min(th, float(rows[k][6])):
            rows[k][3], rows[k][4], rows[k][5], rows[k][6] = (f"{cells[j][1]:.3f}", f"{cells[j][2]:.3f}",
                                                              f"{means[j]:.2f}", f"{sc.min():.3f}")
            fixed.append(k)
        else:
            dropped.append(k)
    for b in range(0, len(fixed), BATCH):
        chunk = fixed[b:b + BATCH]
        q = flood([float(rows[k][3]) for k in chunk], [float(rows[k][4]) for k in chunk], YEAR0, YEAR1,
                  "reyear_" + "_".join(str(k) for k in chunk[:3]) + f"_{len(chunk)}")
        for k, v in zip(chunk, q):
            g["Q"][k] = np.array(v, float)
            g["cell_lat"][k], g["cell_lon"][k] = float(rows[k][3]), float(rows[k][4])
    for k in dropped:
        g["Q"][k] = np.nan
    np.savez(os.path.join(DATA, "glofas_2025.npz"), **g)
    with open(os.path.join(DATA, "glofas_samples.csv"), "w") as f:
        f.write("reach_index,HYRIV_ID,DIS_AV_CMS,cell_lat,cell_lon,sep1_14_mean,abs_log_mismatch,status\n")
        for k, r in enumerate(rows):
            stt = "moved (5x5)" if k in fixed else ("dropped" if k in dropped else "ok (3x3)")
            f.write(",".join(r[:7]) + f",{stt}\n")
    ST["rescreen_result"] = f"fixed {len(fixed)}, dropped {len(dropped)}"
    log(ST["rescreen_result"])


def start_rescreen():
    ST.update({"err": None, "rescreen_done": False})
    threading.Thread(target=_wrap, args=(run_rescreen, "rescreen_done"), daemon=True).start()
