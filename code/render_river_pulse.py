"""
Canada's River Pulse 2025 - rebuild layers from saved data and render frames.

Run inside the QGIS Python console (QGIS 4 / Qt6):

    exec(open(r"D:\\River Pulse\\render_river_pulse.py", encoding="utf-8").read())
    setup()                     # builds the daily series + map layers (about 1 min the first time)
    render_range(0, 364)        # renders frames/day_000.png ... day_364.png  (any range works)
    make_project()              # saves a QGIS project showing the peak day

Inputs (all in data/):
    rivers_canada.gpkg     HydroRIVERS v1.0 reaches (North America + Arctic) inside Canada, DIS_AV_CMS >= 1
    reaches_attr.npz       reach attributes, midpoints, downstream links, sea-outlet flags
    glofas_2025.npz        GloFAS v4 daily discharge 2025 at the sample reaches (Open-Meteo Flood API)
    t2m_2025.npz           NASA POWER (MERRA-2) daily 2 m air temperature, 2-degree grid
Outputs:
    data/discharge_to_sea_2025.csv, frames/day_XXX.png, Canada_River_Pulse_2025_peak_day.qgz
"""
import os, csv, math, builtins, datetime as dt
import numpy as np
from qgis.core import (Qgis, QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry, QgsField, QgsFields,
                       QgsCoordinateReferenceSystem, QgsCoordinateTransform, QgsMapSettings,
                       QgsMapRendererCustomPainterJob, QgsLineSymbol, QgsFillSymbol, QgsSingleSymbolRenderer,
                       QgsGeometryGeneratorSymbolLayer, QgsSimpleFillSymbolLayer, QgsSimpleLineSymbolLayer,
                       QgsSymbolLayer, QgsProperty, QgsFeatureRequest, QgsExpressionContextUtils,
                       QgsVectorFileWriter, QgsRectangle, QgsPointXY)
from qgis.PyQt.QtCore import Qt, QSize, QRectF, QPointF, QMetaType
from qgis.PyQt.QtGui import (QImage, QPainter, QColor, QFont, QPen, QBrush, QLinearGradient, QPainterPath,
                             QFontMetricsF, QTransform)

# ----------------------------------------------------------------- settings
ROOT = r"D:\River Pulse"
DATA = os.path.join(ROOT, "data")
FRAMES = os.path.join(ROOT, "frames")
W, H = 1080, 1350
MAP_CRS = "EPSG:3978"          # Canada Atlas Lambert (the Bangladesh run used UTM 46N)
TILT = 0.62                    # vertical squash for the tilted view
K_H = 4000.0                   # wall height in metres = K_H * sqrt(Q)
MAP_BOTTOM = 885               # pixel row where the southern edge of the land sits
T_MIN, T_MAX = 0.0, 24.0       # legend range, degrees C
DAY0 = dt.date(2025, 1, 1)
NDAYS = 365

BG = QColor("#F3E5EF")
LAND = QColor("#FBF4F8")
LAND_EDGE = QColor("#D8C2D3")
SHADOW = QColor("#E4CFDE")
INK = QColor("#3A2940")
MUTED = QColor("#7D6A84")
ACCENT = QColor("#2E6F95")
RAMP = [(0.0, "#1B9AAA"), (0.3, "#7FC8B8"), (0.55, "#F2D06B"), (0.78, "#F28C38"), (1.0, "#D6311F")]
FONT = "Segoe UI"


def day_date(d):
    return DAY0 + dt.timedelta(days=int(d))


# ----------------------------------------------------------------- daily series
def _km(lo1, la1, lo2, la2):
    p = np.radians
    a = np.sin(p(la2 - la1) / 2) ** 2 + np.cos(p(la1)) * np.cos(p(la2)) * np.sin(p(lo2 - lo1) / 2) ** 2
    return 12742 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def _fill_nan(a):
    a = a.copy()
    x = np.arange(a.shape[1])
    for r in a:
        ok = np.isfinite(r)
        if ok.any() and not ok.all():
            r[~ok] = np.interp(x[~ok], x[ok], r[ok])
    return a


def build_series(force=False):
    cache = os.path.join(DATA, "series_cache.npz")
    if os.path.exists(cache) and not force:
        return dict(np.load(cache))
    z = np.load(os.path.join(DATA, "reaches_attr.npz"))
    g = np.load(os.path.join(DATA, "glofas_2025.npz"))
    t = np.load(os.path.join(DATA, "t2m_2025.npz"))
    n = len(z["ids"])
    Q = _fill_nan(g["Q"][:, :NDAYS])
    mean = np.nanmean(Q, axis=1)
    valid = np.isfinite(Q).all(axis=1) & (mean > 0.01)
    ratio = Q[valid] / mean[valid, None]
    s_reach = g["reach_index"][valid]
    sid_of = np.full(n, -1)
    sid_of[s_reach] = np.arange(len(s_reach))
    # 1) first sampled reach downstream
    down = z["down"]
    assign = np.full(n, -2)
    for i in range(n):
        if assign[i] != -2:
            continue
        path, c = [], i
        while c >= 0 and assign[c] == -2 and sid_of[c] < 0:
            path.append(c)
            c = down[c]
        val = -1 if c < 0 else (sid_of[c] if sid_of[c] >= 0 else assign[c])
        if c >= 0 and sid_of[c] >= 0:
            assign[c] = sid_of[c]
        for p in path:
            assign[p] = val
    via_down = int((assign >= 0).sum())
    # 2) no sample downstream (e.g. the reaches just above a river mouth): use the biggest sampled
    #    reach of the same river system (MAIN_RIV), so a river mouth follows its own main stem
    main = z["main"]
    best_in_basin = {}
    for k, i in enumerate(s_reach):
        m = int(main[i])
        if m not in best_in_basin or z["dis"][i] > z["dis"][s_reach[best_in_basin[m]]]:
            best_in_basin[m] = k
    for i in np.where(assign < 0)[0]:
        k = best_in_basin.get(int(main[i]))
        if k is not None:
            assign[i] = k
    via_basin = int((assign >= 0).sum()) - via_down
    # 3) otherwise the nearest sampled reach
    lon, lat = z["lon"], z["lat"]
    miss = np.where(assign < 0)[0]
    slon, slat = lon[s_reach], lat[s_reach]
    for k in range(0, len(miss), 5000):
        m = miss[k:k + 5000]
        d = _km(lon[m, None], lat[m, None], slon[None, :], slat[None, :])
        assign[m] = np.argmin(d, axis=1)
    # temperature: 7-day trailing mean of air temperature, IDW (3 nearest grid points)
    T = _fill_nan(t["T"])          # days from 2024-12-25
    T7 = np.stack([T[:, d + 1:d + 8].mean(axis=1) for d in range(NDAYS)], axis=1)
    WT = np.maximum(T7, 0.0)       # water cannot be colder than freezing
    tidx = np.zeros((n, 3), int)
    tw = np.zeros((n, 3))
    for k in range(0, n, 5000):
        sl = slice(k, k + 5000)
        d = _km(lon[sl, None], lat[sl, None], t["lon"][None, :], t["lat"][None, :])
        o = np.argsort(d, axis=1)[:, :3]
        dd = np.take_along_axis(d, o, axis=1)
        w = 1.0 / np.maximum(dd, 1.0) ** 2
        tidx[sl] = o
        tw[sl] = w / w.sum(axis=1, keepdims=True)
    out = z["outlet"]
    total = (z["dis"][out, None] * ratio[assign[out]]).sum(axis=0)
    res = dict(assign=assign, ratio=ratio, WT=WT, tidx=tidx, tw=tw, total=total, s_reach=s_reach,
               via_down=np.array(via_down), via_basin=np.array(via_basin))
    np.savez(cache, **res)
    with open(os.path.join(DATA, "discharge_to_sea_2025.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["date", "discharge_to_sea_m3s"])
        for d in range(NDAYS):
            wr.writerow([day_date(d).isoformat(), round(float(total[d]), 1)])
    return res


# ----------------------------------------------------------------- layers
def _tilt_transform():
    return QTransform.fromScale(1.0, TILT)


def _ramp_expr(field):
    stops = ",".join(f"{p},'{c}'" for p, c in RAMP)
    return f"ramp_color(create_ramp(map({stops})), scale_linear({field}, {T_MIN}, {T_MAX}, 0, 1))"


def wall_renderer(h_field='"h"', t_field='"wt"', q_field='"q"'):
    sym = QgsLineSymbol()
    sym.deleteSymbolLayer(0)
    # 1) the wall face
    fill = QgsFillSymbol()
    fl = fill.symbolLayer(0)
    fl.setStrokeStyle(Qt.PenStyle.NoPen)
    fl.setDataDefinedProperty(QgsSymbolLayer.Property.FillColor,
                              QgsProperty.fromExpression(
                                  f"set_color_part({_ramp_expr(t_field)}, 'alpha', "
                                  f"round(scale_linear(log10(max({q_field}, 1)), 0, 3.5, 35, 240)))"))
    gg = QgsGeometryGeneratorSymbolLayer.create({"geometryModifier": f"extrude($geometry, 0, {h_field})"})
    gg.setSymbolType(Qgis.SymbolType.Fill)
    gg.setSubSymbol(fill)
    sym.appendSymbolLayer(gg)
    # 2) a crisp top edge
    top = QgsLineSymbol()
    tl = top.symbolLayer(0)
    tl.setWidth(0.18)
    tl.setWidthUnit(Qgis.RenderUnit.Millimeters)
    tl.setDataDefinedProperty(QgsSymbolLayer.Property.StrokeColor,
                              QgsProperty.fromExpression(
                                  f"set_color_part(darker({_ramp_expr(t_field)}, 135), 'alpha', "
                                  f"230)"))
    gt = QgsGeometryGeneratorSymbolLayer.create(
        {"geometryModifier": f"if({q_field} >= 30, translate($geometry, 0, {h_field}), NULL)"})
    gt.setSymbolType(Qgis.SymbolType.Line)
    gt.setSubSymbol(top)
    sym.appendSymbolLayer(gt)
    r = QgsSingleSymbolRenderer(sym)
    r.setOrderBy(QgsFeatureRequest.OrderBy([QgsFeatureRequest.OrderByClause("ymid", False)]))
    r.setOrderByEnabled(True)
    return r


def _land_layer(name, color, edge, dy):
    lyr = QgsVectorLayer(f"Polygon?crs={MAP_CRS}", name, "memory")
    b = QgsVectorLayer(os.path.join(ROOT, "raw", "CAN_ADM0_simplified.geojson"), "b", "ogr")
    xf = QgsCoordinateTransform(b.crs(), QgsCoordinateReferenceSystem(MAP_CRS), QgsProject.instance())
    feats = []
    for f in b.getFeatures():
        g = QgsGeometry(f.geometry())
        g.transform(xf)
        g = g.simplify(1500)
        g.transform(_tilt_transform())
        g.translate(0, dy)
        nf = QgsFeature()
        nf.setGeometry(g)
        feats.append(nf)
    lyr.dataProvider().addFeatures(feats)
    s = QgsFillSymbol.createSimple({"color": color.name(), "outline_color": edge.name(), "outline_width": "0.25",
                                    "outline_width_unit": "MM"})
    lyr.setRenderer(QgsSingleSymbolRenderer(s))
    return lyr


def prepare_walls(force=False):
    """Reproject, simplify and tilt the river reaches once, with fast QGIS processing tools."""
    import processing
    out = os.path.join(DATA, "walls_prepared.gpkg")
    if os.path.exists(out) and not force:
        return out
    src = os.path.join(DATA, "rivers_canada.gpkg") + "|layername=rivers"
    r1 = processing.run("native:reprojectlayer", {"INPUT": src, "TARGET_CRS": QgsCoordinateReferenceSystem(MAP_CRS),
                                                  "OUTPUT": "TEMPORARY_OUTPUT"})["OUTPUT"]
    r2 = processing.run("native:simplifygeometries", {"INPUT": r1, "METHOD": 0, "TOLERANCE": 400,
                                                      "OUTPUT": "TEMPORARY_OUTPUT"})["OUTPUT"]
    r3 = processing.run("native:affinetransform", {"INPUT": r2, "DELTA_X": 0, "DELTA_Y": 0, "DELTA_Z": 0, "DELTA_M": 0,
                                                   "SCALE_X": 1, "SCALE_Y": TILT, "SCALE_Z": 1, "SCALE_M": 1,
                                                   "ROTATION_Z": 0, "OUTPUT": "TEMPORARY_OUTPUT"})["OUTPUT"]
    processing.run("native:fieldcalculator", {"INPUT": r3, "FIELD_NAME": "ymid", "FIELD_TYPE": 0, "FIELD_LENGTH": 20,
                                              "FIELD_PRECISION": 3, "FORMULA": "y(centroid($geometry))",
                                              "OUTPUT": out})
    return out


def make_layers():
    z = np.load(os.path.join(DATA, "reaches_attr.npz"))
    pos = {int(h): k for k, h in enumerate(z["ids"])}
    src = QgsVectorLayer(prepare_walls() + "|layername=walls_prepared", "src", "ogr")
    req = QgsFeatureRequest().setSubsetOfAttributes(["HYRIV_ID", "DIS_AV_CMS", "ymid"], src.fields())
    lyr = src.materialize(req)
    lyr.setName("River walls")
    pr = lyr.dataProvider()
    pr.addAttributes([QgsField("ridx", QMetaType.Type.Int), QgsField("h", QMetaType.Type.Double),
                      QgsField("wt", QMetaType.Type.Double), QgsField("q", QMetaType.Type.Double)])
    lyr.updateFields()
    ir = lyr.fields().indexOf("ridx")
    fids, ridx = [], []
    for f in lyr.getFeatures(QgsFeatureRequest().setFlags(QgsFeatureRequest.Flag.NoGeometry)):
        fids.append(f.id())
        ridx.append(pos[int(f["HYRIV_ID"])])
    pr.changeAttributeValues({fid: {ir: r} for fid, r in zip(fids, ridx)})
    lyr.setRenderer(wall_renderer())
    land = _land_layer("Canada", LAND, LAND_EDGE, 0)
    shadow = _land_layer("Canada shadow", SHADOW, SHADOW, -22000)
    return dict(walls=lyr, land=land, shadow=shadow, fids=np.array(fids), ridx=np.array(ridx))


def setup(force=False):
    builtins.RP = getattr(builtins, "RP", {})
    RP = builtins.RP
    RP["ser"] = build_series(force)
    RP["z"] = dict(np.load(os.path.join(DATA, "reaches_attr.npz")))
    if force or "walls" not in RP:
        RP.update(make_layers())
    tot = RP["ser"]["total"]
    RP["lo"], RP["hi"] = int(np.argmin(tot)), int(np.argmax(tot))
    ext = QgsRectangle()
    for L in (RP["land"], RP["walls"]):
        e = L.extent()
        ext = e if ext.isEmpty() else (ext.combineExtentWith(e) or ext)
    RP["land_ext"] = RP["land"].extent()
    return RP


# ----------------------------------------------------------------- per-day values
def day_values(d):
    RP = builtins.RP
    s, z = RP["ser"], RP["z"]
    q = z["dis"] * s["ratio"][s["assign"], d]
    wt = (s["WT"][s["tidx"], d] * s["tw"]).sum(axis=1)
    return q, wt


def set_day(d):
    RP = builtins.RP
    q, wt = day_values(d)
    lyr = RP["walls"]
    ih, iw, iq = lyr.fields().indexOf("h"), lyr.fields().indexOf("wt"), lyr.fields().indexOf("q")
    r = RP["ridx"]
    hq = K_H * np.sqrt(np.maximum(q[r], 0))
    ch = {int(f): {ih: float(a), iw: float(b), iq: float(c)} for f, a, b, c in zip(RP["fids"], hq, wt[r], q[r])}
    lyr.dataProvider().changeAttributeValues(ch)


# ----------------------------------------------------------------- map + overlay
def map_settings():
    RP = builtins.RP
    e = RP["land_ext"]
    margin = 22
    s = e.width() / (W - 2 * margin)                  # metres per pixel
    ymin = e.yMinimum() - (H - MAP_BOTTOM) * s
    xmin = e.xMinimum() - margin * s
    ext = QgsRectangle(xmin, ymin, xmin + W * s, ymin + H * s)
    ms = QgsMapSettings()
    ms.setLayers([RP["walls"], RP["land"], RP["shadow"]])
    ms.setDestinationCrs(QgsCoordinateReferenceSystem(MAP_CRS))
    ms.setOutputSize(QSize(W, H))
    ms.setOutputDpi(96)
    ms.setExtent(ext)
    ms.setBackgroundColor(BG)
    ms.setFlag(Qgis.MapSettingsFlag.Antialiasing, True)
    return ms


def _font(size, weight=QFont.Weight.Normal, spacing=0.0):
    f = QFont(FONT)
    f.setPixelSize(int(size))
    f.setWeight(weight)
    if spacing:
        f.setLetterSpacing(QFont.SpacingType.AbsoluteSpacing, spacing)
    return f


def _text(p, x, y, s, font, color, align=Qt.AlignmentFlag.AlignLeft, halo=None):
    p.setFont(font)
    fm = QFontMetricsF(font)
    w = fm.horizontalAdvance(s)
    if align == Qt.AlignmentFlag.AlignRight:
        x -= w
    elif align == Qt.AlignmentFlag.AlignHCenter:
        x -= w / 2
    if halo is not None:
        path = QPainterPath()
        path.addText(QPointF(x, y), font, s)
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.setPen(QPen(halo, 5, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin))
        p.drawPath(path)
    p.setPen(QPen(color))
    p.drawText(QPointF(x, y), s)
    return w


def _fmt(v):
    return f"{v:,.0f}"


def draw_overlay(p, d):
    RP = builtins.RP
    tot = RP["ser"]["total"]
    lo, hi = RP["lo"], RP["hi"]
    L, R = 56, W - 56
    # --- title block
    _text(p, L, 92, "Canada's River Pulse", _font(54, QFont.Weight.Bold), INK)
    _text(p, L, 132, "EVERY RIVER, EVERY DAY OF 2025", _font(19, QFont.Weight.DemiBold, 3.2), MUTED)
    _text(p, L, 164, "Map: Chinmoy Ghosh Shuvo", _font(16), MUTED)
    # --- big date (top right)
    dd = day_date(d)
    _text(p, R, 96, dd.strftime("%d %b").lstrip("0").upper(), _font(64, QFont.Weight.Bold), INK,
          Qt.AlignmentFlag.AlignRight)
    _text(p, R, 132, dd.strftime("%A").upper() + "  ·  2025", _font(17, QFont.Weight.DemiBold, 2.4), MUTED,
          Qt.AlignmentFlag.AlignRight)
    # --- temperature legend (top right, under the date)
    lw, lh = 250, 10
    lx, ly = R - lw, 170
    g = QLinearGradient(lx, 0, lx + lw, 0)
    for pos, c in RAMP:
        g.setColorAt(pos, QColor(c))
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QBrush(g))
    p.drawRoundedRect(QRectF(lx, ly, lw, lh), 5, 5)
    _text(p, lx, ly + 30, f"{T_MIN:.0f}°C", _font(14), MUTED)
    _text(p, lx + lw, ly + 30, f"{T_MAX:.0f}°C+", _font(14), MUTED, Qt.AlignmentFlag.AlignRight)
    _text(p, lx + lw / 2, ly + 30, "cool   ·   warm", _font(14), MUTED, Qt.AlignmentFlag.AlignHCenter)
    _text(p, R, ly - 8, "WALL COLOUR = WATER TEMPERATURE (EST. FROM AIR TEMP.)", _font(11.5, QFont.Weight.DemiBold, 1.2),
          MUTED, Qt.AlignmentFlag.AlignRight)
    # --- discharge panel
    py = 968
    _text(p, L, py, "DISCHARGE TO THE SEA", _font(16, QFont.Weight.DemiBold, 2.6), MUTED)
    vw = _text(p, L, py + 62, _fmt(tot[d]), _font(58, QFont.Weight.Bold), INK)
    _text(p, L + vw + 12, py + 62, "m³/s", _font(24, QFont.Weight.DemiBold), MUTED)
    _text(p, R, py, "WALL HEIGHT = THAT DAY'S DISCHARGE", _font(11.5, QFont.Weight.DemiBold, 1.2), MUTED,
          Qt.AlignmentFlag.AlignRight)
    fold = tot[hi] / tot[lo]
    _text(p, R, py + 62, f"{fold:.1f}×", _font(46, QFont.Weight.Bold), ACCENT, Qt.AlignmentFlag.AlignRight)
    _text(p, R, py + 86, "LOW-TO-PEAK SWING", _font(12, QFont.Weight.DemiBold, 1.6), MUTED, Qt.AlignmentFlag.AlignRight)
    # sparkline
    sx0, sx1, sy0, sy1 = L, R, py + 125, py + 268
    vmax = tot.max() * 1.08
    X = lambda i: sx0 + (sx1 - sx0) * i / (NDAYS - 1)
    Y = lambda v: sy1 - (sy1 - sy0) * v / vmax
    p.setPen(QPen(QColor(INK.red(), INK.green(), INK.blue(), 40), 1))
    p.drawLine(QPointF(sx0, sy1), QPointF(sx1, sy1))
    ghost = QPainterPath(QPointF(X(0), Y(tot[0])))
    for i in range(1, NDAYS):
        ghost.lineTo(QPointF(X(i), Y(tot[i])))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.setPen(QPen(QColor(INK.red(), INK.green(), INK.blue(), 45), 1.4, Qt.PenStyle.DashLine))
    p.drawPath(ghost)
    area = QPainterPath(QPointF(X(0), sy1))
    for i in range(0, d + 1):
        area.lineTo(QPointF(X(i), Y(tot[i])))
    area.lineTo(QPointF(X(d), sy1))
    area.closeSubpath()
    ga = QLinearGradient(0, sy0, 0, sy1)
    ga.setColorAt(0, QColor(46, 111, 149, 110))
    ga.setColorAt(1, QColor(46, 111, 149, 20))
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QBrush(ga))
    p.drawPath(area)
    line = QPainterPath(QPointF(X(0), Y(tot[0])))
    for i in range(1, d + 1):
        line.lineTo(QPointF(X(i), Y(tot[i])))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.setPen(QPen(ACCENT, 2.4))
    p.drawPath(line)
    # month ticks
    for m in range(1, 13):
        i = (dt.date(2025, m, 1) - DAY0).days
        _text(p, X(i) + 2, sy1 + 18, dt.date(2025, m, 1).strftime("%b")[0], _font(12), MUTED)
    # low & peak markers
    for idx, lab in ((lo, "LOW"), (hi, "PEAK")):
        c = MUTED if idx > d else INK
        p.setPen(QPen(c, 1.2))
        p.setBrush(QBrush(BG))
        p.drawEllipse(QPointF(X(idx), Y(tot[idx])), 4, 4)
        s = f"{lab} {_fmt(tot[idx])} · {day_date(idx).strftime('%d %b').lstrip('0')}"
        al = Qt.AlignmentFlag.AlignLeft if X(idx) < W / 2 else Qt.AlignmentFlag.AlignRight
        off = 8 if al == Qt.AlignmentFlag.AlignLeft else -8
        ty = Y(tot[idx]) - 10 if lab == "PEAK" else Y(tot[idx]) - 20
        _text(p, X(idx) + off, ty, s, _font(13, QFont.Weight.DemiBold), c, al, halo=BG)
    # current-day dot
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QBrush(ACCENT))
    p.drawEllipse(QPointF(X(d), Y(tot[d])), 5, 5)
    # --- data credit
    _text(p, W / 2, H - 22,
          "Data: GloFAS v4 (Copernicus EMS) via Open-Meteo · HydroRIVERS v1.0 · NASA POWER air temp. · geoBoundaries"
          "  |  After Aaron J. Becker", _font(12.5), MUTED, Qt.AlignmentFlag.AlignHCenter)


def render_day(d, path=None):
    RP = builtins.RP
    set_day(d)
    img = QImage(W, H, QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(BG)
    p = QPainter(img)
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    p.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    job = QgsMapRendererCustomPainterJob(map_settings(), p)
    job.renderSynchronously()
    draw_overlay(p, d)
    p.end()
    path = path or os.path.join(FRAMES, f"day_{d:03d}.png")
    img.save(path)
    return path


def render_range(a=0, b=NDAYS - 1, time_budget=None):
    """Render days a..b inclusive. With time_budget (s) it stops early and returns the next day to do."""
    import time
    t0 = time.time()
    os.makedirs(FRAMES, exist_ok=True)
    for d in range(a, b + 1):
        render_day(d)
        if time_budget and time.time() - t0 > time_budget:
            return d + 1
    return b + 1


# ----------------------------------------------------------------- QGIS project (peak day)
def make_project(d=None, path=None):
    RP = builtins.RP
    d = RP["hi"] if d is None else d
    set_day(d)
    gpkg = os.path.join(DATA, "river_pulse_layers.gpkg")
    opts = QgsVectorFileWriter.SaveVectorOptions()
    opts.driverName = "GPKG"
    ctx = QgsProject.instance().transformContext()
    saved = []
    for key, name in (("walls", f"walls_{day_date(d).isoformat()}"), ("land", "canada_land"),
                      ("shadow", "canada_shadow")):
        opts.layerName = name
        opts.actionOnExistingFile = (QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteFile if not saved
                                     else QgsVectorFileWriter.ActionOnExistingFile.CreateOrOverwriteLayer)
        QgsVectorFileWriter.writeAsVectorFormatV3(RP[key], gpkg, ctx, opts)
        saved.append((key, name))
    prj = QgsProject.instance()
    prj.clear()
    prj.setCrs(QgsCoordinateReferenceSystem(MAP_CRS))
    prj.setBackgroundColor(BG)
    prj.setTitle(f"Canada's River Pulse - peak day {day_date(d).isoformat()}")
    for key, name in reversed(saved):
        L = QgsVectorLayer(f"{gpkg}|layername={name}", name, "ogr")
        L.setRenderer(RP[key].renderer().clone())
        prj.addMapLayer(L)
    path = path or os.path.join(ROOT, "Canada_River_Pulse_2025_peak_day.qgz")
    prj.write(path)
    return path


def render_background(a=0, b=NDAYS - 1):
    """Render days a..b in a background thread so QGIS stays responsive. Progress: builtins.RP_RENDER."""
    import threading, time, traceback
    st = builtins.RP_RENDER = {"next": a, "last": b, "done": False, "err": None, "t0": time.time()}

    def work():
        try:
            for d in range(a, b + 1):
                if st.get("stop"):
                    break
                render_day(d)
                st["next"] = d + 1
        except Exception:
            st["err"] = traceback.format_exc()
        st["done"] = True
    threading.Thread(target=work, daemon=True).start()
    return st
