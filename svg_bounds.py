r"""svg_bounds.py — 调查 results/ 下所有 SVG 的内容边界 vs viewBox。

零依赖：xml.etree 解析，估算 text 宽度（font-size * 0.55 * 字符数，
text-anchor=start/end/middle 处理方向）。
用法：python svg_bounds.py
"""
import os
import re
import xml.etree.ElementTree as ET

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
SVGNS = "{http://www.w3.org/2000/svg}"


def est_text_right(x, anchor, fs, s):
    w = fs * 0.55 * len(s)
    if anchor == "middle":
        return x + w / 2
    if anchor == "end":
        return x
    return x + w


for fn in sorted(os.listdir(RES)):
    if not fn.endswith(".svg"):
        continue
    path = os.path.join(RES, fn)
    tree = ET.parse(path)
    root = tree.getroot()
    vb = root.get("viewBox", "")
    vbw = float(vb.split()[2]) if vb else 0.0
    vbh = float(vb.split()[3]) if vb else 0.0
    maxr, offender = 0.0, ""
    for el in root.iter():
        tag = el.tag.replace(SVGNS, "")
        if tag == "text" and el.text:
            x = float(el.get("x", 0))
            fs = float(el.get("font-size", 12))
            anchor = el.get("text-anchor", "start")
            r = est_text_right(x, anchor, fs, el.text)
            if r > maxr:
                maxr, offender = r, "text x=%g fs=%g %r" % (x, fs, el.text[:30])
        elif tag == "rect":
            r = float(el.get("x", 0)) + float(el.get("width", 0))
            if r > maxr:
                maxr, offender = r, "rect x=%s w=%s" % (el.get("x", "?"), el.get("width", "?"))
        elif tag == "line":
            r = max(float(el.get("x1", 0)), float(el.get("x2", 0)))
            if r > maxr:
                maxr, offender = r, "line x2=%g" % el.get("x2")
    flag = "  <-- OVERFLOW (need %.0f)" % (maxr + 8) if maxr > vbw else ""
    print("%-24s viewBox=%gx%g  content-right=%.0f  %s%s"
          % (fn, vbw, vbh, maxr, offender[:40], flag))
