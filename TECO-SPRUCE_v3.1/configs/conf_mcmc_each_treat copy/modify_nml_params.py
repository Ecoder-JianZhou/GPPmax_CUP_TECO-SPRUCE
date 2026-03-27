from pathlib import Path
import re
import copy

# ------------------------------------------------------------
# Input / Output
# ------------------------------------------------------------
input_file = Path("input.nml")          # 改成你的原始 nml 文件名
output_file = Path("modified_params.nml")

text = input_file.read_text()

# ------------------------------------------------------------
# Parsers
# ------------------------------------------------------------
def parse_site_blocks(section_text):
    pattern = re.compile(
        r"(?P<indent>\s*)st\((?P<idx>\d+)\)%parname\s*=\s*'(?P<name>[^']+)'\s*\n"
        r"(?P=indent)st\((?P=idx)\)%parval\s*=\s*(?P<parval>.+?)\s*\n"
        r"(?P=indent)st\((?P=idx)\)%parmin\s*=\s*(?P<parmin>.+?)\s*\n"
        r"(?P=indent)st\((?P=idx)\)%parmax\s*=\s*(?P<parmax>.+?)\s*(?=\n|$)",
        re.S
    )
    blocks = []
    for m in pattern.finditer(section_text):
        d = m.groupdict()
        blocks.append({
            "indent": d["indent"],
            "name": d["name"],
            "parval": d["parval"].strip(),
            "parmin": d["parmin"].strip(),
            "parmax": d["parmax"].strip(),
        })
    return blocks


def parse_sp_blocks(section_text):
    pattern = re.compile(
        r"(?P<indent>\s*)sp\((?P<sp>\d+)\)%var\((?P<idx>\d+)\)%parname\s*=\s*'(?P<name>[^']+)'\s*\n"
        r"(?P=indent)sp\((?P=sp)\)%var\((?P=idx)\)%parval\s*=\s*(?P<parval>.+?)\s*\n"
        r"(?P=indent)sp\((?P=sp)\)%var\((?P=idx)\)%parmin\s*=\s*(?P<parmin>.+?)\s*\n"
        r"(?P=indent)sp\((?P=sp)\)%var\((?P=idx)\)%parmax\s*=\s*(?P<parmax>.+?)\s*(?=\n|$)",
        re.S
    )
    blocks = []
    for m in pattern.finditer(section_text):
        d = m.groupdict()
        blocks.append({
            "indent": d["indent"],
            "sp": int(d["sp"]),
            "name": d["name"],
            "parval": d["parval"].strip(),
            "parmin": d["parmin"].strip(),
            "parmax": d["parmax"].strip(),
        })
    return blocks


# ------------------------------------------------------------
# Block utilities
# ------------------------------------------------------------
def duplicate_insert_after(blocks, src_name, new_name, anchor_name):
    src_block = next((b for b in blocks if b["name"] == src_name), None)
    if src_block is None:
        raise ValueError(f"Source block not found: {src_name}")

    anchor_index = next((i for i, b in enumerate(blocks) if b["name"] == anchor_name), None)
    if anchor_index is None:
        raise ValueError(f"Anchor block not found: {anchor_name}")

    new_block = copy.deepcopy(src_block)
    new_block["name"] = new_name
    blocks.insert(anchor_index + 1, new_block)
    return blocks


def rebuild_site(section_header, blocks):
    lines = [section_header]
    for i, b in enumerate(blocks, start=1):
        indent = b["indent"] or "    "
        lines.extend([
            f"{indent}st({i})%parname = '{b['name']}'",
            f"{indent}st({i})%parval = {b['parval']}",
            f"{indent}st({i})%parmin = {b['parmin']}",
            f"{indent}st({i})%parmax = {b['parmax']}",
        ])
    lines.append("/")
    return "\n".join(lines)


def rebuild_sp(section_header, blocks):
    lines = [section_header]
    species_ids = sorted(set(b["sp"] for b in blocks))
    for spid in species_ids:
        sub = [b for b in blocks if b["sp"] == spid]
        for i, b in enumerate(sub, start=1):
            indent = b["indent"] or "    "
            lines.extend([
                f"{indent}sp({spid})%var({i})%parname = '{b['name']}'",
                f"{indent}sp({spid})%var({i})%parval = {b['parval']}",
                f"{indent}sp({spid})%var({i})%parmin = {b['parmin']}",
                f"{indent}sp({spid})%var({i})%parmax = {b['parmax']}",
            ])
    lines.append("/")
    return "\n".join(lines)


# ------------------------------------------------------------
# Extract sections
# ------------------------------------------------------------
site_match = re.search(r"&sitedaparams\s*(.*?)\n/\s*\n\s*&spdaparams", text, re.S)
sp_match = re.search(r"&spdaparams\s*(.*?)\n/\s*$", text, re.S)

if not site_match:
    raise ValueError("Cannot find &sitedaparams section.")
if not sp_match:
    raise ValueError("Cannot find &spdaparams section.")

site_blocks = parse_site_blocks(site_match.group(1))
sp_blocks = parse_sp_blocks(sp_match.group(1))

# ------------------------------------------------------------
# Modify site params
# ------------------------------------------------------------
site_blocks = duplicate_insert_after(site_blocks, "Tpro_me_1", "Tpro_me_9", "Tpro_me_8")
site_blocks = duplicate_insert_after(site_blocks, "Q10rh_7", "Q10rh_9", "Q10rh_8")
site_blocks = duplicate_insert_after(site_blocks, "Q10pro_7", "Q10pro_9", "Q10pro_8")

# ------------------------------------------------------------
# Modify species params (apply to each species block)
# ------------------------------------------------------------
sp_new = []
for spid in sorted(set(b["sp"] for b in sp_blocks)):
    sub = [copy.deepcopy(b) for b in sp_blocks if b["sp"] == spid]

    sub = duplicate_insert_after(sub, "s_vea_7", "s_vea_9", "s_vea_8")
    sub = duplicate_insert_after(sub, "Entrpy_7", "Entrpy_9", "Entrpy_8")
    sub = duplicate_insert_after(sub, "Vcmax0_7", "Vcmax0_9", "Vcmax0_8")
    sub = duplicate_insert_after(sub, "Q10_7", "Q10_9", "Q10_8")
    sub = duplicate_insert_after(sub, "f_rg_7", "f_rg_9", "f_rg_8")
    sub = duplicate_insert_after(sub, "fn2r_1", "fn2r_9", "fn2r_8")
    sub = duplicate_insert_after(sub, "fn2l_1", "fn2l_9", "fn2l_8")

    sp_new.extend(sub)

# ------------------------------------------------------------
# Rebuild sections
# ------------------------------------------------------------
new_site = rebuild_site("&sitedaparams", site_blocks)
new_sp = rebuild_sp("&spdaparams", sp_new)

new_text = re.sub(
    r"&sitedaparams\s*.*?\n/\s*\n\s*&spdaparams",
    new_site + "\n\n&spdaparams",
    text,
    flags=re.S
)
new_text = re.sub(
    r"&spdaparams\s*.*?\n/\s*$",
    new_sp,
    new_text,
    flags=re.S
)

output_file.write_text(new_text)
print(f"Done. Output written to: {output_file}")
