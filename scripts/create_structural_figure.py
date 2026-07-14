#!/usr/bin/env python3
"""Render the manuscript structural panels from archived receptor and pose files."""
from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from pymol import cmd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "manuscript"
OUT.mkdir(parents=True, exist_ok=True)
PANEL_DIR = OUT / "structural_panels"
PANEL_DIR.mkdir(parents=True, exist_ok=True)

KEY_RESIDUES = "resi 745+790+793+855"
WIDTH, HEIGHT = 1500, 1050


def reset_scene() -> None:
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("cartoon_transparency", 0.12)
    cmd.set("stick_radius", 0.18)
    cmd.set("sphere_scale", 0.22)


def protein_style(object_name: str, ligand_selection: str) -> None:
    cmd.hide("everything", object_name)
    cmd.show("cartoon", object_name)
    cmd.color("gray80", object_name)
    pocket = f"({object_name} and byres ({object_name} within 5 of ({ligand_selection})))"
    cmd.show("sticks", pocket)
    cmd.color("gray60", pocket)
    cmd.show("sticks", f"{object_name} and ({KEY_RESIDUES})")
    cmd.color("tv_blue", f"{object_name} and ({KEY_RESIDUES})")


def render(path: Path, selection: str) -> None:
    cmd.orient(selection)
    cmd.zoom(selection, 7.0)
    cmd.turn("x", -12)
    cmd.turn("y", 8)
    cmd.ray(WIDTH, HEIGHT)
    cmd.png(str(path), dpi=300)


def panel_native() -> Path:
    reset_scene()
    receptor = ROOT / "data/external/structural_examples/1xkk_a_fmm_91/receptor_clean.pdb"
    ligand = ROOT / "data/external/structural_examples/1xkk_a_fmm_91/native_ligand.pdb"
    cmd.load(str(receptor), "protein")
    cmd.load(str(ligand), "native")
    protein_style("protein", "native")
    cmd.show("sticks", "native")
    cmd.color("orange", "native and elem C")
    cmd.color("blue", "native and elem N")
    cmd.color("red", "native and elem O")
    cmd.label("protein and name CA and resi 745+790+793+855", 'resn+resi')
    cmd.set("label_size", 18)
    cmd.set("label_color", "black")
    path = PANEL_DIR / "panel_a_native.png"
    render(path, "native or (protein and resi 745+790+793+855)")
    return path


def panel_redocking() -> Path:
    reset_scene()
    root = ROOT / "data/external/structural_examples/redocking_4hjo"
    receptor = root / "receptor.pdb"
    reference = root / "reference.pdbqt"
    rank1 = root / "rank1_pose.pdbqt"
    native_like = root / "native_like_pose.pdbqt"
    cmd.load(str(receptor), "protein")
    cmd.load(str(reference), "reference")
    cmd.load(str(rank1), "rank1")
    cmd.load(str(native_like), "native_like")
    protein_style("protein", "reference")
    for name, color in [("reference", "gray35"), ("rank1", "red"), ("native_like", "forest")]:
        cmd.show("sticks", name)
        cmd.color(color, f"{name} and elem C")
        cmd.color("blue", f"{name} and elem N")
        cmd.color("red", f"{name} and elem O")
    path = PANEL_DIR / "panel_b_redocking.png"
    render(path, "reference or rank1 or native_like")
    return path


def panel_enrichment(
    receptor: Path,
    ligand: Path,
    output_name: str,
    ligand_color: str,
) -> Path:
    reset_scene()
    cmd.load(str(receptor), "protein")
    cmd.load(str(ligand), "pose")
    protein_style("protein", "pose")
    cmd.show("sticks", "pose")
    cmd.color(ligand_color, "pose and elem C")
    cmd.color("blue", "pose and elem N")
    cmd.color("red", "pose and elem O")
    cmd.show("surface", "protein within 4.5 of pose")
    cmd.set("transparency", 0.62, "protein within 4.5 of pose")
    path = PANEL_DIR / output_name
    render(path, "pose or (protein within 5 of pose)")
    return path


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    return ImageFont.truetype(name, size)


def assemble(paths: list[Path]) -> None:
    images = []
    for path in paths:
        source = Image.open(path).convert("RGBA")
        background = Image.new("RGBA", source.size, "white")
        images.append(Image.alpha_composite(background, source).convert("RGB"))
    tile_w, tile_h = 1450, 1015
    margin, header = 55, 120
    canvas = Image.new("RGB", (2 * tile_w + 3 * margin, 2 * (tile_h + header) + 3 * margin), "white")
    titles = [
        ("A", "Native-derived interaction prior", "1XKK/FMM; key EGFR residues highlighted"),
        ("B", "Sampling and ranking are distinct", "AQ4: rank 1, 5.81 Å (red); rank 11, 0.83 Å (green); crystal pose gray"),
        ("C", "High neural score, weak prior recovery", "DUD-E decoy C49121561; CNNscore 0.951, recall 0.158; rank 57 to 820"),
        ("D", "Interaction-consistent pose promoted", "DUD-E active 475036; CNNscore 0.895, recall 0.421; rank 519 to 27"),
    ]
    draw = ImageDraw.Draw(canvas)
    for index, (image, (letter, title, subtitle)) in enumerate(zip(images, titles)):
        row, col = divmod(index, 2)
        x = margin + col * (tile_w + margin)
        y = margin + row * (tile_h + header + margin)
        image.thumbnail((tile_w, tile_h), Image.Resampling.LANCZOS)
        px = x + (tile_w - image.width) // 2
        py = y + header + (tile_h - image.height) // 2
        canvas.paste(image, (px, py))
        draw.text((x, y), letter, fill="#172126", font=font(44, True))
        draw.text((x + 60, y + 2), title, fill="#172126", font=font(34, True))
        draw.text((x + 60, y + 54), subtitle, fill="#41535c", font=font(24))
    canvas.save(OUT / "figure1_structural.png", dpi=(300, 300))


def main() -> None:
    paths = [
        panel_native(),
        panel_redocking(),
        panel_enrichment(
            ROOT / "data/external/structural_examples/4hjo_a_aq4_1001/receptor.h.pdb",
            ROOT / "data/external/structural_examples/4hjo_a_aq4_1001/C49121561_top.posed_h.pdb",
            "panel_c_low_recall.png",
            "red",
        ),
        panel_enrichment(
            ROOT / "data/external/structural_examples/1xkk_a_fmm_91/receptor.h.pdb",
            ROOT / "data/external/structural_examples/1xkk_a_fmm_91/475036_top.posed_h.pdb",
            "panel_d_promoted.png",
            "forest",
        ),
    ]
    assemble(paths)


if __name__ == "__main__":
    main()
