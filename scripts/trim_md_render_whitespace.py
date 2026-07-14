#!/usr/bin/env python3
"""Trim transparent canvas margins from MD renders without trimming molecular content."""

from pathlib import Path

from PIL import Image


RENDER_DIR = Path("figures/manuscript")
FILENAMES = (
    "md_control_002_rep01_start.png",
    "md_control_002_rep01_end.png",
    "md_misdocked_control_rep01_start.png",
    "md_misdocked_control_rep01_end.png",
)
PADDING_PX = 80


def main() -> None:
    for filename in FILENAMES:
        path = RENDER_DIR / filename
        image = Image.open(path).convert("RGBA")
        alpha = image.getchannel("A")
        bounds = alpha.getbbox()
        if bounds is None:
            raise RuntimeError(f"No rendered content found in {path}")
        left, top, right, bottom = bounds
        left = max(0, left - PADDING_PX)
        top = max(0, top - PADDING_PX)
        right = min(image.width, right + PADDING_PX)
        bottom = min(image.height, bottom + PADDING_PX)
        image.crop((left, top, right, bottom)).save(path)


if __name__ == "__main__":
    main()
