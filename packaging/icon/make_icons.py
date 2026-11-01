"""Build the Dougalien Systems letter marks (D, R, M): a letter cut from layered rock with a small fault.
The layers and the fault are copied exactly from dougalien_mark.svg; only the letter outline changes."""
import io, sys, cairosvg
from PIL import Image

LAYERS = ('<rect x="0" y="0" width="62" height="30" fill="#D6B27C"/><rect x="0" y="30" width="62" height="22" fill="#3E5C6B"/>'
 '<rect x="0" y="52" width="62" height="26" fill="#B4532A"/><rect x="0" y="78" width="62" height="18" fill="#D6B27C"/>'
 '<rect x="0" y="96" width="62" height="24" fill="#3E5C6B"/><rect x="65" y="-12" width="55" height="30" fill="#D6B27C"/>'
 '<rect x="65" y="18" width="55" height="22" fill="#3E5C6B"/><rect x="65" y="40" width="55" height="26" fill="#B4532A"/>'
 '<rect x="65" y="66" width="55" height="18" fill="#D6B27C"/><rect x="65" y="84" width="55" height="24" fill="#3E5C6B"/>'
 '<rect x="65" y="108" width="55" height="12" fill="#D6B27C"/>')

LETTERS = {
    "D": "M22 10 H58 A50 50 0 0 1 58 110 H22 Z",
    "R": "M22 10 H74 A29 29 0 0 1 90.3 63 L110 110 H82 L64.4 68 H50 V110 H22 Z M50 31 V50 H73 A9.5 9.5 0 0 0 73 31 Z",
    "M": "M22 110 V10 H48 L63.5 50 L79 10 H105 V110 H81 V56 L68 88 H59 L46 56 V110 Z",
}
LABEL = {"D": "Dougalien Systems", "R": "Reaction", "M": "Magma"}

def svg(letter):
    return (f'<svg role="img" aria-label="{LABEL[letter]}" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 120 120">'
            f'<clipPath id="m"><path clip-rule="evenodd" d="{LETTERS[letter]}"/></clipPath>'
            f'<g clip-path="url(#m)">{LAYERS}</g></svg>\n')

def png(letter, size):
    return Image.open(io.BytesIO(cairosvg.svg2png(bytestring=svg(letter).encode(), output_width=size, output_height=size))).convert("RGBA")

def ico(letter, path):
    sizes = [16, 24, 32, 48, 64, 128, 256]
    big = png(letter, 256)
    big.save(path, format="ICO", sizes=[(s, s) for s in sizes],
             append_images=[png(letter, s) for s in sizes[:-1]])

if __name__ == "__main__":
    for L, name in (("R", "reaction"), ("M", "magma")):
        open(f"{name}_mark.svg", "w").write(svg(L))
        ico(L, f"{name}.ico")
        png(L, 512).save(f"{name}_512.png")
    # preview sheet: D, R, M at 256 and at 32 and 16 (scaled up so small sizes can be judged)
    sheet = Image.new("RGBA", (3 * 300, 300 + 140), "white")
    for i, L in enumerate("DRM"):
        sheet.alpha_composite(png(L, 256), (i * 300 + 22, 20))
        sheet.alpha_composite(png(L, 32).resize((96, 96), Image.NEAREST), (i * 300 + 40, 310))
        sheet.alpha_composite(png(L, 16).resize((96, 96), Image.NEAREST), (i * 300 + 160, 310))
    sheet.save("preview.png")
