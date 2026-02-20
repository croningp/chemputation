"""Code to display reactions in a DataFrame with images of molecules
   
   Authors: Sebastian Pagel
   

Adapted from:
https://github.com/pagel-s/rxn_frame    
    
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D, rdDepictor
from rdkit.Chem import Draw
from rdkit.Chem import rdCoordGen
import pandas as pd
from PIL import Image, ImageChops
from io import BytesIO
import base64
import re

from IPython.core.display import HTML

class RXN_frame:
    def __init__(
        self,
        df,
        mol_cols,
        rxn_cols,
        filter_criteria=None,
        size=(200, 100),
        img_only=False,
        rxn_render_scale=4,
        mol_render_scale=4,
        prefer_svg=True,
        **kwargs,
    ):
        self.df = df
        self.mol_cols = mol_cols
        self.rxn_cols = rxn_cols
        self.filter_criteria = filter_criteria
        self.size = size
        self.dopts = None
        self.drawer = None
        self.options = {}
        self.img_only = img_only
        self.rxn_render_scale = max(1, int(rxn_render_scale))
        self.mol_render_scale = max(1, int(mol_render_scale))
        self.prefer_svg = prefer_svg
        self.kwargs = kwargs

        if filter_criteria is not None:
            self.df.query(filter_criteria, inplace=True)

        print(f"Number of rows: {self.df.shape[0]}")
        self.df = self.df.drop(
            self.df.filter(regex="^Unnamed").columns, axis=1
        )

    def _set_draw_options(self, template):
        if isinstance(template, str):
            template = Chem.MolFromSmiles(template)
        rdDepictor.Compute2DCoords(template)
        drawer = rdMolDraw2D.MolDraw2DCairo(-1, -1, self.size[0], self.size[1])
        dopts = drawer.drawOptions()
        dopts.bondLineWidth = 1.0  # default is 2.0
        Draw.SetACS1996Mode(dopts, Draw.MeanBondLength(template) * 0.70)
        dopts.bondLineWidth = 1.5
        dopts.baseFontSize = 12

        # Reduce excess whitespace around drawings and keep reactions readable.
        # Some RDKit builds expose additional padding knobs; set them when present.
        if hasattr(dopts, "padding"):
            dopts.padding = 0.01
        if hasattr(dopts, "reagentPadding"):
            dopts.reagentPadding = 0.0
        if hasattr(dopts, "componentPadding"):
            dopts.componentPadding = 0.0
        if hasattr(dopts, "centreMoleculesBeforeDrawing"):
            dopts.centreMoleculesBeforeDrawing = True

        self.dopts = dopts
        self.drawer = drawer

    def _set_formatter(self):
        # formatters = {
        #     col: self.molecule_to_image_html for col in self.mol_cols
        # }
        formatters = {}
        formatters.update(
            {col: self.reaction_to_image_html for col in self.rxn_cols}
        )
        return formatters

    # image generation functions
    def _mol_to_image(self, mol_smiles, scale=4):
        """Generate a higher resolution molecule image"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if self.dopts is None:
            self._set_draw_options(mol)

        # render larger (scale multiplier) for better quality
        big_size = (self.size[0] * scale, self.size[1] * scale)

        return Draw.MolToImage(
            mol,
            size=big_size,
            kekulize=True,
            wedgeBonds=True,
            imageType="png",
            options=self.dopts,
        )

    def _rxn_to_image(self, rxn, scale=1):
        """Generate a higher resolution reaction image"""
        rxn = self._prepare_reaction(rxn)

        if self.dopts is None:
            template = rxn.GetReactants()[0]
            self._set_draw_options(template)

        # scale up the drawing surface
        drawer = rdMolDraw2D.MolDraw2DCairo(
            self.size[0] * scale, self.size[1] * scale
        )
        drawer.SetDrawOptions(self.dopts)
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        img = drawer.GetDrawingText()
        return Image.open(BytesIO(img))

    def _prepare_reaction(self, rxn):
        rxn = AllChem.ReactionFromSmarts(rxn, useSmiles=True)
        new_rxn = AllChem.ChemicalReaction()
        for mol in rxn.GetReactants():
            mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
            rdCoordGen.AddCoords(mol)
            new_rxn.AddReactantTemplate(mol)
        for mol in rxn.GetProducts():
            rdCoordGen.AddCoords(mol)
            mol = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol))
            new_rxn.AddProductTemplate(mol)

        return new_rxn

    def _rxn_to_svg(self, rxn):
        """Generate a vector reaction drawing (SVG)."""
        rxn = self._prepare_reaction(rxn)

        if self.dopts is None:
            template = rxn.GetReactants()[0]
            self._set_draw_options(template)

        drawer = rdMolDraw2D.MolDraw2DSVG(
            self.size[0], self.size[1]
        )
        drawer.SetDrawOptions(self.dopts)
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    
    # image processing functions
    def make_transparent(self, img):
        img = img.convert("RGBA")
        pixdata = img.load()
        width, height = img.size
        for y in range(height):
            for x in range(width):
                if pixdata[x, y] == (255, 255, 255, 255):
                    pixdata[x, y] = (255, 255, 255, 0)
        return img

    def autocrop(self, img, bgcolor="white"):
        if img.mode != "RGB":
            img = img.convert("RGB")
        bg = Image.new("RGB", img.size, bgcolor)
        diff = ImageChops.difference(img, bg)
        bbox = diff.getbbox()
        if bbox:
            return img.crop(bbox)
        return None  # no contents

    def _encode_mol(self, smiles):
        img = self._mol_to_image(smiles, scale=self.mol_render_scale)
        img_file = self.autocrop(img) #self.make_transparent(self.autocrop(img))
        buffered = BytesIO()
        img_file.save(buffered, format="PNG")
        b64 = base64.b64encode(buffered.getvalue())
        b64 = b64.decode()
        return b64

    def _encode_rxn(self, rxn):
        transparent_pixel = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNkYAAAAAYAAjCB0C8AAAAASUVORK5CYII="

        if not rxn or pd.isna(rxn) or (isinstance(rxn, str) and not rxn.strip()):
            return transparent_pixel
        
        img = self._rxn_to_image(rxn, scale=self.rxn_render_scale)
        img_file = img#self.make_transparent(self.autocrop(img))
        buffered = BytesIO()
        img_file.save(buffered, format="PNG")
        b64 = base64.b64encode(buffered.getvalue())
        b64 = b64.decode()
        return b64

    def molecule_to_image_html(self, mol, rxn=False):
        return f"""<img {self.options} src="data:image/png;base64,{self._encode_mol(mol)}" alt="Mol"/>"""

    def _reaction_to_svg_html(self, rxn):
        if not rxn or pd.isna(rxn) or (isinstance(rxn, str) and not rxn.strip()):
            return f'<div class="rxn-svg" style="width:{self.size[0]}px;height:{self.size[1]}px;"></div>'

        svg = self._rxn_to_svg(rxn)
        if svg.startswith("<?xml"):
            svg = svg.split("?>", 1)[1]

        svg = self._tighten_svg_viewbox(svg, x_frac=0.14, y_frac=0.09)

        svg = svg.replace(
            "<svg ",
            '<svg width="100%" height="100%" preserveAspectRatio="xMidYMid meet" ',
            1,
        )
        return f'<div class="rxn-svg" style="width:{self.size[0]}px;height:{self.size[1]}px;">{svg}</div>'

    def _tighten_svg_viewbox(self, svg, x_frac=0.12, y_frac=0.08):
        """Apply symmetric center-crop on SVG viewBox to reduce border whitespace."""
        viewbox_match = re.search(
            r"\bviewBox\s*=\s*['\"]\s*([\-\d\.eE]+)\s+([\-\d\.eE]+)\s+([\-\d\.eE]+)\s+([\-\d\.eE]+)\s*['\"]",
            svg,
        )

        if viewbox_match:
            x0 = float(viewbox_match.group(1))
            y0 = float(viewbox_match.group(2))
            width = float(viewbox_match.group(3))
            height = float(viewbox_match.group(4))
        else:
            width_match = re.search(r"\bwidth\s*=\s*['\"]([\-\d\.eE]+)", svg)
            height_match = re.search(r"\bheight\s*=\s*['\"]([\-\d\.eE]+)", svg)
            if not width_match or not height_match:
                return svg
            x0, y0 = 0.0, 0.0
            width = float(width_match.group(1))
            height = float(height_match.group(1))

        if width <= 0 or height <= 0:
            return svg

        x_pad = width * max(0.0, min(0.25, x_frac))
        y_pad = height * max(0.0, min(0.25, y_frac))

        new_x = x0 + x_pad
        new_y = y0 + y_pad
        new_w = width - (2 * x_pad)
        new_h = height - (2 * y_pad)

        if new_w <= 0 or new_h <= 0:
            return svg

        new_viewbox = f'{new_x:.3f} {new_y:.3f} {new_w:.3f} {new_h:.3f}'

        if viewbox_match:
            return re.sub(
                r"\bviewBox\s*=\s*['\"][^'\"]+['\"]",
                f'viewBox="{new_viewbox}"',
                svg,
                count=1,
            )

        return re.sub(r"<svg\b", f'<svg viewBox="{new_viewbox}"', svg, count=1)

    def reaction_to_image_html(self, rxn):
        if self.prefer_svg:
            return self._reaction_to_svg_html(rxn)
        return f"""<img {self.options} src="data:image/png;base64,{self._encode_rxn(rxn)}" alt="Rxn"/>"""

    def _save_pdf(self, html_content, pdf_path):
        """Save HTML as PDF via Edge/Chrome headless (no native libs needed).

        Strategy (in order):
          1. Microsoft Edge headless  – always present on Windows 10/11
          2. Google Chrome headless   – if installed
          3. pdfkit                   – if installed (needs wkhtmltopdf)
          4. Save as .html fallback   – open in browser → File → Print → Save as PDF
        """
        import os
        import subprocess
        import tempfile
        from pathlib import Path

        pdf_path = str(Path(pdf_path).resolve())

        def _is_valid_pdf(path):
            if not os.path.exists(path):
                return False
            if os.path.getsize(path) < 1024:
                return False
            try:
                with open(path, "rb") as f:
                    return f.read(5) == b"%PDF-"
            except OSError:
                return False

        # Write HTML to a temp file so the browser can load it as a file:// URL
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".html", encoding="utf-8", delete=False
        ) as tmp:
            tmp.write(html_content)
            tmp_html = tmp.name

        try:
            # Remove stale/old file so success checks are always for fresh output.
            if os.path.exists(pdf_path):
                os.unlink(pdf_path)

            # ── 1 & 2: Edge / Chrome headless ────────────────────────────────
            candidates = [
                r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
                r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
                r"C:\Program Files\Google\Chrome\Application\chrome.exe",
                r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
            ]
            html_url = Path(tmp_html).as_uri()

            arg_sets = [
                ["--headless=new"],
                ["--headless"],
            ]

            for browser in candidates:
                if not os.path.exists(browser):
                    continue
                for headless_args in arg_sets:
                    cmd = [
                        browser,
                        *headless_args,
                        "--disable-gpu",
                        "--allow-file-access-from-files",
                        "--run-all-compositor-stages-before-draw",
                        "--virtual-time-budget=5000",
                        f"--print-to-pdf={pdf_path}",
                        "--print-to-pdf-no-header",
                        "--no-margins",
                        html_url,
                    ]
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        timeout=90,
                    )
                    if result.returncode == 0 and _is_valid_pdf(pdf_path):
                        print(f"PDF saved to {pdf_path}")
                        return

            # ── 3: pdfkit ────────────────────────────────────────────────────
            try:
                import pdfkit
                pdfkit.from_file(tmp_html, pdf_path)
                if _is_valid_pdf(pdf_path):
                    print(f"PDF saved to {pdf_path}")
                    return
            except ImportError:
                pass

            # ── 4: HTML fallback ─────────────────────────────────────────────
            html_out = pdf_path.replace(".pdf", ".html")
            with open(html_out, "w", encoding="utf-8") as f:
                f.write(html_content)
            print(
                f"Could not generate PDF automatically.\n"
                f"HTML saved to: {html_out}\n"
                "Open it in a browser and use File → Print → Save as PDF."
            )
        except Exception as e:
            print(f"Could not save PDF: {e}")
        finally:
            try:
                os.unlink(tmp_html)
            except OSError:
                pass

    def display_rxns(self, pdf_path=None):
        self.options = "style='display:inline; margin:5px;'"
        formatters = self._set_formatter()

        display_df = self.df.copy()
        display_df.index = range(1, len(display_df) + 1)

        styled_df = (
            display_df[self.rxn_cols + self.mol_cols].style.format(formatters)
            if self.img_only
            else display_df.style.format(formatters)
        )
        

        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
        <meta charset="utf-8">
        <style>
        @page {{
            size: A4 landscape;
            margin: 20mm;
        }}
        table, th, td {{
            border: 1px solid black;
            border-collapse: collapse;
            padding: 6px;
        }}
        th {{
            background-color: #f2f2f2;
            text-align: center;
        }}
        td {{
            vertical-align: middle;
        }}
        table {{
            width: 100%;
            table-layout: auto;
        }}
        img {{
            max-width: {self.size[0]}px;
            max-height: {self.size[1]}px;
            width: auto;
            height: auto;
        }}
        .rxn-svg {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            vertical-align: middle;
            overflow: hidden;
        }}
        .rxn-svg svg {{
            width: 100%;
            height: 100%;
            display: block;
            margin: 0 auto;
        }}
        </style>
        </head>
        <body>
        {styled_df.to_html()}
        </body>
        </html>
        """

        if pdf_path:
            self._save_pdf(html_content, pdf_path)

        return HTML(html_content)

    def display_rxns_grouped(self, group_col, rowspan_cols=None, pdf_path=None):
        """Display reactions with one step per row, grouped by group_col.

        Uses HTML rowspan to visually merge cells that share the same group,
        making it clear which steps belong to the same multi-step reaction.

        Args:
            group_col:    Column whose consecutive equal values define a group.
            rowspan_cols: Columns to merge via rowspan within each group.
                          Defaults to [group_col].
            pdf_path:     Optional path to save the result as a PDF.
        """
        self.options = "style='display:inline; margin:5px;'"

        if rowspan_cols is None:
            rowspan_cols = [group_col]

        df = self.df.copy().reset_index(drop=True)
        all_cols = list(df.columns)
        rxn_mol_cols = self.rxn_cols + self.mol_cols
        non_rxn_cols = [c for c in all_cols if c not in rxn_mol_cols]
        header_cols = non_rxn_cols + rxn_mol_cols

        # ------------------------------------------------------------------
        # Find consecutive group boundaries
        # ------------------------------------------------------------------
        group_boundaries = []
        i = 0
        while i < len(df):
            gval = df.at[i, group_col]
            j = i + 1
            while j < len(df) and df.at[j, group_col] == gval:
                j += 1
            group_boundaries.append((i, j - 1))
            i = j

        # Per-row metadata: alternating background colour + last-in-group flag
        row_meta = {}
        for g_idx, (start, end) in enumerate(group_boundaries):
            color = "#ffffff" if g_idx % 2 == 0 else "#eef2f8"
            for k in range(start, end + 1):
                row_meta[k] = {"color": color, "is_last": k == end}

        # Rowspan map: (row_idx, col) -> span  |  0 = cell is covered, skip
        rowspan_map = {}
        for start, end in group_boundaries:
            span = end - start + 1
            for col in rowspan_cols:
                rowspan_map[(start, col)] = span
                for k in range(start + 1, end + 1):
                    rowspan_map[(k, col)] = 0

        # ------------------------------------------------------------------
        # Build HTML table
        # ------------------------------------------------------------------
        headers = "".join(f"<th>{col}</th>" for col in header_cols)

        GROUP_SEP = "border-bottom: 2px solid #555;"

        rows_html = []
        for i, row in df.iterrows():
            meta = row_meta[i]
            bg = meta["color"]
            is_last = meta["is_last"]
            row_sep = GROUP_SEP if is_last else ""

            cells = []

            for col in non_rxn_cols:
                if col in rowspan_cols:
                    span = rowspan_map.get((i, col), 1)
                    if span == 0:
                        continue  # covered by a rowspan from above
                    val = "" if pd.isnull(row[col]) else row[col]
                    # rowspan cells always bottom-align with the group border
                    cells.append(
                        f'<td rowspan="{span}" style="background:{bg};'
                        f" vertical-align:middle; text-align:center;"
                        f' {GROUP_SEP}">{val}</td>'
                    )
                else:
                    val = "" if pd.isnull(row[col]) else row[col]
                    cells.append(
                        f'<td style="background:{bg}; vertical-align:middle;'
                        f' text-align:center; {row_sep}">{val}</td>'
                    )

            for col in rxn_mol_cols:
                if col in self.rxn_cols:
                    img_html = self.reaction_to_image_html(row[col])
                else:
                    img_html = self.molecule_to_image_html(row[col])
                cells.append(
                    f'<td style="background:{bg}; vertical-align:middle;'
                    f' text-align:center; {row_sep}">{img_html}</td>'
                )

            rows_html.append(f"<tr>{''.join(cells)}</tr>")

        rows_body = "".join(rows_html)
        table_html = (
            f"<table>"
            f"<thead><tr>{headers}</tr></thead>"
            f"<tbody>{rows_body}</tbody>"
            f"</table>"
        )

        html_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<style>
@page {{
    size: A4 landscape;
    margin: 20mm;
}}
table, th, td {{
    border: 1px solid #ccc;
    border-collapse: collapse;
    padding: 8px;
}}
th {{
    background-color: #d0d8e8;
    text-align: center;
    font-weight: bold;
}}
table {{
    width: 100%;
    table-layout: auto;
}}
img {{
    max-width: {self.size[0]}px;
    max-height: {self.size[1]}px;
    width: auto;
    height: auto;
}}
.rxn-svg {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    vertical-align: middle;
    overflow: hidden;
}}
.rxn-svg svg {{
    width: 100%;
    height: 100%;
    display: block;
    margin: 0 auto;
}}
</style>
</head>
<body>
{table_html}
</body>
</html>"""

        if pdf_path:
            self._save_pdf(html_content, pdf_path)
        return HTML(html_content)
