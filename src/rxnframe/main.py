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
import IPython

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
        self.dopts = dopts
        self.drawer = drawer

    def _set_formatter(self):
        # formatters = {
        #     col: self.molecule_to_image_html for col in self.mol_cols
        # }
        formattrs = {}
        formatters = {}
        formatters.update(
            {col: self.reaction_to_image_html for col in self.rxn_cols}
        )
        return formatters

    # image generation functions
    def _mol_to_image(self, mol_smiles, scale=2):
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

    def _rxn_to_image(self, rxn, scale=2):
        """Generate a higher resolution reaction image"""
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

        if self.dopts is None:
            template = new_rxn.GetReactants()[0]
            self._set_draw_options(template)

        # scale up the drawing surface
        drawer = rdMolDraw2D.MolDraw2DCairo(
            self.size[0] * scale, self.size[1] * scale
        )
        drawer.SetDrawOptions(self.dopts)
        drawer.DrawReaction(new_rxn)
        drawer.FinishDrawing()
        img = drawer.GetDrawingText()
        return Image.open(BytesIO(img))
    
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
        img = self._mol_to_image(smiles)
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
        
        img = self._rxn_to_image(rxn)
        img_file = self.make_transparent(self.autocrop(img))
        buffered = BytesIO()
        img_file.save(buffered, format="PNG")
        b64 = base64.b64encode(buffered.getvalue())
        b64 = b64.decode()
        return b64

    def molecule_to_image_html(self, mol, rxn=False):
        return f"""<img {self.options} src="data:image/png;base64,{self._encode_mol(mol)}" alt="Mol"/>"""

    def reaction_to_image_html(self, rxn):
        return f"""<img {self.options} src="data:image/png;base64,{self._encode_rxn(rxn)}" alt="Rxn"/>"""

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
        }}
        </style>
        </head>
        <body>
        {styled_df.to_html()}
        </body>
        </html>
        """

        if pdf_path:
            try:
                from weasyprint import HTML as WeasyHTML
                WeasyHTML(string=html_content, encoding="utf-8").write_pdf(pdf_path)
            except Exception as e:
                print(f"Could not save PDF to {pdf_path}: {e}")

        return HTML(html_content)
