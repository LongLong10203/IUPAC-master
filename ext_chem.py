import random
from rdkit import Chem
from rdkit.Chem import Draw
from pubchempy import get_compounds, PubChemPyError

import io
import base64
from PIL import Image

def random_smiles() -> str:
    smiles = "C"

    for _ in range(random.randint(1, 5)):
        rnd = random.randint(1, 4)

        # functional group
        if rnd == 1:
            rnd2 = random.randint(1, 2)
            if rnd2 == 1: # -ene
                smiles += "C=C"
            elif rnd2 == 2: # -ol
                smiles += "C(O)"
        
        # halogen
        elif rnd == 2:
            rnd2 = random.randint(1, 5)
            if rnd2 == 1:
                smiles += "(Br)"
            elif rnd2 == 2:
                smiles += "(Cl)"
            elif rnd2 == 3:
                smiles += "(F)"
            elif rnd2 == 4:
                smiles += "(C)"
            else:
                smiles += "(CC)"

        else:
            smiles += "C"

    # hydroxyl group should not collapse with carboxyl group (out-syl)
    if "C(O)" not in smiles:
        # 1/4 probability -oic group in front
        if random.randint(1, 4) == 1:
            smiles = "OC(=O)" + smiles
        
        # 1/4 probability -oic group at back
        if random.randint(1, 4) == 1:
            smiles += "C(=O)O"

    # i have no idea how SMILES work
    if valid(smiles):
        return smiles
    else:
        return random_smiles()

def iupac_name(smiles: str) -> str:
    try:
        compounds = get_compounds(smiles, namespace="smiles")
        if compounds[0]:
            return compounds[0].iupac_name.replace("(", "").replace(")", "") # to prevent something like dibromo(chloro)methane
        else:
            return None
    except PubChemPyError as e:
        print(e)
        return None

def valid(smiles: str):
    return Chem.MolFromSmiles(smiles) is not None and iupac_name(smiles) is not None

def generate_image(smiles: str) -> Image.Image:
    return Draw.MolToImage(Chem.MolFromSmiles(smiles))

def image_to_base64(image: Image.Image) -> str:
    img_byte_arr = io.BytesIO()
    image.save(img_byte_arr, format="PNG")
    img_byte_arr = img_byte_arr.getvalue()
    img_base64 = base64.b64encode(img_byte_arr).decode("utf-8")
    return img_base64

def generate_base64_image(smiles: str) -> str:
    return image_to_base64(generate_image(smiles))

class OrganicCompound:
    def __init__(self):
        self.smiles = random_smiles()
        self.iupac = iupac_name(self.smiles)
        self.img = generate_image(self.smiles)
        self.img_base64 = image_to_base64(self.img)

# DEBUG
if __name__ == "__main__":
    # compound = OrganicCompound()
    # print(compound.smiles)
    # print(compound.iupac)
    # compound.img.show()
    iupac_name("idk")