import random
from rdkit import Chem
from rdkit.Chem import Draw
from pubchempy import get_compounds, PubChemPyError

import io
import base64
from PIL import Image

def random_smiles_easy() -> str:
    # meth ~ oct, normal distribution
    longest_carbon_chain = random.choices([i for i in range(1, 9)], weights=[min([i, 8 - i, 3]) for i in range(1, 9)])[0]
    smiles = ["C"] * longest_carbon_chain

    cyclic = False
    # 1/4 probability cyclic
    if longest_carbon_chain >= 3 and random.randint(1, 4) == 1:
        cyclic = True
        smiles.insert(1, "1")
        smiles.append("1")

    # alkene (at least 2 C atoms) 1/3 probability to include
    if longest_carbon_chain >= 2 and random.randint(1, 3) == 1:
        no_of_alkenes = random.randint(1, longest_carbon_chain-1)
        for _ in range(no_of_alkenes):
            pos = random.randint(0, len(smiles)-2)
            # check invalid, also lower the probability of more alkenes
            if smiles[pos] == "=" or smiles[pos+1] == "=":
                continue
            # extract the selected atoms
            first = smiles.pop(pos)
            second = smiles.pop(pos)
            # insert alkene group at selected place
            smiles.insert(pos, f"{first}={second}")

    # 0 ~ 4, normal distribution
    no_of_substituents = random.choices([i for i in range(0, 5)], weights=[min(i, 4 - i) for i in range(0, 5)])[0]
    substituents = [
        "Br", "Cl", "F", "C", "CC"
    ]
    for _ in range(no_of_substituents):
        substituent = random.choice(substituents)
        smiles.insert(random.randint(0 if 1 > longest_carbon_chain-1 else 1, longest_carbon_chain-1), f"({substituent})")

    rnd = random.randint(1, 4)
    # either hydroxyl group or carboxyl group
    if rnd == 1:
        smiles.insert(random.randint(0, len(smiles)-1), "(O)")
    elif rnd == 2 and not cyclic:
        if random.randint(1, 2) == 1: # at front
            smiles.insert(0, "OC(=O)")
        else: # at back
            smiles.append("C(=O)O")

    smiles = "".join(smiles)

    # i have no idea how SMILES work
    if valid(smiles):
        return smiles
    else:
        return random_smiles_easy()

# TODO
def random_smiles_hard() -> str:
    ...

def iupac_name(smiles: str) -> str:
    try:
        compounds = get_compounds(smiles, namespace="smiles")
        if compounds[0].iupac_name is not None:
            return compounds[0].iupac_name.replace("(", "").replace(")", "") # to prevent something like dibromo(chloro)methane
        else:
            # print(f"Failed: {smiles}")
            return None
    except PubChemPyError as e:
        # print(e)
        return None

def valid(smiles: str):
    return iupac_name(smiles) is not None and Chem.MolFromSmiles(smiles) is not None

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

# DEBUG
if __name__ == "__main__":
    ...

    print(iupac_name("C(=O)C"))