import random
from rdkit import Chem
from rdkit.Chem import Draw
from pubchempy import get_compounds, PubChemPyError

import io
import base64
from PIL import Image

# for debug
from colorama import init, Fore, Style
init()
red = lambda txt: Fore.RED + txt + Style.RESET_ALL
blue = lambda txt: Fore.BLUE + txt + Style.RESET_ALL
green = lambda txt: Fore.GREEN + txt + Style.RESET_ALL

def general_append(smiles: list[str], longest_carbon_chain: int) -> None:
    # alkene (at least 2 C atoms) 1/3 probability to include
    if longest_carbon_chain >= 2 and random.randint(1, 3) == 1:
        no_of_alkenes = random.randint(1, longest_carbon_chain-1)
        for _ in range(no_of_alkenes):
            pos = random.randint(0, len(smiles)-2)
            # check invalid, also lower the probability of more alkenes
            if smiles[pos] == "C=C" or smiles[pos+1] == "C=C": # to avoid cumulene (consecutive double bonds)
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
        pos = random.randint(0 if 1 > longest_carbon_chain-1 else 1, longest_carbon_chain-1)
        smiles.insert(pos, f"({substituent})")

def random_smiles_easy() -> str:
    # meth ~ oct, normal distribution, max = 3
    longest_carbon_chain = random.choices([i for i in range(1, 9)], weights=[min([i, 9 - i, 3]) for i in range(1, 9)])[0]
    smiles = ["C"] * longest_carbon_chain

    cyclic = False
    # 1/4 probability cyclic
    if longest_carbon_chain >= 3 and random.randint(1, 4) == 1:
        cyclic = True
        smiles.insert(1, "1")
        smiles.append("1")

    general_append(smiles, longest_carbon_chain)    

    rnd = random.randint(1, 4)
    # either hydroxyl group or carboxyl group
    if rnd == 1:
        smiles.insert(random.randint(0, len(smiles)-1), "(O)")
    elif rnd == 2 and not cyclic:
        if random.randint(1, 2) == 1: # front
            smiles.insert(0, "OC(=O)")
        else: # back
            smiles.append("C(=O)O")

    smiles = "".join(smiles)

    # i have no idea how SMILES work
    if valid(smiles):
        return smiles
    else:
        # print("Invalid:", rnd, smiles)
        return random_smiles_easy()

def random_smiles_hard() -> str:
    # meth ~ oct, distribution prone to more carbons, no limit/max
    longest_carbon_chain = random.choices([i for i in range(1, 9)], weights=[min([i, 11 - i]) for i in range(1, 9)])[0]
    smiles = ["C"] * longest_carbon_chain

    # alkene 2/3 probability to include
    if longest_carbon_chain >= 2 and random.randint(1, 3) in (1, 2):
        # print(blue("Alkene:"))
        no_of_alkenes = random.randint(1, min(2, longest_carbon_chain-1)) # at most 2 alkenes
        for _ in range(no_of_alkenes):
            pos = random.randint(0, len(smiles)-2)
            if smiles[pos] == "C=C" or smiles[pos+1] == "C=C": # to avoid cumulene (consecutive double bonds)
                continue
            first = smiles.pop(pos)
            second = smiles.pop(pos)
            smiles.insert(pos, f"{first}={second}")

    if longest_carbon_chain > 2 and len(smiles) > 2:
        rnd = random.randint(1, 4)
    else: # no ketone
        rnd = random.randint(1, 3)
    # print(rnd)
    
    if rnd == 1: # aldehyde
        if smiles[0] == "C=C": # to avoid O=C=C, which is invalid
            # break the double bond
            smiles.pop(0)
            smiles = ["C", "C"] + smiles
        smiles.insert(0, "O=") # front
        if random.randint(1, 2) == 1: # 1/2 probability at the back
            if smiles[-1] == "C=C": # same logic
                smiles.pop()
                smiles = smiles + ["C", "C"]
            smiles.insert(len(smiles), "=O")
    
    elif rnd == 2: # amide
        smiles.pop()
        smiles.insert(0, "O=C(N)") # front
        if random.randint(1, 2) == 1: # 1/2 probability at the back
            smiles.insert(len(smiles), "(=O)N")
    
    else:
        # easier to draw less number of functional groups
        no_of_func_grps = random.choices([1, 2, 3], weights=[3, 2, 1])[0]

        for _ in range(no_of_func_grps):
            if rnd == 3: # amine
                smiles.insert(random.randint(1, len(smiles)), "(N)")
            elif rnd == 4: # ketone
                for _ in range(20): # run some time to avoid infinite loop, in case no valid position can be found
                    pos = random.randint(2, len(smiles)-1)
                    if smiles[pos-1] != "(=O)" and smiles[pos] != "(=O)": # avoid double ketone sticking to each other
                        break
                smiles.insert(pos, "(=O)")

    smiles = "".join(smiles)

    if valid(smiles):
        return smiles
    else:
        # print(red(f"Invalid: {rnd = }, {smiles = }"))
        return random_smiles_hard()

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
    # both modules are ready
    return iupac_name(smiles) is not None and Chem.MolFromSmiles(smiles) is not None

def generate_image(smiles: str) -> Image.Image:
    return Draw.MolToImage(Chem.MolFromSmiles(smiles))

def generate_base64_image(smiles: str) -> str:
    return image_to_base64(generate_image(smiles))

def image_to_base64(image: Image.Image) -> str:
    img_byte_arr = io.BytesIO()
    image.save(img_byte_arr, format="PNG")
    img_byte_arr = img_byte_arr.getvalue()
    img_base64 = base64.b64encode(img_byte_arr).decode("utf-8")
    return img_base64

# DEBUG
if __name__ == "__main__":
    ...

    for _ in range(10):
        smiles = random_smiles_hard()
        print(smiles, iupac_name(smiles))