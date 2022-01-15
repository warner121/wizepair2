from rdkit.Chem import MolFromSmiles, SaltRemover, RemoveStereochemistry, MolToSmiles

def strip_salts(smiles):
    
    # parse smiles as rdkit molecule
    try: mol = MolFromSmiles(smiles)
    except TypeError: return None
    
    # remove salts
    remover = SaltRemover.SaltRemover()
    mol, salts = remover.StripMolWithDeleted(mol)
    smiles = MolToSmiles(mol)
    
    # return
    return smiles

def strip_stereo(smiles):
    
    # parse smiles as rdkit molecule
    try: mol = MolFromSmiles(smiles)
    except TypeError: return None
    
    # remove stereochemistry
    RemoveStereochemistry(mol)
    
    #return
    return MolToSmiles(mol)