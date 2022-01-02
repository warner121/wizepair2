from rdkit.Chem import MolFromSmiles, SaltRemover, MolToSmiles

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