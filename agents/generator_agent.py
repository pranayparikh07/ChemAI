import numpy as np
import random
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit import DataStructs
import warnings

# Suppress RDKit and scikit-learn warnings
warnings.filterwarnings('ignore')
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# Suppress sklearn parallel warning
import sklearn.utils.parallel
original_warn = warnings.warn
def custom_warn(*args, **kwargs):
    if 'sklearn.utils.parallel.delayed' not in str(args):
        original_warn(*args, **kwargs)
warnings.warn = custom_warn


class GeneratorAgent:
    def __init__(self, seed_molecules=None):
        self.seed_molecules = seed_molecules or []
        self.mutation_operations = [
            'add_atom',
            'remove_atom',
            'change_atom',
            'add_bond',
            'remove_bond',
            'add_ring',
            'add_functional_group'
        ]
        
        self.functional_groups = {
            'hydroxyl': 'O',
            'amine': 'N',
            'carboxyl': 'C(=O)O',
            'methyl': 'C',
            'fluorine': 'F',
            'chlorine': 'Cl',
            'cyano': 'C#N',
            'methoxy': 'OC'
        }
        
        print(f"GeneratorAgent initialized with {len(self.seed_molecules)} seed molecules")
    
    def set_seed_molecules(self, smiles_list):
        valid_smiles = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                valid_smiles.append(smiles)
        self.seed_molecules = valid_smiles
        print(f"Set {len(valid_smiles)} valid seed molecules")
    
    def mutate_smiles(self, smiles, num_mutations=1):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        for _ in range(num_mutations):
            operation = random.choice(self.mutation_operations)
            mol = self._apply_mutation(mol, operation)
            if mol is None:
                return None
        
        try:
            Chem.SanitizeMol(mol)
            new_smiles = Chem.MolToSmiles(mol)
            # Validate the generated SMILES
            if Chem.MolFromSmiles(new_smiles) is not None:
                return new_smiles
        except:
            pass
        
        return None
    
    def _apply_mutation(self, mol, operation):
        try:
            rw_mol = Chem.RWMol(mol)
            
            if operation == 'add_atom':
                atom_types = [6, 7, 8, 9, 16, 17]
                new_atom_idx = rw_mol.AddAtom(Chem.Atom(random.choice(atom_types)))
                if rw_mol.GetNumAtoms() > 1:
                    existing_idx = random.randint(0, rw_mol.GetNumAtoms() - 2)
                    rw_mol.AddBond(existing_idx, new_atom_idx, Chem.BondType.SINGLE)
            
            elif operation == 'remove_atom' and rw_mol.GetNumAtoms() > 3:
                atom_idx = random.randint(0, rw_mol.GetNumAtoms() - 1)
                atom = rw_mol.GetAtomWithIdx(atom_idx)
                if atom.GetDegree() <= 1:
                    rw_mol.RemoveAtom(atom_idx)
            
            elif operation == 'change_atom':
                atom_types = [6, 7, 8, 16]
                atom_idx = random.randint(0, rw_mol.GetNumAtoms() - 1)
                rw_mol.GetAtomWithIdx(atom_idx).SetAtomicNum(random.choice(atom_types))
            
            elif operation == 'add_functional_group':
                fg_name = random.choice(list(self.functional_groups.keys()))
                fg_smiles = self.functional_groups[fg_name]
                fg_mol = Chem.MolFromSmiles(fg_smiles)
                if fg_mol:
                    combined = Chem.CombineMols(rw_mol, fg_mol)
                    rw_mol = Chem.RWMol(combined)
                    rw_mol.AddBond(0, mol.GetNumAtoms(), Chem.BondType.SINGLE)
            
            Chem.SanitizeMol(rw_mol)
            return rw_mol.GetMol()
        
        except:
            return mol
    
    def crossover(self, smiles1, smiles2):
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        if mol1 is None or mol2 is None:
            return None
        
        try:
            combined = Chem.CombineMols(mol1, mol2)
            rw_mol = Chem.RWMol(combined)
            
            n1 = mol1.GetNumAtoms()
            n2 = mol2.GetNumAtoms()
            
            if n1 > 0 and n2 > 0:
                idx1 = random.randint(0, n1 - 1)
                idx2 = random.randint(n1, n1 + n2 - 1)
                rw_mol.AddBond(idx1, idx2, Chem.BondType.SINGLE)
            
            Chem.SanitizeMol(rw_mol)
            return Chem.MolToSmiles(rw_mol.GetMol())
        except:
            return None
    
    def generate_analogs(self, smiles, num_analogs=10, similarity_threshold=0.5):
        analogs = []
        attempts = 0
        max_attempts = num_analogs * 10
        
        original_mol = Chem.MolFromSmiles(smiles)
        if original_mol is None:
            return []
        
        original_fp = AllChem.GetMorganFingerprintAsBitVect(original_mol, 2, nBits=1024)
        
        while len(analogs) < num_analogs and attempts < max_attempts:
            attempts += 1
            
            num_mutations = random.randint(1, 3)
            new_smiles = self.mutate_smiles(smiles, num_mutations)
            
            if new_smiles and new_smiles != smiles and new_smiles not in analogs:
                new_mol = Chem.MolFromSmiles(new_smiles)
                if new_mol:
                    new_fp = AllChem.GetMorganFingerprintAsBitVect(new_mol, 2, nBits=1024)
                    similarity = DataStructs.TanimotoSimilarity(original_fp, new_fp)
                    
                    if similarity >= similarity_threshold:
                        analogs.append({
                            'smiles': new_smiles,
                            'similarity': round(similarity, 3)
                        })
        
        return analogs
    
    def generate_from_seeds(self, num_per_seed=5):
        if not self.seed_molecules:
            print("No seed molecules available")
            return []
        
        all_generated = []
        
        for seed in self.seed_molecules:
            analogs = self.generate_analogs(seed, num_analogs=num_per_seed)
            for analog in analogs:
                analog['parent'] = seed
            all_generated.extend(analogs)
        
        print(f"Generated {len(all_generated)} molecules from {len(self.seed_molecules)} seeds")
        return all_generated
    
    def scaffold_hop(self, smiles, target_scaffold_smiles=None):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []
        
        try:
            from rdkit.Chem.Scaffolds import MurckoScaffold
            core = MurckoScaffold.GetScaffoldForMol(mol)
            core_smiles = Chem.MolToSmiles(core)
            
            hopped = []
            for _ in range(5):
                new_smiles = self.mutate_smiles(smiles, num_mutations=2)
                if new_smiles:
                    hopped.append(new_smiles)
            
            return hopped
        except:
            return []
    
    def generate_diverse_set(self, num_molecules=50, min_diversity=0.3):
        if not self.seed_molecules:
            print("No seed molecules for diverse generation")
            return []
        
        diverse_set = []
        all_fps = []
        
        while len(diverse_set) < num_molecules:
            seed = random.choice(self.seed_molecules)
            new_smiles = self.mutate_smiles(seed, num_mutations=random.randint(1, 4))
            
            if new_smiles:
                new_mol = Chem.MolFromSmiles(new_smiles)
                if new_mol:
                    new_fp = AllChem.GetMorganFingerprintAsBitVect(new_mol, 2, nBits=1024)
                    
                    is_diverse = True
                    for existing_fp in all_fps:
                        sim = DataStructs.TanimotoSimilarity(new_fp, existing_fp)
                        if sim > (1 - min_diversity):
                            is_diverse = False
                            break
                    
                    if is_diverse:
                        diverse_set.append(new_smiles)
                        all_fps.append(new_fp)
        
        print(f"Generated {len(diverse_set)} diverse molecules")
        return diverse_set


if __name__ == "__main__":
    agent = GeneratorAgent()
    
    seeds = [
        "CC(=O)Oc1ccccc1C(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    ]
    agent.set_seed_molecules(seeds)
    
    print("\n" + "=" * 60)
    print("Generator Agent Test")
    print("=" * 60)
    
    print("\n--- Generating analogs ---")
    analogs = agent.generate_analogs(seeds[0], num_analogs=5)
    for a in analogs:
        print(f"  {a['smiles']} (similarity: {a['similarity']})")
    
    print("\n--- Generating from all seeds ---")
    generated = agent.generate_from_seeds(num_per_seed=3)
    print(f"Total generated: {len(generated)}")
    
    print("\n--- Generating diverse set ---")
    diverse = agent.generate_diverse_set(num_molecules=10, min_diversity=0.3)
    print(f"Diverse molecules: {len(diverse)}")
