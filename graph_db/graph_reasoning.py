"""
Graph-based Reasoning Engine for Drug Discovery
Enables explanation and prediction using knowledge graph reasoning
"""

from typing import List, Dict, Tuple, Optional
from graph_algorithms import GraphAlgorithms


class ReasoningEngine:
    """Explains predictions and generates new hypotheses using graph reasoning"""
    
    def __init__(self, driver=None):
        """
        Initialize reasoning engine
        
        Args:
            driver: Neo4j driver instance
        """
        self.driver = driver
        self.algorithms = GraphAlgorithms(driver)
    
    # ============================================================
    # EXPLAINABILITY METHODS
    # ============================================================
    
    def explain_drug_target_prediction(self, molecule_id: str, target_id: str) -> Dict:
        """
        Explain why a molecule might target a protein
        
        Use case: Validate ML predictions with mechanistic reasoning
        
        Args:
            molecule_id: Molecule ChEMBL ID or SMILES
            target_id: Protein ChEMBL ID
            
        Returns:
            Dictionary with reasoning explanation
        """
        explanation = {
            'molecule_id': molecule_id,
            'target_id': target_id,
            'evidence': [],
            'confidence': 0.0,
            'mechanistic_reasoning': ''
        }
        
        if not self.driver:
            return explanation
        
        with self.driver.session() as session:
            # Check structural similarity to known binders
            similarity_result = session.run("""
                MATCH (m_query:Molecule)-[s:SIMILAR_TO]->(m_known:Molecule)
                      -[t:TARGETS]->(p:Protein)
                WHERE (m_query.smiles = $mol_id OR m_query.chembl_id = $mol_id)
                AND p.chembl_id = $target_id
                RETURN m_known.chembl_id as similar_mol, 
                       s.similarity_score as similarity,
                       t.pIC50 as known_activity
                ORDER BY similarity DESC
                LIMIT 5
            """, mol_id=molecule_id, target_id=target_id)
            
            for record in similarity_result:
                explanation['evidence'].append({
                    'type': 'structural_analogy',
                    'description': f"Similar to known binder {record['similar_mol']} " +
                                  f"(similarity: {record['similarity']:.2f}, " +
                                  f"pIC50: {record['known_activity']:.1f})"
                })
        
            # Check pathway participation
            pathway_result = session.run("""
                MATCH (m:Molecule)-[t:TARGETS]->(p1:Protein)
                      -[pp:PARTICIPATES_IN]->(pw:Pathway)
                      <-[ppp:PARTICIPATES_IN]-(p2:Protein)
                WHERE (m.smiles = $mol_id OR m.chembl_id = $mol_id)
                AND p2.chembl_id = $target_id
                RETURN pw.name as pathway, count(DISTINCT p1) as interconnection
                ORDER BY interconnection DESC
            """, mol_id=molecule_id, target_id=target_id)
            
            for record in pathway_result:
                explanation['evidence'].append({
                    'type': 'pathway_connectivity',
                    'description': f"Participates in {record['pathway']} pathway " +
                                  f"(interconnection degree: {record['interconnection']})"
                })
            
            # Calculate confidence
            explanation['confidence'] = min(0.95, len(explanation['evidence']) * 0.25)
            
            # Generate reasoning
            if explanation['evidence']:
                explanation['mechanistic_reasoning'] = (
                    f"Based on {len(explanation['evidence'])} pieces of evidence: "
                    f"(1) Structural analogy to known binders, "
                    f"(2) Pathway-level connectivity through related proteins"
                )
            else:
                explanation['mechanistic_reasoning'] = (
                    "No direct mechanistic evidence found in knowledge graph"
                )
        
        return explanation
    
    def explain_bioactivity_prediction(self, molecule_id: str) -> Dict:
        """
        Explain predicted bioactivity/potency based on graph structure
        
        Args:
            molecule_id: Molecule ChEMBL ID or SMILES
            
        Returns:
            Explanation with contributing factors
        """
        explanation = {
            'molecule_id': molecule_id,
            'factors': [],
            'predicted_mechanism': ''
        }
        
        if not self.driver:
            return explanation
        
        with self.driver.session() as session:
            # Find similar active compounds
            similar = session.run("""
                MATCH (m_query:Molecule)-[s:SIMILAR_TO]->(m_known:Molecule)
                      -[t:TARGETS]->(p:Protein)
                WHERE m_query.smiles = $mol_id OR m_query.chembl_id = $mol_id
                RETURN m_known.chembl_id as similar_mol,
                       avg(t.pIC50) as avg_activity,
                       s.similarity_score as similarity,
                       count(DISTINCT p) as target_count
                ORDER BY similarity DESC
                LIMIT 10
            """, mol_id=molecule_id)
            
            activities = []
            for record in similar:
                activities.append(record['avg_activity'])
                explanation['factors'].append({
                    'factor': f"Similar to active compound {record['similar_mol']}",
                    'similarity': record['similarity'],
                    'reference_activity': record['avg_activity'],
                    'target_count': record['target_count']
                })
            
            if activities:
                avg_activity = sum(activities) / len(activities)
                explanation['predicted_mechanism'] = (
                    f"High bioactivity expected (estimated pIC50: {avg_activity:.1f}) "
                    f"based on strong structural analogs in the database"
                )
        
        return explanation
    
    # ============================================================
    # HYPOTHESIS GENERATION METHODS
    # ============================================================
    
    def generate_repurposing_hypothesis(self, molecule_id: str) -> Dict:
        """
        Generate drug repurposing hypothesis for a molecule
        
        Use case: Identify new indications for existing drugs
        
        Args:
            molecule_id: Molecule (drug) ChEMBL ID
            
        Returns:
            Repurposing hypothesis with paths
        """
        hypothesis = {
            'molecule_id': molecule_id,
            'known_indications': [],
            'repurposing_candidates': [],
            'mechanistic_support': []
        }
        
        if not self.driver:
            return hypothesis
        
        with self.driver.session() as session:
            # Get known targets
            targets = session.run("""
                MATCH (m:Molecule)-[t:TARGETS]->(p:Protein)
                WHERE m.smiles = $mol_id OR m.chembl_id = $mol_id
                RETURN p.chembl_id as target_id, p.name, t.pIC50
                ORDER BY t.pIC50 DESC
                LIMIT 5
            """, mol_id=molecule_id)
            
            target_list = []
            for record in targets:
                target_list.append(record['target_id'])
                
                # Find diseases for this target
                diseases = session.run("""
                    MATCH (p:Protein)-[td:TARGETS_DISEASE]->(d:Disease)
                    WHERE p.chembl_id = $target_id
                    RETURN d.name as disease, td.association_strength as strength
                    ORDER BY strength DESC
                """, target_id=record['target_id'])
                
                for disease in diseases:
                    hypothesis['known_indications'].append({
                        'disease': disease['disease'],
                        'target': record['name'],
                        'target_strength': record['pIC50'],
                        'association_strength': disease['strength']
                    })
            
            # Find multi-hop disease opportunities
            if target_list:
                repurpose = session.run("""
                    MATCH (m:Molecule)-[t:TARGETS]->(p1:Protein)
                          -[pp:PARTICIPATES_IN]->(pw:Pathway)
                          <-[ppp:PARTICIPATES_IN]-(p2:Protein)
                          -[td:TARGETS_DISEASE]->(d:Disease)
                    WHERE (m.smiles = $mol_id OR m.chembl_id = $mol_id)
                    AND NOT (m)-[:TARGETS_DISEASE]->(d)
                    WITH d.name as disease, p2.name as new_target, 
                         pw.name as pathway, count(DISTINCT p1) as path_count
                    RETURN disease, collect(new_target) as targets, 
                           pathway, path_count
                    ORDER BY path_count DESC
                    LIMIT 5
                """, mol_id=molecule_id)
                
                for record in repurpose:
                    hypothesis['repurposing_candidates'].append({
                        'disease': record['disease'],
                        'targets': record['targets'],
                        'pathway': record['pathway'],
                        'mechanistic_support': f"Multi-target coverage of {record['pathway']} pathway"
                    })
        
        return hypothesis
    
    def generate_combination_therapy_hypothesis(self, molecule_ids: List[str]) -> Dict:
        """
        Generate combination therapy hypothesis for multiple molecules
        
        Use case: Synergistic drug combinations
        
        Args:
            molecule_ids: List of molecule IDs or SMILES
            
        Returns:
            Combination therapy hypothesis
        """
        hypothesis = {
            'molecules': molecule_ids,
            'target_coverage': [],
            'pathway_coverage': [],
            'synergy_potential': 0.0,
            'mechanism': ''
        }
        
        if not self.driver or len(molecule_ids) < 2:
            return hypothesis
        
        with self.driver.session() as session:
            # Compute coverage
            for mol_id in molecule_ids:
                targets = session.run("""
                    MATCH (m:Molecule)-[t:TARGETS]->(p:Protein)
                    WHERE m.smiles = $mol_id OR m.chembl_id = $mol_id
                    RETURN p.chembl_id as target_id, t.pIC50
                    ORDER BY t.pIC50 DESC
                """, mol_id=mol_id)
                
                mol_targets = [record['target_id'] for record in targets]
                hypothesis['target_coverage'].append(len(mol_targets))
            
            # Find synergistic pathways
            if len(molecule_ids) >= 2:
                synergy = session.run("""
                    MATCH (m1:Molecule)-[t1:TARGETS]->(p1:Protein)
                          -[pp:PARTICIPATES_IN]->(pw:Pathway)
                          <-[ppp:PARTICIPATES_IN]-(p2:Protein)
                          <-[t2:TARGETS]-(m2:Molecule)
                    WHERE (m1.smiles = $mol1 OR m1.chembl_id = $mol1)
                    AND (m2.smiles = $mol2 OR m2.chembl_id = $mol2)
                    WITH pw.name as pathway, 
                         collect(DISTINCT p1.name) as targets1,
                         collect(DISTINCT p2.name) as targets2,
                         count(*) as connections
                    RETURN pathway, targets1, targets2, connections
                    ORDER BY connections DESC
                """, mol1=molecule_ids[0], mol2=molecule_ids[1])
                
                for record in synergy:
                    hypothesis['pathway_coverage'].append({
                        'pathway': record['pathway'],
                        'molecule1_targets': record['targets1'],
                        'molecule2_targets': record['targets2']
                    })
            
            # Calculate synergy potential
            hypothesis['synergy_potential'] = min(
                0.95, 
                (len(hypothesis['target_coverage']) * 0.3 + 
                 len(hypothesis['pathway_coverage']) * 0.4)
            )
            
            if hypothesis['pathway_coverage']:
                pathway_list = [p['pathway'] for p in hypothesis['pathway_coverage']]
                hypothesis['mechanism'] = (
                    f"Multi-target synergy via {', '.join(pathway_list[:3])} pathways"
                )
        
        return hypothesis
    
    # ============================================================
    # VALIDATION METHODS
    # ============================================================
    
    def validate_prediction_with_graph(self, prediction: Dict) -> Dict:
        """
        Validate ML prediction using graph structure
        
        Args:
            prediction: ML model prediction (molecule, target, activity)
            
        Returns:
            Validation result with confidence adjustment
        """
        validation = {
            'original_prediction': prediction,
            'graph_evidence': [],
            'confidence_adjustment': 0.0,
            'final_confidence': prediction.get('confidence', 0.5)
        }
        
        if not self.driver:
            return validation
        
        molecule_id = prediction.get('molecule_id')
        target_id = prediction.get('target_id')
        predicted_activity = prediction.get('predicted_pIC50', 0.0)
        
        # Get explanations
        explanation = self.explain_drug_target_prediction(molecule_id, target_id)
        
        validation['graph_evidence'] = explanation['evidence']
        validation['confidence_adjustment'] = explanation['confidence'] - prediction.get('confidence', 0.5)
        validation['final_confidence'] = min(
            0.99,
            prediction.get('confidence', 0.5) + validation['confidence_adjustment']
        )
        
        return validation
    
    # ============================================================
    # UTILITY METHODS
    # ============================================================
    
    def get_graph_statistics(self) -> Dict:
        """Get knowledge graph statistics"""
        stats = {
            'molecules': 0,
            'proteins': 0,
            'diseases': 0,
            'pathways': 0,
            'interactions': 0
        }
        
        if not self.driver:
            return stats
        
        with self.driver.session() as session:
            for node_type in ['Molecule', 'Protein', 'Disease', 'Pathway']:
                result = session.run(f"MATCH (n:{node_type}) RETURN count(n) as count")
                count = result.single()['count']
                stats[node_type.lower() + 's'] = count
            
            # Count interactions
            result = session.run("""
                MATCH ()-[r:TARGETS]->()
                RETURN count(r) as count
            """)
            stats['interactions'] = result.single()['count']
        
        return stats
    
    def close(self):
        """Clean up resources"""
        self.algorithms.close()
        if self.driver:
            self.driver.close()


if __name__ == "__main__":
    print("Reasoning Engine for Drug Discovery")
    print("Provides explainability and hypothesis generation")
