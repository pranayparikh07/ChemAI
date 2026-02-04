"""
Graph Algorithms for Drug Discovery Reasoning
Implements graph-based algorithms for molecule discovery and target ranking
"""

from typing import List, Dict, Tuple, Set
import numpy as np
from collections import defaultdict, deque
import heapq


class GraphAlgorithms:
    """Core graph algorithms for knowledge graph reasoning"""
    
    def __init__(self, driver=None):
        """
        Initialize algorithms
        
        Args:
            driver: Neo4j driver instance
        """
        self.driver = driver
    
    # ============================================================
    # SHORTEST PATH ALGORITHMS (Drug Repurposing)
    # ============================================================
    
    def find_shortest_path(self, start_node: str, end_node: str, 
                          max_depth: int = 5) -> List[Dict]:
        """
        Find shortest path between two nodes
        
        Use case: Drug repurposing - find path from known drug to new disease
        
        Args:
            start_node: Starting node ID (e.g., molecule ChEMBL ID)
            end_node: Target node ID (e.g., disease ID)
            max_depth: Maximum path length
            
        Returns:
            List of shortest paths
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run("""
                MATCH p = shortestPath((n1)-[*..{max_depth}]-(n2))
                WHERE n1.chembl_id = $start OR n1.smiles = $start
                AND (n2.disease_id = $end OR n2.mesh_id = $end)
                RETURN p
                LIMIT 10
            """, start=start_node, end=end_node, max_depth=max_depth)
            
            return [record['p'] for record in result]
    
    def find_multi_hop_targets(self, molecule_id: str) -> Dict:
        """
        Find targets via multiple paths (protein→disease, protein→pathway, etc.)
        
        Use case: Understanding drug mechanism of action
        
        Args:
            molecule_id: Molecule ChEMBL ID or SMILES
            
        Returns:
            Dictionary of paths and targets
        """
        if not self.driver:
            return {}
        
        results = {
            'direct_targets': [],
            'targets_via_pathways': [],
            'targets_via_diseases': [],
            'potential_side_effects': []
        }
        
        with self.driver.session() as session:
            # Direct targets
            direct = session.run("""
                MATCH (m:Molecule)-[r:TARGETS]->(p:Protein)
                WHERE m.smiles = $mol_id OR m.chembl_id = $mol_id
                RETURN p.chembl_id as target_id, p.name, r.pIC50, r.ic50_nM
                ORDER BY r.pIC50 DESC
            """, mol_id=molecule_id)
            
            results['direct_targets'] = [dict(record) for record in direct]
            
            # Targets via pathways
            pathways = session.run("""
                MATCH (m:Molecule)-[t:TARGETS]->(p:Protein)
                      -[pp:PARTICIPATES_IN]->(pw:Pathway)
                      <-[tp:PARTICIPATES_IN]-(p2:Protein)
                WHERE m.smiles = $mol_id OR m.chembl_id = $mol_id
                AND p <> p2
                RETURN DISTINCT p2.chembl_id as target_id, p2.name, 
                       pw.name as pathway, t.pIC50, count(*) as strength
                ORDER BY strength DESC
            """, mol_id=molecule_id)
            
            results['targets_via_pathways'] = [dict(record) for record in pathways]
            
            # Potential side effects via protein similarity
            side_effects = session.run("""
                MATCH (m:Molecule)-[t:TARGETS]->(p:Protein)
                      <-[se:CAUSES]-(s:SideEffect)
                WHERE m.smiles = $mol_id OR m.chembl_id = $mol_id
                RETURN se.name as side_effect, s.severity, s.frequency, 
                       p.name as protein
                ORDER BY s.severity DESC
            """, mol_id=molecule_id)
            
            results['potential_side_effects'] = [dict(record) for record in side_effects]
        
        return results
    
    # ============================================================
    # CENTRALITY ALGORITHMS (Target Importance Ranking)
    # ============================================================
    
    def compute_degree_centrality(self, node_type: str = "Protein") -> List[Tuple]:
        """
        Compute degree centrality for nodes
        
        Use case: Identify hub proteins (druggable, promiscuous targets)
        
        Args:
            node_type: Type of node (Protein, Gene, Disease, Pathway)
            
        Returns:
            List of (node_id, node_name, degree) sorted by degree
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run(f"""
                MATCH (n:{node_type})
                WITH n, size(()-[]-(n)) as degree
                RETURN n.chembl_id as node_id, n.name as node_name, degree
                ORDER BY degree DESC
                LIMIT 100
            """)
            
            return [(record['node_id'], record['node_name'], record['degree']) 
                    for record in result]
    
    def compute_betweenness_centrality(self, node_type: str = "Protein") -> List[Tuple]:
        """
        Compute betweenness centrality for nodes
        
        Use case: Identify bottleneck proteins connecting pathways
        
        Args:
            node_type: Type of node (Protein, Gene, etc.)
            
        Returns:
            List of (node_id, node_name, betweenness_score) sorted by score
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run(f"""
                MATCH (n:{node_type})
                WITH n, 
                     size((m:Molecule)-[]-(n)) * size((d:Disease)-[]-(n)) 
                     as betweenness_score
                RETURN n.chembl_id as node_id, n.name as node_name, 
                       betweenness_score
                ORDER BY betweenness_score DESC
                LIMIT 100
            """)
            
            return [(record['node_id'], record['node_name'], 
                    record['betweenness_score']) for record in result]
    
    def rank_targets_by_importance(self, disease_id: str) -> List[Dict]:
        """
        Rank protein targets by importance for a disease
        
        Use case: Identify best targets for drug development
        
        Args:
            disease_id: Disease ID or name
            
        Returns:
            Ranked list of targets with scores
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run("""
                MATCH (d:Disease)<-[r:TARGETS_DISEASE]-(p:Protein)
                WHERE d.disease_id = $disease_id OR d.name = $disease_id
                WITH p, r.association_strength as strength,
                     size((m:Molecule)-[t:TARGETS]->(p)) as num_ligands,
                     size((p)-[]-(:Pathway)) as pathway_count
                WITH p, strength, num_ligands, pathway_count,
                     strength * (1 + log10(num_ligands + 1)) * 
                     log10(pathway_count + 1) as importance_score
                RETURN p.chembl_id as target_id, p.name as target_name,
                       strength, num_ligands, pathway_count, importance_score
                ORDER BY importance_score DESC
                LIMIT 20
            """, disease_id=disease_id)
            
            return [dict(record) for record in result]
    
    # ============================================================
    # SIMILARITY ALGORITHMS (Compound/Target Similarity)
    # ============================================================
    
    def find_similar_compounds(self, molecule_id: str, 
                              min_similarity: float = 0.7,
                              limit: int = 20) -> List[Dict]:
        """
        Find similar compounds based on structural similarity
        
        Use case: SAR (structure-activity relationship) analysis
        
        Args:
            molecule_id: Query molecule ID or SMILES
            min_similarity: Minimum Tanimoto similarity threshold
            limit: Maximum results to return
            
        Returns:
            List of similar molecules with similarity scores
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run("""
                MATCH (m1:Molecule {smiles: $mol_id})-[s:SIMILAR_TO]->(m2:Molecule)
                WHERE s.similarity_score >= $min_sim
                RETURN m2.chembl_id as chembl_id, m2.smiles, m2.name,
                       s.similarity_score as similarity,
                       size((m2)-[t:TARGETS]->(:Protein)) as num_targets
                ORDER BY similarity DESC
                LIMIT $limit
            """, mol_id=molecule_id, min_sim=min_similarity, limit=limit)
            
            return [dict(record) for record in result]
    
    def find_target_analogs(self, target_id: str, 
                           similarity_measure: str = "sequence") -> List[Dict]:
        """
        Find similar targets (orthologs, paralogs, or functionally similar)
        
        Use case: Cross-species drug efficacy prediction
        
        Args:
            target_id: Protein ChEMBL ID
            similarity_measure: 'sequence' or 'functional'
            
        Returns:
            List of similar targets
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run("""
                MATCH (p1:Protein {chembl_id: $target_id})-[s:SIMILAR_TO]->(p2:Protein)
                WITH p1, p2, s.similarity_score as similarity
                RETURN p2.chembl_id as target_id, p2.name, p2.species,
                       similarity,
                       size((m:Molecule)-[t:TARGETS]->(p2)) as num_ligands
                ORDER BY similarity DESC, num_ligands DESC
                LIMIT 20
            """, target_id=target_id)
            
            return [dict(record) for record in result]
    
    # ============================================================
    # CLUSTERING ALGORITHMS (Disease Modules, Drug Classes)
    # ============================================================
    
    def find_disease_modules(self, min_module_size: int = 3) -> List[Dict]:
        """
        Find disease modules (protein clusters associated with diseases)
        
        Use case: Systems pharmacology, disease mechanism understanding
        
        Args:
            min_module_size: Minimum proteins in a module
            
        Returns:
            List of disease modules with proteins and connections
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run("""
                MATCH (d:Disease)<-[r:TARGETS_DISEASE]-(p:Protein)
                MATCH (p)-[]->(pw:Pathway)<-[]-(p2:Protein)
                WHERE p <> p2
                WITH d.name as disease, collect(DISTINCT p.chembl_id) as proteins,
                     count(DISTINCT pw) as pathway_count, count(DISTINCT p2) as degree
                WHERE size(proteins) >= $min_size
                RETURN disease, proteins, pathway_count, degree
                ORDER BY degree DESC
            """, min_size=min_module_size)
            
            return [dict(record) for record in result]
    
    # ============================================================
    # POLYPHARMACOLOGY ALGORITHMS (Multi-target Effects)
    # ============================================================
    
    def find_drug_synergies(self, molecule_id: str) -> List[Dict]:
        """
        Find synergistic drug combinations based on pathway coverage
        
        Use case: Combination therapy discovery
        
        Args:
            molecule_id: Query molecule ID or SMILES
            
        Returns:
            List of synergistic molecule combinations
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run("""
                MATCH (m1:Molecule {smiles: $mol_id})-[t1:TARGETS]->(p1:Protein)
                      -[pp1:PARTICIPATES_IN]->(pw:Pathway)
                      <-[pp2:PARTICIPATES_IN]-(p2:Protein)
                      <-[t2:TARGETS]-(m2:Molecule)
                WHERE m1 <> m2
                MATCH (p1)-[pp:PARTICIPATES_IN]->(pw2:Pathway)
                      <-[ppp:PARTICIPATES_IN]-(p2)
                WITH m1, m2, pw, pw2, 
                     count(DISTINCT p1) + count(DISTINCT p2) as coverage,
                     count(DISTINCT pw) + count(DISTINCT pw2) as pathway_count
                RETURN m2.chembl_id as combo_molecule, m2.name, m2.smiles,
                       coverage, pathway_count,
                       coverage * pathway_count as synergy_score
                ORDER BY synergy_score DESC
                LIMIT 20
            """, mol_id=molecule_id)
            
            return [dict(record) for record in result]
    
    # ============================================================
    # COMMUNITY DETECTION (Drug Classes)
    # ============================================================
    
    def find_compound_clusters(self, min_cluster_size: int = 5) -> List[Dict]:
        """
        Find clusters of similar compounds based on structure and targets
        
        Use case: Drug class identification, SAR analysis
        
        Args:
            min_cluster_size: Minimum compounds in cluster
            
        Returns:
            List of compound clusters
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run("""
                MATCH (m1:Molecule)-[s1:SIMILAR_TO]->(m2:Molecule)
                      -[s2:SIMILAR_TO]->(m3:Molecule)
                WHERE m1 <> m3
                WITH m1, collect(DISTINCT m2.chembl_id) as neighbors
                WHERE size(neighbors) >= $min_size
                RETURN m1.chembl_id as cluster_center, m1.name,
                       neighbors as cluster_members,
                       size(neighbors) as cluster_size
                ORDER BY cluster_size DESC
            """, min_size=min_cluster_size)
            
            return [dict(record) for record in result]
    
    # ============================================================
    # RECOMMENDATION ALGORITHMS
    # ============================================================
    
    def recommend_new_targets(self, molecule_id: str, 
                             disease_id: str = None) -> List[Dict]:
        """
        Recommend unexplored targets for a molecule based on graph structure
        
        Use case: Virtual target screening
        
        Args:
            molecule_id: Molecule ID or SMILES
            disease_id: Optional disease context
            
        Returns:
            Ranked list of recommended targets
        """
        if not self.driver:
            return []
        
        with self.driver.session() as session:
            result = session.run("""
                MATCH (m:Molecule)-[t:TARGETS]->(p1:Protein)
                WHERE m.smiles = $mol_id OR m.chembl_id = $mol_id
                MATCH (p1)-[pp:PARTICIPATES_IN]->(pw:Pathway)
                      <-[ppp:PARTICIPATES_IN]-(p2:Protein)
                WHERE NOT (m)-[:TARGETS]->(p2)
                WITH p2, count(DISTINCT pw) as pathway_overlap,
                     size((m2:Molecule)-[t2:TARGETS]->(p2)) as existing_ligands
                RETURN p2.chembl_id as target_id, p2.name, p2.gene_name,
                       pathway_overlap, existing_ligands,
                       pathway_overlap * (1 + log10(existing_ligands + 1)) as score
                ORDER BY score DESC
                LIMIT 20
            """, mol_id=molecule_id)
            
            return [dict(record) for record in result]
    
    def close(self):
        """Clean up resources"""
        if self.driver:
            self.driver.close()


if __name__ == "__main__":
    print("Graph Algorithms Module")
    print("Provides algorithms for drug discovery reasoning")
