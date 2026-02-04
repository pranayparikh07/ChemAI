#!/usr/bin/env python3
"""
Graph-based Reasoning Engine for ChemAI
Provides path-based explanations and mechanistic inference
"""

import json
from typing import Dict, List, Tuple
from collections import defaultdict, deque


class KnowledgeGraphReasoning:
    """Reasoning engine for knowledge graph queries and explanations"""
    
    def __init__(self):
        self.graph = defaultdict(list)  # adjacency list
        self.nodes = {}
        self.edges = {}
    
    def add_node(self, node_id, node_data):
        """Add node to graph"""
        self.nodes[node_id] = node_data
    
    def add_edge(self, source, target, edge_data):
        """Add edge to graph"""
        self.graph[source].append((target, edge_data))
        self.edges[f"{source}->{target}"] = edge_data
    
    def find_paths(self, start, end, max_length=3):
        """Find all paths from start to end node"""
        paths = []
        visited = set()
        
        def dfs(current, target, path, depth):
            if depth > max_length:
                return
            
            if current == target:
                paths.append(path)
                return
            
            for neighbor, edge_data in self.graph[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    new_path = path + [(current, neighbor, edge_data)]
                    dfs(neighbor, target, new_path, depth + 1)
                    visited.remove(neighbor)
        
        visited.add(start)
        dfs(start, end, [], 0)
        return paths
    
    def explain_bioactivity(self, molecule_id, protein_id):
        """Generate mechanistic explanation for bioactivity"""
        explanation = {
            'molecule': molecule_id,
            'protein': protein_id,
            'reasoning_paths': [],
            'similar_compounds': [],
            'mechanisms': []
        }
        
        # Find direct path
        paths = self.find_paths(molecule_id, protein_id, max_length=2)
        for path in paths:
            path_description = self._describe_path(path)
            explanation['reasoning_paths'].append(path_description)
        
        # Find similar molecules with known activity
        similar = self._find_similar_molecules(molecule_id)
        explanation['similar_compounds'] = similar
        
        # Identify mechanisms
        mechanisms = self._identify_mechanisms(molecule_id, protein_id)
        explanation['mechanisms'] = mechanisms
        
        return explanation
    
    def _describe_path(self, path):
        """Convert path to human-readable description"""
        description = []
        for source, target, edge_data in path:
            edge_type = edge_data.get('type', 'UNKNOWN')
            confidence = edge_data.get('confidence', 0.5)
            description.append(
                f"{source} --[{edge_type}]--> {target} (confidence: {confidence:.2f})"
            )
        return " ".join(description)
    
    def _find_similar_molecules(self, molecule_id, top_k=5):
        """Find similar molecules with known bioactivity"""
        similar = []
        for neighbor, edge_data in self.graph[molecule_id]:
            if edge_data.get('type') == 'SIMILAR_TO':
                similarity = edge_data.get('similarity', 0)
                similar.append({
                    'molecule': neighbor,
                    'similarity': similarity,
                    'bioactivity': 'N/A (would query targets)'
                })
        return sorted(similar, key=lambda x: x['similarity'], reverse=True)[:top_k]
    
    def _identify_mechanisms(self, molecule_id, protein_id):
        """Identify bioactivity mechanisms"""
        mechanisms = []
        
        # Based on molecular properties
        mol_data = self.nodes.get(molecule_id, {})
        if 'QED' in mol_data.get('properties', {}):
            qed = mol_data['properties']['QED']
            if qed > 0.7:
                mechanisms.append("High drug-likeness favors favorable pharmacokinetics")
        
        # Based on protein properties
        prot_data = self.nodes.get(protein_id, {})
        if 'target_type' in prot_data.get('properties', {}):
            mechanisms.append(f"Target type: {prot_data['properties']['target_type']}")
        
        return mechanisms
    
    def predict_with_reasoning(self, molecule_id, objective='bioactivity'):
        """Make prediction with reasoning chain"""
        result = {
            'molecule': molecule_id,
            'objective': objective,
            'prediction': 'N/A',
            'reasoning': [],
            'confidence': 0.0
        }
        
        mol_data = self.nodes.get(molecule_id, {})
        if not mol_data:
            return result
        
        # Rule-based prediction based on properties
        props = mol_data.get('properties', {})
        
        if objective == 'bioactivity':
            # High QED → better bioactivity
            qed = props.get('QED', 0)
            mw = props.get('MW', 0)
            logp = props.get('LogP', 0)
            
            result['reasoning'].append(f"QED score: {qed:.3f}")
            result['reasoning'].append(f"Molecular weight: {mw:.1f} Da")
            result['reasoning'].append(f"LogP: {logp:.2f}")
            
            # Simple heuristic
            if qed > 0.6 and 200 < mw < 600 and -1 < logp < 5:
                result['prediction'] = "HIGH bioactivity expected"
                result['confidence'] = 0.75
            else:
                result['prediction'] = "MODERATE bioactivity expected"
                result['confidence'] = 0.55
        
        return result


class PathBasedExplainer:
    """Generate textual explanations for predictions"""
    
    @staticmethod
    def explain_prediction(prediction_data, kg_reasoning):
        """Generate human-readable explanation"""
        explanation = {
            'summary': '',
            'details': [],
            'confidence': prediction_data.get('confidence', 0),
            'factors': []
        }
        
        # Generate summary
        molecule = prediction_data.get('molecule', 'Unknown')
        objective = prediction_data.get('objective', 'bioactivity')
        prediction = prediction_data.get('prediction', 'N/A')
        
        explanation['summary'] = f"{molecule}: {prediction}"
        
        # Add reasoning factors
        for reason in prediction_data.get('reasoning', []):
            explanation['factors'].append(reason)
        
        # Add details
        explanation['details'] = [
            f"Prediction confidence: {100*explanation['confidence']:.1f}%",
            f"Based on: {objective}",
        ]
        
        return explanation


def load_graph_from_csv(molecules_csv, proteins_csv, bioactivity_csv, similarity_csv):
    """Load KG from CSV files"""
    kg = KnowledgeGraphReasoning()
    
    # Load molecules
    import pandas as pd
    mol_df = pd.read_csv(molecules_csv)
    for idx, row in mol_df.iterrows():
        node_data = {
            'type': 'Molecule',
            'properties': row.to_dict()
        }
        kg.add_node(row['id'], node_data)
    
    # Load proteins
    prot_df = pd.read_csv(proteins_csv)
    for idx, row in prot_df.iterrows():
        node_data = {
            'type': 'Protein',
            'properties': row.to_dict()
        }
        kg.add_node(row['id'], node_data)
    
    # Load bioactivity edges
    bio_df = pd.read_csv(bioactivity_csv)
    for idx, row in bio_df.iterrows():
        edge_data = {
            'type': 'TARGETS',
            'pic50': row.get('pic50'),
            'confidence': row.get('confidence', 0.85),
        }
        kg.add_edge(row['source'], row['target'], edge_data)
    
    # Load similarity edges
    sim_df = pd.read_csv(similarity_csv)
    for idx, row in sim_df.iterrows():
        edge_data = {
            'type': 'SIMILAR_TO',
            'similarity': row.get('similarity'),
        }
        kg.add_edge(row['source'], row['target'], edge_data)
    
    return kg


if __name__ == "__main__":
    # Example usage
    kg = load_graph_from_csv(
        'd:\\ChemAI\\TEAM_SHREYA\\data\\molecules_kg.csv',
        'd:\\ChemAI\\TEAM_SHREYA\\data\\proteins_kg.csv',
        'd:\\ChemAI\\TEAM_SHREYA\\data\\bioactivity_edges.csv',
        'd:\\ChemAI\\TEAM_SHREYA\\data\\similarity_edges.csv'
    )
    
    print("✓ Knowledge graph loaded and ready for reasoning!")
    print(f"  Nodes: {len(kg.nodes)}")
    print(f"  Edges: {len(kg.edges)}")
