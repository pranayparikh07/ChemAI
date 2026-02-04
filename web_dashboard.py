    """
Web Dashboard for ChemAI Drug Discovery Results
Displays molecule discovery results in a beautiful HTML interface with real-time tracking
"""

import json
from datetime import datetime
from pathlib import Path


class WebDashboard:
    def __init__(self, output_file='drug_discovery_dashboard.html'):
        self.output_file = output_file
        self.molecules = []
        self.start_time = datetime.now()
        self.generation_count = 0
        
    def add_molecule(self, smiles, name=None, properties=None, predictions=None, generation=1):
        """Add a discovered molecule to the dashboard"""
        # Flatten nested predictions structure
        flat_predictions = {}
        if predictions:
            # Extract nested values from bioactivity and druglikeness
            if 'bioactivity' in predictions and isinstance(predictions['bioactivity'], dict):
                flat_predictions['pIC50'] = predictions['bioactivity'].get('pIC50', None)
            if 'druglikeness' in predictions and isinstance(predictions['druglikeness'], dict):
                flat_predictions['qed_score'] = predictions['druglikeness'].get('qed_score', None)
                flat_predictions['is_druglike'] = predictions['druglikeness'].get('is_druglike', None)
            # Preserve any top-level keys
            for key, value in predictions.items():
                if key not in ['bioactivity', 'druglikeness', 'toxicity', 'properties']:
                    flat_predictions[key] = value
        
        mol_data = {
            'id': len(self.molecules) + 1,
            'smiles': smiles,
            'name': name or f'Compound_{len(self.molecules) + 1}',
            'generation': generation,
            'timestamp': datetime.now().isoformat(),
            'properties': properties or {},
            'predictions': flat_predictions,
            'full_predictions': predictions or {},  # Keep full structure for reference
        }
        self.molecules.append(mol_data)
        
        # Print to console
        print(f"✓ Molecule Created: {mol_data['name']}")
        print(f"  SMILES: {smiles}")
        if flat_predictions.get('pIC50') is not None:
            pic50_val = flat_predictions['pIC50']
            print(f"  pIC50: {pic50_val if isinstance(pic50_val, str) else f'{pic50_val:.3f}'}")
        if flat_predictions.get('qed_score') is not None:
            qed_val = flat_predictions['qed_score']
            print(f"  QED: {qed_val if isinstance(qed_val, str) else f'{qed_val:.3f}'}")
        print()
        
        return mol_data
    
    def set_generation(self, generation):
        """Set current generation number"""
        self.generation_count = generation
    
    def generate_html(self):
        """Generate HTML dashboard"""
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChemAI Drug Discovery Dashboard</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        
        .header {{
            background: white;
            border-radius: 10px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
        }}
        
        .header h1 {{
            color: #667eea;
            margin-bottom: 10px;
            font-size: 2.5em;
        }}
        
        .header p {{
            color: #666;
            font-size: 1.1em;
        }}
        
        .stats {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        
        .stat-card h3 {{
            font-size: 2em;
            margin-bottom: 5px;
        }}
        
        .stat-card p {{
            opacity: 0.9;
        }}
        
        .molecules-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .molecule-card {{
            background: white;
            border-radius: 10px;
            padding: 20px;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
        }}
        
        .molecule-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 15px 30px rgba(0,0,0,0.2);
        }}
        
        .molecule-id {{
            background: #667eea;
            color: white;
            padding: 8px 12px;
            border-radius: 5px;
            display: inline-block;
            margin-bottom: 10px;
            font-weight: bold;
        }}
        
        .molecule-name {{
            font-size: 1.5em;
            color: #333;
            margin-bottom: 10px;
            font-weight: bold;
        }}
        
        .smiles {{
            background: #f5f5f5;
            padding: 10px;
            border-radius: 5px;
            font-family: monospace;
            word-break: break-all;
            font-size: 0.9em;
            margin: 10px 0;
            border-left: 4px solid #667eea;
        }}
        
        .properties {{
            margin: 15px 0;
        }}
        
        .properties h4 {{
            color: #667eea;
            margin-bottom: 8px;
            font-size: 0.95em;
        }}
        
        .prop-row {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 10px;
            margin-bottom: 8px;
        }}
        
        .prop-item {{
            background: #f9f9f9;
            padding: 8px;
            border-radius: 5px;
            font-size: 0.9em;
        }}
        
        .prop-label {{
            color: #666;
            font-weight: 600;
        }}
        
        .prop-value {{
            color: #333;
        }}
        
        .generation {{
            background: #e8f5e9;
            color: #2e7d32;
            padding: 5px 10px;
            border-radius: 5px;
            font-size: 0.85em;
            display: inline-block;
            margin-top: 10px;
        }}
        
        .timestamp {{
            color: #999;
            font-size: 0.85em;
            margin-top: 10px;
        }}
        
        .predictions {{
            background: #fff3e0;
            padding: 15px;
            border-radius: 5px;
            margin-top: 10px;
        }}
        
        .predictions h4 {{
            color: #e65100;
            margin-bottom: 8px;
        }}
        
        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 3px;
            font-size: 0.85em;
            margin-right: 5px;
            margin-bottom: 5px;
        }}
        
        .badge-success {{
            background: #c8e6c9;
            color: #2e7d32;
        }}
        
        .badge-warning {{
            background: #ffe0b2;
            color: #e65100;
        }}
        
        .badge-danger {{
            background: #ffcdd2;
            color: #c62828;
        }}
        
        .summary {{
            background: white;
            border-radius: 10px;
            padding: 20px;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }}
        
        .summary h2 {{
            color: #667eea;
            margin-bottom: 15px;
        }}
        
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
        }}
        
        .summary-table th {{
            background: #667eea;
            color: white;
            padding: 12px;
            text-align: left;
        }}
        
        .summary-table td {{
            padding: 12px;
            border-bottom: 1px solid #eee;
        }}
        
        .summary-table tr:hover {{
            background: #f5f5f5;
        }}
        
        .footer {{
            text-align: center;
            color: white;
            margin-top: 30px;
            padding: 20px;
        }}
        
        @media (max-width: 768px) {{
            .molecules-grid {{
                grid-template-columns: 1fr;
            }}
            
            .header h1 {{
                font-size: 1.8em;
            }}
            
            .prop-row {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 ChemAI Drug Discovery Dashboard</h1>
            <p>Real-time molecule creation and analysis</p>
            
            <div class="stats">
                <div class="stat-card">
                    <h3>{len(self.molecules)}</h3>
                    <p>Molecules Created</p>
                </div>
                <div class="stat-card">
                    <h3>{self.generation_count}</h3>
                    <p>Generations</p>
                </div>
                <div class="stat-card">
                    <h3>{self._get_avg_qed():.2f}</h3>
                    <p>Avg Drug-likeness</p>
                </div>
                <div class="stat-card">
                    <h3>{self._get_avg_potency():.2f}</h3>
                    <p>Avg Potency (pIC50)</p>
                </div>
            </div>
        </div>
        
        <h2 style="color: white; margin-bottom: 20px;">Discovered Molecules</h2>
        
        <div class="molecules-grid">
            {self._generate_molecule_cards()}
        </div>
        
        <div class="summary">
            <h2>Summary</h2>
            <table class="summary-table">
                <thead>
                    <tr>
                        <th>#</th>
                        <th>Molecule Name</th>
                        <th>SMILES</th>
                        <th>Generation</th>
                        <th>QED Score</th>
                        <th>pIC50</th>
                        <th>Created</th>
                    </tr>
                </thead>
                <tbody>
                    {self._generate_summary_rows()}
                </tbody>
            </table>
        </div>
        
        <div class="footer">
            <p>ChemAI Drug Discovery System • {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>Powered by RDKit, scikit-learn, and ChEMBL</p>
        </div>
    </div>
</body>
</html>
        """
        return html_content
    
    def _generate_molecule_cards(self):
        """Generate HTML cards for each molecule"""
        cards = []
        for mol in self.molecules:
            predictions = mol.get('predictions', {})
            
            qed = predictions.get('qed_score', 'N/A')
            pic50 = predictions.get('pIC50', 'N/A')
            druglike = predictions.get('is_druglike', 'Unknown')
            
            qed_badge = 'badge-success' if isinstance(qed, (int, float)) and qed > 0.6 else 'badge-warning'
            pic50_badge = 'badge-success' if isinstance(pic50, (int, float)) and pic50 > 6 else 'badge-warning'
            
            card = f"""
            <div class="molecule-card">
                <div class="molecule-id">#{mol['id']}</div>
                <div class="molecule-name">{mol['name']}</div>
                
                <div class="smiles">{mol['smiles']}</div>
                
                <div class="predictions">
                    <h4>Predictions</h4>
                    <div>
                        <span class="badge {qed_badge}">QED: {qed if isinstance(qed, str) else f'{qed:.3f}'}</span>
                        <span class="badge {pic50_badge}">pIC50: {pic50 if isinstance(pic50, str) else f'{pic50:.3f}'}</span>
                        <span class="badge {'badge-success' if druglike else 'badge-danger'}">Drug-like: {'Yes' if druglike else 'No'}</span>
                    </div>
                </div>
                
                <div class="generation">Generation {mol['generation']}</div>
                <div class="timestamp">{mol['timestamp']}</div>
            </div>
            """
            cards.append(card)
        
        return '\n'.join(cards) if cards else '<p style="grid-column: 1/-1; text-align: center; color: white;">No molecules created yet...</p>'
    
    def _generate_summary_rows(self):
        """Generate summary table rows"""
        rows = []
        for mol in self.molecules:
            predictions = mol.get('predictions', {})
            qed = predictions.get('qed_score', 'N/A')
            pic50 = predictions.get('pIC50', 'N/A')
            
            qed_str = f"{qed:.3f}" if isinstance(qed, (int, float)) else str(qed)
            pic50_str = f"{pic50:.3f}" if isinstance(pic50, (int, float)) else str(pic50)
            
            row = f"""
            <tr>
                <td>{mol['id']}</td>
                <td>{mol['name']}</td>
                <td style="font-family: monospace; font-size: 0.9em;">{mol['smiles'][:40]}...</td>
                <td>{mol['generation']}</td>
                <td>{qed_str}</td>
                <td>{pic50_str}</td>
                <td>{mol['timestamp'][:19]}</td>
            </tr>
            """
            rows.append(row)
        
        return '\n'.join(rows) if rows else '<tr><td colspan="7" style="text-align: center;">No molecules yet</td></tr>'
    
    def _get_avg_qed(self):
        """Calculate average QED score"""
        qeds = [m.get('predictions', {}).get('qed_score', 0) 
                for m in self.molecules if isinstance(m.get('predictions', {}).get('qed_score'), (int, float))]
        return sum(qeds) / len(qeds) if qeds else 0.0
    
    def _get_avg_potency(self):
        """Calculate average pIC50"""
        pic50s = [m.get('predictions', {}).get('pIC50', 0) 
                  for m in self.molecules if isinstance(m.get('predictions', {}).get('pIC50'), (int, float))]
        return sum(pic50s) / len(pic50s) if pic50s else 0.0
    
    def save_html(self):
        """Save HTML dashboard to file"""
        html = self.generate_html()
        with open(self.output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"\n📊 Dashboard saved: {self.output_file}")
        return self.output_file
    
    def open_in_browser(self):
        """Open dashboard in default browser"""
        import webbrowser
        file_path = Path(self.output_file).resolve()
        webbrowser.open(f'file://{file_path}')


# Global dashboard instance
dashboard = None

def init_dashboard(output_file='drug_discovery_dashboard.html'):
    """Initialize the global dashboard"""
    global dashboard
    dashboard = WebDashboard(output_file)
    return dashboard

def add_molecule_to_dashboard(smiles, name=None, properties=None, predictions=None, generation=1):
    """Add molecule to dashboard"""
    global dashboard
    if dashboard is None:
        dashboard = init_dashboard()
    return dashboard.add_molecule(smiles, name, properties, predictions, generation)

def finalize_dashboard():
    """Save and open dashboard"""
    global dashboard
    if dashboard is not None:
        dashboard.save_html()
        dashboard.open_in_browser()
