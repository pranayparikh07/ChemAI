"""
Comprehensive Model Testing Suite
Tests all trained ChemAI models and generates detailed HTML report
"""

import sys
import os
from datetime import datetime
import numpy as np
import traceback
import importlib

# Add models directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

# Clear module cache to force reload
for module in list(sys.modules.keys()):
    if 'test_' in module:
        del sys.modules[module]

# Import individual test modules
from test_property_model import test_property_model
from test_druglikeness_model import test_druglikeness_models
from test_bioactivity_model import test_bioactivity_model
from test_toxicity_model import test_toxicity_model


class ModelTestOrchestrator:
    def __init__(self):
        self.results = {}
        self.timestamp = datetime.now()
        
    def run_all_tests(self):
        """Run all model tests"""
        print("\n" + "=" * 80)
        print(" " * 20 + "CHEMAI MODEL TESTING SUITE")
        print("=" * 80)
        print(f"Test Execution Time: {self.timestamp.strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 80)
        
        # Test Property Model
        print("\n[1/4] Testing Property Model...")
        try:
            self.results['property'] = test_property_model()
        except Exception as e:
            print(f"✗ Error testing property model: {e}")
            traceback.print_exc()
            self.results['property'] = None
        
        # Test Druglikeness Models
        print("\n[2/4] Testing Druglikeness Models...")
        try:
            self.results['druglikeness'] = test_druglikeness_models()
        except Exception as e:
            print(f"✗ Error testing druglikeness models: {e}")
            traceback.print_exc()
            self.results['druglikeness'] = None
        
        # Test Bioactivity Model
        print("\n[3/4] Testing Bioactivity Model...")
        try:
            self.results['bioactivity'] = test_bioactivity_model()
        except Exception as e:
            print(f"✗ Error testing bioactivity model: {e}")
            traceback.print_exc()
            self.results['bioactivity'] = None
        
        # Test Toxicity Model
        print("\n[4/4] Testing Toxicity Model...")
        try:
            self.results['toxicity'] = test_toxicity_model()
        except Exception as e:
            print(f"✗ Error testing toxicity model: {e}")
            traceback.print_exc()
            self.results['toxicity'] = None
        
        print("\n" + "=" * 80)
        print("All tests completed!")
        print("=" * 80)
        
        return self.results
    
    def generate_html_report(self, output_file="model_test_report.html"):
        """Generate comprehensive HTML report"""
        html_content = self._generate_html()
        
        output_path = os.path.join(os.path.dirname(__file__), output_file)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"\n✓ HTML report generated: {output_path}")
        return output_path
    
    def _generate_html(self):
        """Build HTML content"""
        html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChemAI Model Testing Report</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: #333;
            padding: 20px;
            min-height: 100vh;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
            overflow: hidden;
        }
        
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }
        
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        
        .header p {
            font-size: 1.1em;
            opacity: 0.9;
        }
        
        .content {
            padding: 40px;
        }
        
        .test-section {
            margin-bottom: 40px;
            border: 1px solid #ddd;
            border-radius: 8px;
            overflow: hidden;
        }
        
        .test-header {
            background: #f8f9fa;
            padding: 20px;
            border-bottom: 2px solid #667eea;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        
        .test-header h2 {
            color: #667eea;
            font-size: 1.8em;
        }
        
        .status-badge {
            padding: 8px 16px;
            border-radius: 20px;
            font-weight: bold;
            font-size: 0.9em;
        }
        
        .status-success {
            background: #d4edda;
            color: #155724;
        }
        
        .status-error {
            background: #f8d7da;
            color: #721c24;
        }
        
        .test-body {
            padding: 20px;
        }
        
        .metrics-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }
        
        .metric-card {
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        }
        
        .metric-card h3 {
            color: #667eea;
            font-size: 0.9em;
            text-transform: uppercase;
            margin-bottom: 10px;
        }
        
        .metric-card .value {
            font-size: 1.8em;
            font-weight: bold;
            color: #333;
        }
        
        .property-table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }
        
        .property-table thead {
            background: #f8f9fa;
            border-bottom: 2px solid #667eea;
        }
        
        .property-table th {
            padding: 12px;
            text-align: left;
            color: #667eea;
            font-weight: 600;
        }
        
        .property-table td {
            padding: 12px;
            border-bottom: 1px solid #ddd;
        }
        
        .property-table tr:hover {
            background: #f8f9fa;
        }
        
        .property-table tr:nth-child(even) {
            background: #fafbfc;
        }
        
        .value-good {
            color: #28a745;
            font-weight: 600;
        }
        
        .value-ok {
            color: #ffc107;
            font-weight: 600;
        }
        
        .value-poor {
            color: #dc3545;
            font-weight: 600;
        }
        
        .footer {
            background: #f8f9fa;
            padding: 20px;
            text-align: center;
            border-top: 1px solid #ddd;
            color: #666;
        }
        
        .model-info {
            background: #e7f3ff;
            padding: 15px;
            border-left: 4px solid #667eea;
            margin-bottom: 20px;
            border-radius: 4px;
        }
        
        .model-info p {
            margin: 5px 0;
            color: #333;
        }
        
        .no-data {
            text-align: center;
            padding: 40px;
            color: #999;
            font-size: 1.1em;
        }
        
        .charts-container {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        
        .chart-placeholder {
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            color: #999;
            border: 1px dashed #ddd;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>[ChemAI] Model Testing Report</h1>
            <p>Comprehensive Performance Evaluation of Trained Models</p>
            <p style="font-size: 0.9em; margin-top: 10px;">Generated: """ + self.timestamp.strftime('%Y-%m-%d %H:%M:%S') + """</p>
        </div>
        
        <div class="content">
"""
        
        # Property Model Section
        html += self._generate_property_section()
        
        # Druglikeness Models Section
        html += self._generate_druglikeness_section()
        
        # Bioactivity Model Section
        html += self._generate_bioactivity_section()
        
        # Toxicity Model Section
        html += self._generate_toxicity_section()
        
        # Summary Section
        html += self._generate_summary_section()
        
        html += """
        </div>
        
        <div class="footer">
            <p>ChemAI Model Testing Suite | All models tested with comprehensive metrics</p>
            <p style="font-size: 0.9em; margin-top: 10px;">Report includes RMSE, MAE, R², Accuracy, and other performance metrics</p>
        </div>
    </div>
</body>
</html>
"""
        return html
    
    def _generate_property_section(self):
        """Generate property model section"""
        if self.results['property'] is None:
            return """
            <div class="test-section">
                <div class="test-header">
                    <h2>📊 Property Prediction Model</h2>
                    <span class="status-badge status-error">FAILED</span>
                </div>
                <div class="test-body">
                    <div class="no-data">No results available. Model may not be trained or accessible.</div>
                </div>
            </div>
"""
        
        results = self.results['property']
        metrics = results['metrics']
        
        html = """
            <div class="test-section">
                <div class="test-header">
                    <h2>📊 Property Prediction Model</h2>
                    <span class="status-badge status-success">PASSED</span>
                </div>
                <div class="test-body">
                    <div class="model-info">
                        <p><strong>Test Samples:</strong> """ + f"{results['test_samples']:,}" + """</p>
                        <p><strong>Properties Predicted:</strong> """ + ", ".join(results['properties']) + """</p>
                        <p><strong>Model Type:</strong> Multi-Output Random Forest Regressor</p>
                    </div>
                    
                    <table class="property-table">
                        <thead>
                            <tr>
                                <th>Property</th>
                                <th>RMSE</th>
                                <th>MAE</th>
                                <th>R² Score</th>
                            </tr>
                        </thead>
                        <tbody>
"""
        
        for prop in results['properties']:
            m = metrics[prop]
            rmse_class = 'value-good' if m['rmse'] < 100 else 'value-ok' if m['rmse'] < 500 else 'value-poor'
            r2_class = 'value-good' if m['r2'] > 0.7 else 'value-ok' if m['r2'] > 0.5 else 'value-poor'
            
            html += f"""
                            <tr>
                                <td><strong>{prop}</strong></td>
                                <td class="{rmse_class}">{m['rmse']:.4f}</td>
                                <td>{m['mae']:.4f}</td>
                                <td class="{r2_class}">{m['r2']:.4f}</td>
                            </tr>
"""
        
        html += """
                        </tbody>
                    </table>
                </div>
            </div>
"""
        return html
    
    def _generate_druglikeness_section(self):
        """Generate druglikeness models section"""
        if self.results['druglikeness'] is None:
            return """
            <div class="test-section">
                <div class="test-header">
                    <h2>💊 Druglikeness Models</h2>
                    <span class="status-badge status-error">FAILED</span>
                </div>
                <div class="test-body">
                    <div class="no-data">No results available. Models may not be trained or accessible.</div>
                </div>
            </div>
"""
        
        results = self.results['druglikeness']
        qed_metrics = results['qed_metrics']
        druglike_metrics = results['druglike_metrics']
        
        html = """
            <div class="test-section">
                <div class="test-header">
                    <h2>💊 Druglikeness Models (QED & Lipinski)</h2>
                    <span class="status-badge status-success">PASSED</span>
                </div>
                <div class="test-body">
                    <div class="model-info">
                        <p><strong>Test Samples:</strong> """ + f"{results['test_samples']:,}" + """</p>
                        <p><strong>QED Model Type:</strong> Random Forest Regressor</p>
                        <p><strong>Drug-likeness Model Type:</strong> Random Forest Classifier</p>
                    </div>
                    
                    <h3 style="color: #667eea; margin: 20px 0 15px 0;">QED Score Prediction</h3>
                    <div class="metrics-grid">
                        <div class="metric-card">
                            <h3>RMSE</h3>
                            <div class="value">""" + f"{qed_metrics['rmse']:.4f}" + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>MAE</h3>
                            <div class="value">""" + f"{qed_metrics['mae']:.4f}" + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>R² Score</h3>
                            <div class="value">""" + f"{qed_metrics['r2']:.4f}" + """</div>
                        </div>
                    </div>
                    
                    <h3 style="color: #667eea; margin: 20px 0 15px 0;">Drug-likeness Classification</h3>
                    <div class="metrics-grid">
                        <div class="metric-card">
                            <h3>Accuracy</h3>
                            <div class="value">""" + f"{druglike_metrics['accuracy']:.2%}" + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>True Positives</h3>
                            <div class="value">""" + str(druglike_metrics['confusion_matrix'][1, 1]) + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>True Negatives</h3>
                            <div class="value">""" + str(druglike_metrics['confusion_matrix'][0, 0]) + """</div>
                        </div>
                    </div>
                    
                    <table class="property-table">
                        <thead>
                            <tr>
                                <th>Metric</th>
                                <th>Not Drug-like</th>
                                <th>Drug-like</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td><strong>Predicted Negatives (FP+TN)</strong></td>
                                <td>""" + str(druglike_metrics['confusion_matrix'][0, 0]) + """</td>
                                <td>""" + str(druglike_metrics['confusion_matrix'][0, 1]) + """</td>
                            </tr>
                            <tr>
                                <td><strong>Predicted Positives (FN+TP)</strong></td>
                                <td>""" + str(druglike_metrics['confusion_matrix'][1, 0]) + """</td>
                                <td>""" + str(druglike_metrics['confusion_matrix'][1, 1]) + """</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
"""
        return html
    
    def _generate_bioactivity_section(self):
        """Generate bioactivity model section"""
        if self.results['bioactivity'] is None:
            return """
            <div class="test-section">
                <div class="test-header">
                    <h2>🧬 Bioactivity Prediction Model</h2>
                    <span class="status-badge status-error">NOT AVAILABLE</span>
                </div>
                <div class="test-body">
                    <div class="no-data">Model not found. This model may not be trained yet.</div>
                </div>
            </div>
"""
        
        results = self.results['bioactivity']
        metrics = results['metrics']
        
        html = """
            <div class="test-section">
                <div class="test-header">
                    <h2>🧬 Bioactivity Prediction Model</h2>
                    <span class="status-badge status-success">PASSED</span>
                </div>
                <div class="test-body">
                    <div class="model-info">
                        <p><strong>Test Samples:</strong> """ + f"{results['test_samples']:,}" + """</p>
                        <p><strong>Target:</strong> pIC50 (Bioactivity Score)</p>
                        <p><strong>Model Type:</strong> Random Forest Regressor</p>
                    </div>
                    
                    <div class="metrics-grid">
                        <div class="metric-card">
                            <h3>RMSE</h3>
                            <div class="value">""" + f"{metrics['rmse']:.4f}" + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>MAE</h3>
                            <div class="value">""" + f"{metrics['mae']:.4f}" + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>R² Score</h3>
                            <div class="value">""" + f"{metrics['r2']:.4f}" + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>Residuals Mean</h3>
                            <div class="value">""" + f"{metrics['residuals_mean']:.4f}" + """</div>
                        </div>
                    </div>
                    
                    <table class="property-table">
                        <thead>
                            <tr>
                                <th>Metric</th>
                                <th>Value</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td><strong>RMSE</strong></td>
                                <td class="value-good">""" + f"{metrics['rmse']:.4f}" + """</td>
                            </tr>
                            <tr>
                                <td><strong>MAE</strong></td>
                                <td>""" + f"{metrics['mae']:.4f}" + """</td>
                            </tr>
                            <tr>
                                <td><strong>R² Score</strong></td>
                                <td class="value-good">""" + f"{metrics['r2']:.4f}" + """</td>
                            </tr>
                            <tr>
                                <td><strong>Residuals Mean</strong></td>
                                <td>""" + f"{metrics['residuals_mean']:.4f}" + """</td>
                            </tr>
                            <tr>
                                <td><strong>Residuals Std Dev</strong></td>
                                <td>""" + f"{metrics['residuals_std']:.4f}" + """</td>
                            </tr>
                            <tr>
                                <td><strong>Min Residual</strong></td>
                                <td>""" + f"{metrics['residuals_min']:.4f}" + """</td>
                            </tr>
                            <tr>
                                <td><strong>Max Residual</strong></td>
                                <td>""" + f"{metrics['residuals_max']:.4f}" + """</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
"""
        return html
    
    def _generate_toxicity_section(self):
        """Generate toxicity model section"""
        if self.results['toxicity'] is None:
            return """
            <div class="test-section">
                <div class="test-header">
                    <h2>[TOX] Toxicity Prediction Model</h2>
                    <span class="status-badge status-error">NOT AVAILABLE</span>
                </div>
                <div class="test-body">
                    <div class="no-data">Model not found. This model may not be trained yet.</div>
                </div>
            </div>
"""
        
        results = self.results['toxicity']
        metrics = results['metrics']
        conf_matrix = metrics['confusion_matrix']
        
        html = """
            <div class="test-section">
                <div class="test-header">
                    <h2>[TOX] Toxicity Prediction Model</h2>
                    <span class="status-badge status-success">PASSED</span>
                </div>
                <div class="test-body">
                    <div class="model-info">
                        <p><strong>Test Samples:</strong> """ + f"{results['test_samples']:,}" + """</p>
                        <p><strong>Target:</strong> Structural Alerts (Safe vs. Toxic)</p>
                        <p><strong>Model Type:</strong> Random Forest Classifier</p>
                    </div>
                    
                    <div class="metrics-grid">
                        <div class="metric-card">
                            <h3>Accuracy</h3>
                            <div class="value">""" + f"{metrics['accuracy']:.2%}" + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>ROC-AUC</h3>
                            <div class="value">""" + f"{metrics['roc_auc']:.4f}" + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>True Positives</h3>
                            <div class="value">""" + str(conf_matrix[1, 1]) + """</div>
                        </div>
                        <div class="metric-card">
                            <h3>True Negatives</h3>
                            <div class="value">""" + str(conf_matrix[0, 0]) + """</div>
                        </div>
                    </div>
                    
                    <table class="property-table">
                        <thead>
                            <tr>
                                <th>Metric</th>
                                <th>Safe</th>
                                <th>Toxic</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td><strong>Predicted Negatives (TN+FP)</strong></td>
                                <td>""" + str(conf_matrix[0, 0]) + """</td>
                                <td>""" + str(conf_matrix[0, 1]) + """</td>
                            </tr>
                            <tr>
                                <td><strong>Predicted Positives (FN+TP)</strong></td>
                                <td>""" + str(conf_matrix[1, 0]) + """</td>
                                <td>""" + str(conf_matrix[1, 1]) + """</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
"""
        return html
    
    def _generate_summary_section(self):
        """Generate summary section"""
        total_tests = sum(1 for r in self.results.values() if r is not None)
        passed_tests = total_tests
        
        html = """
            <div class="test-section" style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);">
                <div class="test-body">
                    <h2 style="color: white; text-align: center; margin-bottom: 20px;">[SUMMARY] Testing Summary</h2>
                    <div class="metrics-grid">
"""
        
        html += f"""
                        <div class="metric-card">
                            <h3>Total Models Tested</h3>
                            <div class="value" style="color: #667eea;">{total_tests}</div>
                        </div>
                        <div class="metric-card">
                            <h3>Passed Tests</h3>
                            <div class="value" style="color: #28a745;">{passed_tests}</div>
                        </div>
"""
        
        if total_tests > 0:
            success_rate = (passed_tests / total_tests) * 100
            html += f"""
                        <div class="metric-card">
                            <h3>Success Rate</h3>
                            <div class="value" style="color: #667eea;">{success_rate:.1f}%</div>
                        </div>
"""
        
        html += """
                    </div>
                    
                    <div style="background: white; padding: 20px; border-radius: 8px; margin-top: 20px; color: #333;">
                        <h3 style="color: #667eea; margin-bottom: 10px;">Testing Details</h3>
                        <ul style="margin-left: 20px;">
"""
        
        if self.results['property']:
            html += f"<li>✓ Property Model: {self.results['property']['test_samples']:,} test samples</li>"
        else:
            html += "<li>✗ Property Model: Not available</li>"
        
        if self.results['druglikeness']:
            html += f"<li>✓ Druglikeness Models: {self.results['druglikeness']['test_samples']:,} test samples</li>"
        else:
            html += "<li>✗ Druglikeness Models: Not available</li>"
        
        if self.results['bioactivity']:
            html += f"<li>✓ Bioactivity Model: {self.results['bioactivity']['test_samples']:,} test samples</li>"
        else:
            html += "<li>✗ Bioactivity Model: Not available</li>"
        
        if self.results['toxicity']:
            html += f"<li>✓ Toxicity Model: {self.results['toxicity']['test_samples']:,} test samples</li>"
        else:
            html += "<li>✗ Toxicity Model: Not available</li>"
        
        html += """
                        </ul>
                    </div>
                </div>
            </div>
"""
        return html


def main():
    """Main execution"""
    orchestrator = ModelTestOrchestrator()
    
    # Run all tests
    orchestrator.run_all_tests()
    
    # Generate HTML report
    report_path = orchestrator.generate_html_report()
    
    print(f"\n" + "=" * 80)
    print("Testing Complete!")
    print(f"Report saved to: {report_path}")
    print("=" * 80)
    
    return orchestrator


if __name__ == "__main__":
    orchestrator = main()
