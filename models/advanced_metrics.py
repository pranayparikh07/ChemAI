"""
Advanced Metrics Calculator for ChemAI Models
Extended metrics beyond standard regression/classification
"""

import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score


class AdvancedMetrics:
    """Calculate extended performance metrics"""
    
    @staticmethod
    def calculate_regression_metrics(y_true, y_pred):
        """Calculate comprehensive regression metrics"""
        metrics = {}
        
        # Basic metrics
        metrics['rmse'] = np.sqrt(mean_squared_error(y_true, y_pred))
        metrics['mae'] = mean_absolute_error(y_true, y_pred)
        metrics['r2'] = r2_score(y_true, y_pred)
        
        # Additional metrics
        residuals = y_true - y_pred
        metrics['residuals_mean'] = np.mean(residuals)
        metrics['residuals_std'] = np.std(residuals)
        metrics['residuals_median'] = np.median(residuals)
        metrics['residuals_min'] = np.min(residuals)
        metrics['residuals_max'] = np.max(residuals)
        metrics['residuals_iqr'] = np.percentile(residuals, 75) - np.percentile(residuals, 25)
        
        # Error percentiles
        abs_errors = np.abs(residuals)
        metrics['error_p50'] = np.percentile(abs_errors, 50)  # Median error
        metrics['error_p90'] = np.percentile(abs_errors, 90)  # 90th percentile
        metrics['error_p95'] = np.percentile(abs_errors, 95)  # 95th percentile
        
        # Explained variance
        metrics['explained_variance'] = 1 - (np.var(residuals) / np.var(y_true))
        
        # Mean absolute percentage error (MAPE)
        mask = y_true != 0
        if np.any(mask):
            mape = np.mean(np.abs((y_true[mask] - y_pred[mask]) / y_true[mask])) * 100
            metrics['mape'] = mape
        
        # Mean squared log error
        if np.all(y_true > 0) and np.all(y_pred > 0):
            metrics['msle'] = mean_squared_error(np.log1p(y_true), np.log1p(y_pred))
        
        return metrics
    
    @staticmethod
    def calculate_classification_metrics(y_true, y_pred, y_pred_proba=None):
        """Calculate comprehensive classification metrics"""
        from sklearn.metrics import (
            confusion_matrix, accuracy_score, precision_score,
            recall_score, f1_score, roc_auc_score
        )
        
        metrics = {}
        
        # Basic metrics
        cm = confusion_matrix(y_true, y_pred)
        metrics['accuracy'] = accuracy_score(y_true, y_pred)
        
        # For binary classification
        if len(np.unique(y_true)) == 2:
            tn, fp, fn, tp = cm.ravel()
            metrics['true_positives'] = int(tp)
            metrics['true_negatives'] = int(tn)
            metrics['false_positives'] = int(fp)
            metrics['false_negatives'] = int(fn)
            
            # Sensitivity/Specificity
            metrics['sensitivity'] = tp / (tp + fn) if (tp + fn) > 0 else 0
            metrics['specificity'] = tn / (tn + fp) if (tn + fp) > 0 else 0
            
            # Precision/Recall
            metrics['precision'] = precision_score(y_true, y_pred, zero_division=0)
            metrics['recall'] = recall_score(y_true, y_pred, zero_division=0)
            metrics['f1_score'] = f1_score(y_true, y_pred, zero_division=0)
            
            # Matthews correlation coefficient
            denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
            if denominator != 0:
                metrics['mcc'] = (tp * tn - fp * fn) / denominator
            
            # ROC AUC
            if y_pred_proba is not None:
                try:
                    metrics['roc_auc'] = roc_auc_score(y_true, y_pred_proba)
                except:
                    pass
        
        return metrics
    
    @staticmethod
    def generate_metrics_report(y_true, y_pred, model_type='regression'):
        """Generate formatted metrics report"""
        report = {}
        
        if model_type == 'regression':
            metrics = AdvancedMetrics.calculate_regression_metrics(y_true, y_pred)
            
            report['Performance'] = {
                'R² Score': f"{metrics['r2']:.4f}",
                'RMSE': f"{metrics['rmse']:.4f}",
                'MAE': f"{metrics['mae']:.4f}",
            }
            
            report['Variance'] = {
                'Explained Variance': f"{metrics['explained_variance']:.4f}",
                'Residuals Std Dev': f"{metrics['residuals_std']:.4f}",
                'Residuals Mean': f"{metrics['residuals_mean']:.4f}",
            }
            
            report['Error Distribution'] = {
                'Median Error': f"{metrics['error_p50']:.4f}",
                '90th Percentile': f"{metrics['error_p90']:.4f}",
                '95th Percentile': f"{metrics['error_p95']:.4f}",
                'Min Residual': f"{metrics['residuals_min']:.4f}",
                'Max Residual': f"{metrics['residuals_max']:.4f}",
            }
            
            if 'mape' in metrics:
                report['Error Distribution']['MAPE (%)'] = f"{metrics['mape']:.2f}%"
        
        elif model_type == 'classification':
            metrics = AdvancedMetrics.calculate_classification_metrics(y_true, y_pred)
            
            report['Performance'] = {
                'Accuracy': f"{metrics['accuracy']:.4f}",
            }
            
            if 'sensitivity' in metrics:
                report['Classification'] = {
                    'Sensitivity': f"{metrics['sensitivity']:.4f}",
                    'Specificity': f"{metrics['specificity']:.4f}",
                    'Precision': f"{metrics['precision']:.4f}",
                    'Recall': f"{metrics['recall']:.4f}",
                    'F1 Score': f"{metrics['f1_score']:.4f}",
                }
                
                report['Confusion Matrix'] = {
                    'True Positives': str(metrics['true_positives']),
                    'True Negatives': str(metrics['true_negatives']),
                    'False Positives': str(metrics['false_positives']),
                    'False Negatives': str(metrics['false_negatives']),
                }
        
        return report
    
    @staticmethod
    def format_report_text(report):
        """Format report as text"""
        lines = []
        for section, metrics in report.items():
            lines.append(f"\n{section}:")
            lines.append("-" * (len(section) + 1))
            for key, value in metrics.items():
                lines.append(f"  {key:<30} {value}")
        return "\n".join(lines)


# Example usage
if __name__ == "__main__":
    # Example regression
    y_true = np.array([3.0, -0.5, 2.0, 7.0, 4.3])
    y_pred = np.array([2.5, 0.0, 2.0, 8.0, 4.0])
    
    print("Regression Metrics Example")
    print("=" * 50)
    report = AdvancedMetrics.generate_metrics_report(y_true, y_pred, 'regression')
    print(AdvancedMetrics.format_report_text(report))
    
    # Example classification
    print("\n\nClassification Metrics Example")
    print("=" * 50)
    y_true_class = np.array([0, 0, 1, 1, 1, 0, 1, 1])
    y_pred_class = np.array([0, 0, 1, 1, 0, 0, 1, 0])
    
    report_class = AdvancedMetrics.generate_metrics_report(y_true_class, y_pred_class, 'classification')
    print(AdvancedMetrics.format_report_text(report_class))
