#!/usr/bin/env python3
"""
SQLite Database Viewer - Web UI for ChEMBL database
Access at: http://localhost:5000
"""

from flask import Flask, render_template_string, jsonify, request
import sqlite3
import json
from pathlib import Path

app = Flask(__name__)
DB_PATH = "chembl_36/chembl_36_sqlite/chembl_36.db"

def get_connection():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn

HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ChEMBL Database Viewer</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif; background: #f5f5f5; }
        .container { max-width: 1400px; margin: 0 auto; padding: 20px; }
        header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 8px; margin-bottom: 30px; }
        header h1 { font-size: 2.5em; margin-bottom: 10px; }
        header p { font-size: 1.1em; opacity: 0.9; }
        
        .tables-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 15px; margin-bottom: 30px; }
        .table-card { background: white; border-radius: 8px; padding: 20px; cursor: pointer; transition: all 0.3s; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
        .table-card:hover { transform: translateY(-4px); box-shadow: 0 4px 16px rgba(0,0,0,0.15); background: #667eea; color: white; }
        .table-card h3 { margin-bottom: 10px; font-size: 1.2em; }
        .table-card .row-count { font-size: 0.9em; opacity: 0.7; }
        .table-card .col-count { font-size: 0.85em; margin-top: 5px; opacity: 0.7; }
        
        .data-view { background: white; border-radius: 8px; padding: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
        .data-view h2 { margin-bottom: 20px; color: #333; }
        .close-btn { float: right; cursor: pointer; font-size: 1.5em; color: #666; }
        .close-btn:hover { color: #000; }
        
        table { width: 100%; border-collapse: collapse; }
        th { background: #f8f9fa; padding: 12px; text-align: left; font-weight: 600; border-bottom: 2px solid #dee2e6; }
        td { padding: 12px; border-bottom: 1px solid #dee2e6; }
        tr:hover { background: #f8f9fa; }
        
        .pagination { margin-top: 20px; text-align: center; }
        .pagination button { background: #667eea; color: white; border: none; padding: 8px 15px; margin: 0 5px; border-radius: 4px; cursor: pointer; }
        .pagination button:hover { background: #764ba2; }
        .pagination button:disabled { background: #ccc; cursor: not-allowed; }
        .pagination span { margin: 0 10px; }
        
        .loading { text-align: center; padding: 20px; color: #666; }
        .spinner { border: 4px solid #f3f3f3; border-top: 4px solid #667eea; border-radius: 50%; width: 40px; height: 40px; animation: spin 1s linear infinite; margin: 0 auto 10px; }
        @keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }
        
        .error { background: #f8d7da; color: #721c24; padding: 15px; border-radius: 4px; margin-bottom: 20px; }
        .success { background: #d4edda; color: #155724; padding: 15px; border-radius: 4px; margin-bottom: 20px; }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>📊 ChEMBL v36 Database Viewer</h1>
            <p>Interactive SQLite Database Explorer - Click any table to view data</p>
        </header>
        
        <div id="tablesContainer" class="tables-grid"></div>
        <div id="dataView"></div>
    </div>

    <script>
        const API_URL = '/api';
        let currentPage = 1;
        let currentTable = null;
        let currentLimit = 50;

        async function loadTables() {
            try {
                const response = await fetch(`${API_URL}/tables`);
                const tables = await response.json();
                
                const container = document.getElementById('tablesContainer');
                container.innerHTML = '';
                
                tables.forEach(table => {
                    const card = document.createElement('div');
                    card.className = 'table-card';
                    card.innerHTML = `
                        <h3>${table.name}</h3>
                        <div class="row-count">📈 ${table.rows.toLocaleString()} rows</div>
                        <div class="col-count">📋 ${table.columns} columns</div>
                    `;
                    card.onclick = () => loadTableData(table.name);
                    container.appendChild(card);
                });
            } catch (error) {
                console.error('Error loading tables:', error);
            }
        }

        async function loadTableData(tableName, page = 1, limit = 50) {
            currentTable = tableName;
            currentPage = page;
            
            const dataView = document.getElementById('dataView');
            dataView.innerHTML = `
                <div class="loading">
                    <div class="spinner"></div>
                    Loading data from ${tableName}...
                </div>
            `;
            
            try {
                const response = await fetch(`${API_URL}/table/${tableName}?page=${page}&limit=${limit}`);
                const data = await response.json();
                
                let html = `
                    <div class="data-view">
                        <h2>
                            ${tableName}
                            <span class="close-btn" onclick="document.getElementById('dataView').innerHTML=''">✕</span>
                        </h2>
                        <p style="color: #666; margin-bottom: 15px;">
                            Total rows: <strong>${data.total.toLocaleString()}</strong> | 
                            Columns: <strong>${data.columns.length}</strong>
                        </p>
                        
                        <div style="overflow-x: auto;">
                            <table>
                                <thead>
                                    <tr>
                                        ${data.columns.map(col => `<th>${col}</th>`).join('')}
                                    </tr>
                                </thead>
                                <tbody>
                                    ${data.rows.map(row => `
                                        <tr>
                                            ${data.columns.map(col => {
                                                const value = row[col];
                                                let displayValue = value === null ? '<em style="color: #999;">NULL</em>' : String(value).substring(0, 100);
                                                return `<td>${displayValue}</td>`;
                                            }).join('')}
                                        </tr>
                                    `).join('')}
                                </tbody>
                            </table>
                        </div>
                        
                        <div class="pagination">
                            <button onclick="loadTableData('${tableName}', ${Math.max(1, page - 1)}, ${limit})" ${page === 1 ? 'disabled' : ''}>← Previous</button>
                            <span>Page ${page} of ${Math.ceil(data.total / limit)}</span>
                            <button onclick="loadTableData('${tableName}', ${page + 1}, ${limit})" ${page * limit >= data.total ? 'disabled' : ''}>Next →</button>
                            
                            <select onchange="loadTableData('${tableName}', 1, this.value)" style="margin-left: 20px; padding: 6px;">
                                <option value="10" ${limit === 10 ? 'selected' : ''}>Show 10</option>
                                <option value="25" ${limit === 25 ? 'selected' : ''}>Show 25</option>
                                <option value="50" ${limit === 50 ? 'selected' : ''}>Show 50</option>
                                <option value="100" ${limit === 100 ? 'selected' : ''}>Show 100</option>
                            </select>
                        </div>
                    </div>
                `;
                
                dataView.innerHTML = html;
            } catch (error) {
                dataView.innerHTML = `<div class="error">Error loading table data: ${error.message}</div>`;
            }
        }

        // Load tables on page load
        window.onload = loadTables;
    </script>
</body>
</html>
'''

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/api/tables')
def api_tables():
    """Get all tables with metadata"""
    try:
        conn = get_connection()
        cursor = conn.cursor()
        
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
        tables = cursor.fetchall()
        
        result = []
        for table in tables:
            table_name = table[0]
            
            # Get row count
            cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
            row_count = cursor.fetchone()[0]
            
            # Get column count
            cursor.execute(f"PRAGMA table_info({table_name})")
            columns = cursor.fetchall()
            
            result.append({
                'name': table_name,
                'rows': row_count,
                'columns': len(columns)
            })
        
        conn.close()
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/table/<table_name>')
def api_table_data(table_name):
    """Get table data with pagination"""
    try:
        page = int(request.args.get('page', 1))
        limit = int(request.args.get('limit', 50))
        
        conn = get_connection()
        cursor = conn.cursor()
        
        # Get total count
        cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
        total = cursor.fetchone()[0]
        
        # Get columns
        cursor.execute(f"PRAGMA table_info({table_name})")
        columns = [col[1] for col in cursor.fetchall()]
        
        # Get paginated data
        offset = (page - 1) * limit
        cursor.execute(f"SELECT * FROM {table_name} LIMIT {limit} OFFSET {offset}")
        rows = cursor.fetchall()
        
        # Convert rows to dictionaries
        rows_data = [dict(row) for row in rows]
        
        conn.close()
        
        return jsonify({
            'total': total,
            'columns': columns,
            'rows': rows_data
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    print("=" * 80)
    print("🚀 ChEMBL Database Viewer")
    print("=" * 80)
    print()
    print("📍 Server running at: http://localhost:5000")
    print("📊 Database: chembl_36/chembl_36_sqlite/chembl_36.db")
    print()
    print("Press Ctrl+C to stop the server")
    print("=" * 80)
    print()
    
    app.run(debug=False, host='localhost', port=5000)
