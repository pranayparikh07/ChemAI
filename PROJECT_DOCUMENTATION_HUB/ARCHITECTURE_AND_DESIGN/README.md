# 🏛️ ARCHITECTURE AND DESIGN

**System Design, Technical Roadmap, and Data Flow**

---

## 📂 Files in This Folder

| File | Purpose |
|------|---------|
| `TECHNICAL_ROADMAP.md` | Long-term technical plan |
| `PROJECT_PROGRESS_ANALYSIS.md` | Implementation status |
| `IMPLEMENTATION_SUMMARY.txt` | What's been built |
| `VISUAL_FLOW_DIAGRAM.txt` | System flow visualization |
| `COMPREHENSIVE_MONTH1_REPORT.txt` | Month 1 progress report |

---

## 🎯 Architecture Overview

```
┌────────────────────────────────────────────────────┐
│           ChemAI Multi-Agent System                │
└────────────────────────────────────────────────────┘

INPUT: Drug discovery requirements
   ↓
┌─ Generator Agent (GNN, Pranay) ─────────────────┐
│ • Creates novel molecules                        │
│ • Applies chemical constraints                   │
│ • Generates SMILES strings                       │
└──────────────┬────────────────────────────────────┘
               ↓
┌─ Predictor Agent ─────────────────────────────────┐
│ • Evaluates: bioactivity, toxicity, druglikeness│
│ • Returns: multi-objective scores                │
└──────────────┬────────────────────────────────────┘
               ↓
┌─ Optimizer Agent (RL, Vishwa) ────────────────────┐
│ • Modifies molecules via RL                      │
│ • Optimizes for multiple objectives              │
│ • Maximizes reward function                      │
└──────────────┬────────────────────────────────────┘
               ↓
┌─ Ranker Agent ────────────────────────────────────┐
│ • Ranks candidates                               │
│ • Selects top molecules                          │
│ • Filters for diversity                          │
└──────────────┬────────────────────────────────────┘
               ↓
┌─ Orchestrator Agent (Control, Shreya) ────────────┐
│ • Decides next action (generate/optimize/rank)  │
│ • Manages autonomous loop                        │
│ • Logs all decisions                             │
└──────────────┬────────────────────────────────────┘
               ↓
         Neo4j Knowledge Graph
     (Molecule properties & relationships)
```

---

## 📊 Data Flow

```
Experiment Initialization
  ↓
├─→ Load seed molecules
├─→ Setup model predictors
├─→ Initialize RL agent
└─→ Configure parameters
  ↓
┌─ GENERATION PHASE ────────────┐
│ Pranay's GNN generates mols   │
│ Applies valence rules         │
│ Creates novel structures      │
└───────────┬──────────────────┘
            ↓
┌─ EVALUATION PHASE ────────────┐
│ Predictor scores molecules    │
│ Multi-objective scoring       │
│ Stores results in Neo4j       │
└───────────┬──────────────────┘
            ↓
┌─ OPTIMIZATION PHASE ─────────┐
│ Vishwa's RL agent optimizes   │
│ Applies reward function       │
│ Generates improved molecules  │
└───────────┬──────────────────┘
            ↓
┌─ RANKING PHASE ──────────────┐
│ Ranker sorts by score        │
│ Filters diversity            │
│ Selects top candidates       │
└───────────┬──────────────────┘
            ↓
┌─ DECISION PHASE ─────────────┐
│ Shreya's control agent       │
│ Decides: Continue? Restart?  │
│ Logs causal analysis         │
└───────────┬──────────────────┘
            ↓
         Output Results
    (Best molecules found)
```

---

## 🔧 Technical Components

### **1. Generator Agent (Pranay)**
- **Input**: Seed molecules
- **Process**: GNN-based generation with mutations
- **Output**: Novel SMILES strings
- **Constraints**: 
  - Valence rules
  - Ring closures
  - Chemical validity
- **Metrics**: Validity ratio, Novelty, Diversity

### **2. Predictor Agent**
- **Input**: SMILES strings
- **Models**: 
  - Bioactivity predictor
  - Toxicity predictor
  - Drug-likeness scorer
  - Property predictor
- **Output**: Multi-objective scores
- **Location**: `models/train_*.py`

### **3. Optimizer Agent (Vishwa)**
- **Input**: Generated molecules + scores
- **RL Framework**:
  - State: Molecular representation
  - Action: Molecular modifications
  - Reward: Multi-objective function
- **Output**: Optimized molecules
- **Comparison**: Random search vs RL

### **4. Ranker Agent**
- **Input**: All candidates with scores
- **Process**:
  - Sort by multi-objective score
  - Filter for diversity
  - Select top-k candidates
- **Output**: Ranked candidate list

### **5. Orchestrator Agent (Shreya)**
- **Input**: Experiment config + all agent outputs
- **Logic**:
  - Coordinates agent execution
  - Makes autonomous decisions
  - Manages experiment loop
  - Logs causal traces
- **Output**: Final candidates + experiment logs

---

## 💾 Database Schema (Neo4j)

### **Nodes**
```
:Molecule
  - id (unique)
  - smiles
  - name
  - mw (molecular weight)
  - logp
  - qed (drug-likeness)
  - is_drug_like
  - source
  - confidence

:Protein
  - id (unique)
  - uniprot
  - name
  - gene
  - organism
  - function
  - source

:Experiment
  - id
  - timestamp
  - parameters
  - status
```

### **Relationships**
```
(Molecule)-[:TARGETS]->(Protein)
  - pic50
  - activity_type
  - confidence

(Molecule)-[:SIMILAR_TO]->(Molecule)
  - similarity
  - method (Tanimoto)

(Molecule)-[:GENERATED_IN]->(Experiment)
  - generation_number
  - score

(Molecule)-[:OPTIMIZED_TO]->(Molecule)
  - optimization_step
  - reward
```

---

## 🔄 Autonomous Loop Logic

```python
while experiment_running:
    # Decision making (Shreya)
    decision = orchestrator.make_decision(
        previous_scores=scores,
        generation_count=gen_count,
        convergence_check=check_convergence()
    )
    
    if decision == "GENERATE":
        # Generation (Pranay)
        molecules = generator.generate(
            num_molecules=100,
            seed=current_best
        )
    
    elif decision == "OPTIMIZE":
        # Optimization (Vishwa)
        molecules = optimizer.optimize(
            molecules=candidates,
            iterations=30,
            reward_fn=vishwa_reward
        )
    
    elif decision == "RANK":
        # Ranking
        ranked = ranker.rank(molecules, top_k=10)
    
    # Evaluation
    scores = predictor.evaluate(molecules)
    
    # Logging (Shreya)
    orchestrator.log_decision(
        decision=decision,
        input=molecules,
        output=scores,
        timestamp=now(),
        reason=decision_reason
    )
    
    # Check termination
    if should_stop():
        break
```

---

## 📈 Performance Metrics

### **Generator Metrics**
- Validity ratio (% valid molecules)
- Novelty (% new structures)
- Diversity (tanimoto similarity spread)

### **Predictor Metrics**
- R² score (per property)
- RMSE (per model)
- Correlation with ChEMBL

### **Optimizer Metrics**
- Improvement per iteration
- Reward convergence
- Comparison: Random vs RL

### **Overall Metrics**
- Candidates found
- Average score improvement
- Time per generation
- Total experiment time

---

## 🏗️ Implementation Checklist

- [x] Base agent framework created
- [x] Predictor models trained
- [x] Database schema designed
- [ ] **Shreya**: Orchestrator logic
- [ ] **Vishwa**: RL optimizer
- [ ] **Pranay**: GNN generator
- [ ] Full integration & testing
- [ ] Experiment execution & reporting

---

## 🚀 Next Phases

### **Phase 1** (Current)
- Implement core orchestrator
- Integrate individual components
- Test autonomous loop

### **Phase 2**
- Optimize reward function
- Implement advanced RL algorithms
- Scale to larger molecule sets

### **Phase 3**
- Deploy to cloud
- Add REST API
- Create web dashboard

---

## 📖 Related Documentation

- See `../../AGENTS_FRAMEWORK/` for agent details
- See `../../ML_MODELS/` for model specifics
- See `../../DATABASE_AND_GRAPHS/` for Neo4j schema

---

**Master Index**: Go back to `../../INDEX.md`
