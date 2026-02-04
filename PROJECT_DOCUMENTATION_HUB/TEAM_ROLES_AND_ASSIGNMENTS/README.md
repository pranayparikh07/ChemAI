# 👥 TEAM ROLES AND ASSIGNMENTS

**Team Structure, Responsibilities, and Deliverables**

---

## 📂 Files in This Folder

| File | Role/Purpose |
|------|-------------|
| `README_TEAM_ORGANIZATION.txt` | Overall team structure |
| `TEAM_STATUS_AND_REMAINING_WORK.txt` | Current progress & TODOs |
| `SHREYA_KNOWLEDGE_GRAPH_COMPLETION.md` | Shreya's specific deliverables |
| `TEAM_SHREYA/` | Shreya's code folder |
| `TEAM_VISHWA/` | Vishwa's code folder |
| `TEAM_PRANAY/` | Pranay's code folder |

---

## 👨‍💼 Team Breakdown

### **🔵 Shreya – Agent Orchestration & Autonomous Experiment Loop**

**Responsibility**: 
- Coordinate all agents (generator, optimizer, predictor, ranker)
- Implement autonomous decision-making logic
- Create experiment tracking & logging
- Manage Neo4j knowledge graph

**Key Files**:
- `agents/orchestrator_agent.py` (main coordination)
- `TEAM_SHREYA/load_to_neo4j.py` (graph ingestion)
- `TEAM_SHREYA/agent_orchestration.py` (to implement)

**Deliverables**:
- [ ] Complete orchestrator_agent.py
- [ ] Implement control agent logic
- [ ] Add experiment logging
- [ ] Setup Neo4j loader

---

### **🟢 Vishwa – Reinforcement Learning for Molecular Optimization**

**Responsibility**:
- Implement RL agent for molecule optimization
- Design reward function (efficacy + toxicity penalty + synthetic feasibility)
- Create molecular simulation environment
- Optimize molecules based on multiple objectives

**Key Files**:
- `agents/optimizer_agent.py` (base optimizer)
- `TEAM_VISHWA/rl_optimizer.py` (to implement)
- `TEAM_VISHWA/reward_function.py` (to implement)
- `TEAM_VISHWA/molecular_env.py` (to implement)

**Deliverables**:
- [ ] Complete RL agent implementation
- [ ] Design reward function
- [ ] Create simulation environment
- [ ] Add experiment comparisons (random vs RL)

---

### **🟡 Pranay – Graph Neural Network-based Molecular Generator**

**Responsibility**:
- Implement GNN for molecular generation
- Learn node and edge construction policies
- Implement chemical constraints (valence, ring closures)
- Calculate validity, novelty, diversity metrics

**Key Files**:
- `agents/generator_agent.py` (base generator)
- `TEAM_PRANAY/gnn_generator.py` (to implement)
- `TEAM_PRANAY/molecular_graph.py` (to implement)
- `TEAM_PRANAY/chemistry_constraints.py` (to implement)

**Deliverables**:
- [ ] Complete GNN generator
- [ ] Implement constraint checking
- [ ] Add validity/novelty/diversity metrics
- [ ] Integrate with orchestrator

---

## 📊 Integration Points

```
┌─────────────────────────────────────────────────┐
│         SHREYA: Orchestrator Agent              │
│  (Controls flow, makes decisions, logs events)  │
└────────┬──────────────────────────┬────────┬────┘
         │                          │        │
    ┌────▼─────┐          ┌────────▼──┐  ┌─▼─────────┐
    │  PRANAY  │          │  VISHWA   │  │ Predictor │
    │Generator │          │Optimizer  │  │  Agent    │
    │(Creates) │          │(Improves) │  │(Evaluates)│
    └──────────┘          └───────────┘  └───────────┘

Flow:
1. Generator creates molecules
2. Predictor evaluates properties
3. Optimizer refines molecules
4. Orchestrator decides next action
5. Orchestrator logs to Neo4j
```

---

## 🔗 How They Work Together

### **Pipeline Execution**

```python
# Shreya orchestrates this flow:

for generation in range(num_generations):
    # Step 1: Generate molecules (Pranay's code)
    molecules = pranay_generator.generate(num_mols=100)
    
    # Step 2: Evaluate (Predictor)
    scores = predictor.evaluate(molecules)
    
    # Step 3: Optimize (Vishwa's code)
    optimized = vishwa_optimizer.optimize(
        molecules=best_molecules,
        reward_fn=vishwa_reward,
        iterations=30
    )
    
    # Step 4: Log & decide next action (Shreya's logic)
    orchestrator.log_results(optimized, scores)
    should_continue = orchestrator.make_decision(scores)
```

---

## 📋 Status Tracking

### **What's Complete**
- [x] Base agent structure (all 5 agents)
- [x] Predictor agent (property prediction)
- [x] Ranker agent (ranking molecules)
- [x] Generator agent (skeleton with mutations)
- [x] Optimizer agent (skeleton with RL basics)

### **What's In Progress**
- [ ] **Shreya**: Orchestration control logic
- [ ] **Vishwa**: RL implementation & reward function
- [ ] **Pranay**: GNN implementation & constraints

---

## 📁 Folder Structure for Team

```
ChemAI/
├── agents/                          ← Core agent framework
│   ├── orchestrator_agent.py       ← Shreya's main file
│   ├── generator_agent.py          ← Pranay's base
│   ├── optimizer_agent.py          ← Vishwa's base
│   ├── predictor_agent.py
│   └── ranker_agent.py
│
├── TEAM_SHREYA/                    ← Shreya's implementations
│   ├── agent_orchestration.py      ← Control logic
│   ├── load_to_neo4j.py           ← Graph management
│   └── experiment_logger.py        ← Logging system
│
├── TEAM_VISHWA/                    ← Vishwa's implementations
│   ├── rl_optimizer.py            ← RL agent
│   ├── reward_function.py         ← Reward design
│   ├── molecular_env.py           ← Simulation env
│   └── experiments.py             ← RL experiments
│
└── TEAM_PRANAY/                    ← Pranay's implementations
    ├── gnn_generator.py           ← GNN model
    ├── molecular_graph.py         ← Graph operations
    ├── chemistry_constraints.py   ← Validation rules
    └── metrics.py                 ← Validity/novelty/diversity
```

---

## 📞 How to Collaborate

### **Read These First**
1. `README_TEAM_ORGANIZATION.txt` - Team overview
2. `TEAM_STATUS_AND_REMAINING_WORK.txt` - What's left

### **For Shreya**
- Focus on `orchestrator_agent.py`
- Review Neo4j integration
- See `SHREYA_KNOWLEDGE_GRAPH_COMPLETION.md`

### **For Vishwa**
- Implement RL agent in `TEAM_VISHWA/rl_optimizer.py`
- Design reward function
- Run experiments comparing search strategies

### **For Pranay**
- Implement GNN in `TEAM_PRANAY/gnn_generator.py`
- Add chemistry constraints
- Calculate metrics

---

## 🎯 Weekly Sync Points

- **Check**: `TEAM_STATUS_AND_REMAINING_WORK.txt`
- **Update**: Push your code to your team folder
- **Test**: Run `test_all_models.py`
- **Review**: Look at other team's integration points

---

## 🚀 Next Actions

**For Shreya**:
- Start with `SHREYA_KNOWLEDGE_GRAPH_COMPLETION.md`
- Implement decision logic in orchestrator

**For Vishwa**:
- Design RL reward function
- Create molecular simulation environment

**For Pranay**:
- Implement GNN architecture
- Add chemistry constraint validation

---

**Master Index**: Go back to `../INDEX.md`
