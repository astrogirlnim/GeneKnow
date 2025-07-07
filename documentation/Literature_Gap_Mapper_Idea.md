
# 📚 Literature Gap Mapper — Open Source AI Agent

## 🚀 Project Summary

**Literature Gap Mapper** is an AI-powered agent that ingests academic papers in a specific domain and generates a **visual map of existing knowledge** and, more importantly, **gaps in understanding**. This allows researchers, students, or biotech/pharma professionals to rapidly identify underexplored areas, conflicting results, or novel research opportunities.

It's like having an AI research assistant that builds you a conceptual overview of the field — surfacing *what's missing* just as easily as *what's known*.

---

## 🎯 Alignment with the Assignment

This project is directly aligned with the **G2P4 Open Source Week assignment** in the following ways:

### ✅ Context Management
- Maintains a long-term memory of concepts, claims, and paper findings.
- Connects ideas across dozens of documents using embeddings and knowledge graphs.

### ✅ Generation + Verification
- AI proposes claims, contradictions, and underexplored areas.
- Human researchers verify whether they agree or want to refine results.
- Slider to control confidence threshold or degree of suggestion autonomy.

### ✅ Incremental Processing
- Papers are chunked into abstract, method, result, and discussion blocks.
- Each paper is processed individually, then linked across documents.

### ✅ Visual Interface
- Outputs a **dynamic semantic graph** of key claims, topics, and contradictions.
- Users click on nodes to see which paper supports each idea.

### ✅ Partial Autonomy
- Users select how much context the agent should infer or link automatically.
- Slider to control the scope (narrow vs exploratory linking).

---

## 🧠 Why It Matters (Real World Impact)

- **Researchers** spend hours identifying gaps and synthesizing literature — this agent accelerates that.
- **Thesis writers** can instantly understand what's already been covered and where they can contribute.
- **Pharma/biotech** can identify promising angles for drug discovery or repurposing.
- **AI researchers** can use it to discover contradictory benchmarks or evaluations.

---

## 💡 Open Source Fit

- Works entirely on public research papers (e.g. PubMed, arXiv).
- Easy to fork and customize for different fields.
- Enables contributions from devs, scientists, students, and educators.
- Encourages extension: plug in new paper sources, domain-specific rules, or export formats.

---

## 🧩 Core Components (Team Work Division)

| Component | Description | Owner |
|----------|-------------|-------|
| Paper Ingestion | Parse PDFs or URLs (arXiv, PubMed) | Backend |
| LLM/Claim Extraction | Extract key claims, findings, and methods | AI/NLP |
| Concept Graph Builder | Generate knowledge map, identify contradictions | Backend/NLP |
| Visual Interface | Interactive semantic graph + verification UI | Frontend |
| Autonomy Controls | Slider-based tuning of AI inferences | Frontend/UX |

---

## 🧪 Example Use Case

**Goal**: Explore research on "Microbiome and Mental Health"

**Steps**:
1. Upload 20 papers from PubMed.
2. AI extracts key claims from each paper (e.g. “Increased diversity correlates with lower anxiety”).
3. System links overlapping or contradictory findings.
4. Researcher sees a map of clusters (gut bacteria → inflammation → mood).
5. One node highlights *no studies* found on a particular strain in adolescents → a gap.

---

## 🔧 Tech Stack Suggestion

- Python + LangChain / LlamaIndex
- OpenAI / Claude / local LLM support
- Neo4j or NetworkX for graph representation
- Next.js + D3.js for frontend
- GitHub + MIT License for open source

---

## 📍 Next Step: PRD Brainstorm

To build the PRD (Product Requirements Document), we’ll need to define:
- Primary user persona(s)
- Core user journeys
- MVP features vs stretch goals
- Evaluation criteria
- Risks and mitigations

See next steps in chat...
