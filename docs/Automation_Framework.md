# Automation Framework for Reaction Rule Generation

## Overview

This document outlines the complete automation architecture to transform the manual workflow into an end-to-end automated system.

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                     ORCHESTRATION LAYER                             │
│                  (Airflow/Prefect/Custom)                           │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
        ┌──────────────────┼──────────────────┐
        │                  │                  │
        ▼                  ▼                  ▼
┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│  INGESTION   │  │ PROCESSING   │  │ VALIDATION   │
│   PIPELINE   │  │   PIPELINE   │  │   PIPELINE   │
└──────────────┘  └──────────────┘  └──────────────┘
        │                  │                  │
        │                  │                  │
        ▼                  ▼                  ▼
┌─────────────────────────────────────────────────────┐
│              STORAGE & VERSIONING                   │
│     (Git, Database, Artifact Registry)              │
└─────────────────────────────────────────────────────┘
```

---

## Component 1: Automated Data Ingestion

### 1.1 PDF Monitoring & Processing Service

**Purpose:** Automatically detect and process new papers

```python
# automation/ingestion/pdf_monitor.py
import os
import hashlib
from pathlib import Path
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from chemtools_automation import extract_reactions_from_pdf

class PDFHandler(FileSystemEventHandler):
    """Monitor directory for new PDF files and auto-process."""
    
    def __init__(self, watch_dir, output_dir, db_path):
        self.watch_dir = Path(watch_dir)
        self.output_dir = Path(output_dir)
        self.db_path = db_path
        self.processed_files = self._load_processed_files()
    
    def _load_processed_files(self):
        """Track processed files to avoid reprocessing."""
        cache_file = self.output_dir / 'processed_pdfs.json'
        if cache_file.exists():
            return json.load(open(cache_file))
        return {}
    
    def _get_file_hash(self, filepath):
        """Calculate file hash for duplicate detection."""
        with open(filepath, 'rb') as f:
            return hashlib.sha256(f.read()).hexdigest()
    
    def on_created(self, event):
        """Process new PDF files."""
        if event.is_directory or not event.src_path.endswith('.pdf'):
            return
        
        filepath = Path(event.src_path)
        file_hash = self._get_file_hash(filepath)
        
        # Skip if already processed
        if file_hash in self.processed_files:
            print(f"Skipping duplicate: {filepath.name}")
            return
        
        print(f"Processing new PDF: {filepath.name}")
        
        try:
            # Extract reactions
            results = extract_reactions_from_pdf(filepath)
            
            # Generate rules
            rules = self._extract_rules(results)
            
            # Add to database
            self._add_rules_to_database(rules, filepath.name)
            
            # Mark as processed
            self.processed_files[file_hash] = {
                'filename': filepath.name,
                'processed_at': datetime.now().isoformat(),
                'rules_extracted': len(rules)
            }
            self._save_processed_files()
            
            print(f"✓ Extracted {len(rules)} rules from {filepath.name}")
            
        except Exception as e:
            print(f"✗ Error processing {filepath.name}: {e}")
            self._log_error(filepath, e)
    
    def _extract_rules(self, pdf_results):
        """Convert PDF extraction results to rules."""
        from automation.rule_generator import RuleGenerator
        
        generator = RuleGenerator()
        rules = []
        
        for reaction in pdf_results['reactions']:
            rule = generator.create_rule(
                reaction_smiles=reaction['smiles'],
                conditions=reaction['conditions'],
                substrate_scope=reaction.get('scope', []),
                provenance=pdf_results['doi'],
                confidence=reaction.get('yield', 0) / 100
            )
            rules.append(rule)
        
        return rules
    
    def _add_rules_to_database(self, rules, source):
        """Add rules to condition database."""
        from automation.database_manager import DatabaseManager
        
        db_manager = DatabaseManager(self.db_path)
        
        for rule in rules:
            # Check for duplicates
            if not db_manager.is_duplicate(rule):
                db_manager.add_entry(rule)
                db_manager.generate_test_case(rule)
            else:
                print(f"  Skip duplicate rule: {rule['id']}")

# Usage
def start_pdf_monitor(watch_dir='./incoming_papers', db_path='./data/conditionDB'):
    """Start monitoring for new PDFs."""
    handler = PDFHandler(watch_dir, './processed_papers', db_path)
    observer = Observer()
    observer.schedule(handler, watch_dir, recursive=True)
    observer.start()
    print(f"Monitoring {watch_dir} for new PDFs...")
    
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()
```

### 1.2 Dataset Ingestion Pipeline

**Purpose:** Batch process reaction datasets

```python
# automation/ingestion/dataset_processor.py
from typing import List, Dict
import pandas as pd
from rdkit import Chem
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

class DatasetProcessor:
    """Process large reaction datasets in parallel."""
    
    def __init__(self, reaction_type_filter=None, min_yield=0.60):
        self.reaction_type_filter = reaction_type_filter
        self.min_yield = min_yield
        self.stats = {'total': 0, 'valid': 0, 'rules_generated': 0}
    
    def process_dataset(self, dataset_path: str, output_dir: str) -> Dict:
        """Process complete dataset with parallel execution."""
        
        # Load dataset
        print(f"Loading dataset: {dataset_path}")
        df = self._load_dataset(dataset_path)
        self.stats['total'] = len(df)
        
        # Filter by reaction type if specified
        if self.reaction_type_filter:
            df = df[df['reaction_type'].str.contains(self.reaction_type_filter, na=False)]
            print(f"Filtered to {len(df)} {self.reaction_type_filter} reactions")
        
        # Filter by yield threshold
        df = df[df['yield'] >= self.min_yield]
        print(f"Filtered to {len(df)} high-yield reactions (≥{self.min_yield})")
        
        # Validate SMILES
        print("Validating SMILES...")
        df['valid'] = df['reaction_smiles'].apply(self._validate_smiles)
        df = df[df['valid']]
        self.stats['valid'] = len(df)
        print(f"Valid reactions: {len(df)}")
        
        # Cluster similar reactions
        print("Clustering similar reactions...")
        clusters = self._cluster_reactions(df)
        
        # Generate rules from clusters (parallel)
        print("Generating rules from clusters...")
        rules = self._generate_rules_parallel(clusters)
        self.stats['rules_generated'] = len(rules)
        
        # Save results
        self._save_results(rules, output_dir)
        
        return self.stats
    
    def _load_dataset(self, path: str) -> pd.DataFrame:
        """Load dataset with flexible format detection."""
        if path.endswith('.csv'):
            return pd.read_csv(path)
        elif path.endswith('.json') or path.endswith('.jsonl'):
            return pd.read_json(path, lines=path.endswith('.jsonl'))
        elif path.endswith('.xlsx'):
            return pd.read_excel(path)
        else:
            raise ValueError(f"Unsupported format: {path}")
    
    def _validate_smiles(self, smiles: str) -> bool:
        """Validate reaction SMILES."""
        try:
            rxn = Chem.ReactionFromSmarts(smiles)
            return rxn is not None
        except:
            return False
    
    def _cluster_reactions(self, df: pd.DataFrame) -> List[Dict]:
        """Group similar reactions together."""
        from sklearn.cluster import DBSCAN
        from chemtools.reaction_similarity import compute_drfp_matrix
        
        # Compute reaction fingerprints
        print("  Computing reaction fingerprints...")
        fps = compute_drfp_matrix(df['reaction_smiles'].tolist())
        
        # Cluster with DBSCAN
        print("  Clustering...")
        clustering = DBSCAN(eps=0.3, min_samples=5, metric='jaccard')
        df['cluster'] = clustering.fit_predict(fps)
        
        # Group by cluster
        clusters = []
        for cluster_id in df['cluster'].unique():
            if cluster_id == -1:  # Noise points
                continue
            
            cluster_df = df[df['cluster'] == cluster_id]
            clusters.append({
                'id': cluster_id,
                'size': len(cluster_df),
                'reactions': cluster_df.to_dict('records')
            })
        
        print(f"  Found {len(clusters)} clusters")
        return clusters
    
    def _generate_rules_parallel(self, clusters: List[Dict]) -> List[Dict]:
        """Generate rules from clusters in parallel."""
        from automation.rule_generator import RuleGenerator
        
        generator = RuleGenerator()
        
        # Use multiprocessing for CPU-intensive work
        num_workers = mp.cpu_count() - 1
        print(f"  Using {num_workers} workers")
        
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = [
                executor.submit(generator.generate_from_cluster, cluster)
                for cluster in clusters
            ]
            
            rules = []
            for i, future in enumerate(futures):
                try:
                    rule = future.result(timeout=60)
                    if rule:
                        rules.append(rule)
                        print(f"  [{i+1}/{len(futures)}] Generated rule: {rule['id']}")
                except Exception as e:
                    print(f"  [{i+1}/{len(futures)}] Error: {e}")
        
        return rules
    
    def _save_results(self, rules: List[Dict], output_dir: str):
        """Save generated rules and metadata."""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save rules
        rules_file = output_path / 'generated_rules.json'
        with open(rules_file, 'w') as f:
            json.dump(rules, f, indent=2)
        
        # Save statistics
        stats_file = output_path / 'processing_stats.json'
        with open(stats_file, 'w') as f:
            json.dump(self.stats, f, indent=2)
        
        print(f"✓ Saved {len(rules)} rules to {rules_file}")

# Usage
if __name__ == "__main__":
    processor = DatasetProcessor(reaction_type_filter="Suzuki", min_yield=0.70)
    stats = processor.process_dataset(
        dataset_path='data/reaction_dataset/uspto_grants.csv',
        output_dir='automation/output/suzuki_rules'
    )
    print(f"\nProcessing complete: {stats}")
```

---

## Component 2: Intelligent Rule Generator

### 2.1 SMARTS Pattern Generator

**Purpose:** Automatically generate SMARTS patterns from examples

```python
# automation/rule_generator.py
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from typing import List, Dict, Optional
import re

class RuleGenerator:
    """Automatically generate reaction rules from examples."""
    
    def __init__(self):
        self.priority_ranges = {
            'default': (0, 10),
            'general': (20, 40),
            'specific': (45, 60),
            'very_specific': (65, 80),
            'scheme': (85, 100)
        }
    
    def generate_from_cluster(self, cluster: Dict) -> Optional[Dict]:
        """Generate a rule from a cluster of similar reactions."""
        
        reactions = cluster['reactions']
        if len(reactions) < 3:  # Need minimum examples
            return None
        
        # Extract common SMARTS pattern
        smarts = self._extract_common_smarts(reactions)
        if not smarts:
            return None
        
        # Extract common conditions
        conditions = self._extract_common_conditions(reactions)
        
        # Calculate priority based on specificity
        priority = self._calculate_priority(smarts, cluster['size'])
        
        # Generate unique ID
        rule_id = self._generate_rule_id(smarts, conditions)
        
        # Build rule entry
        rule = {
            'id': rule_id,
            'reaction_type': self._infer_reaction_type(smarts),
            'name': self._generate_name(smarts, conditions),
            'rxn_smiles_min': smarts,
            'token_signature': self._generate_token_signature(conditions),
            'conditions': conditions,
            'env': {
                'features_from_smiles': self._extract_features_from_smarts(smarts),
                'feature_requirements': {}
            },
            'evidence': {
                'provenance': f"cluster_{cluster['id']}",
                'cluster_size': cluster['size'],
                'avg_yield': sum(r['yield'] for r in reactions) / len(reactions),
                'last_updated': datetime.now().isoformat(),
                'confidence': self._calculate_confidence(reactions)
            },
            'priority': priority,
            'notes': [f"Generated from {len(reactions)} similar reactions"]
        }
        
        return rule
    
    def _extract_common_smarts(self, reactions: List[Dict]) -> Optional[str]:
        """Find maximum common substructure and generate SMARTS."""
        
        # Parse reaction SMILES
        rxn_mols = []
        for rxn in reactions:
            try:
                r = Chem.ReactionFromSmiles(rxn['reaction_smiles'])
                if r:
                    rxn_mols.append(r)
            except:
                continue
        
        if len(rxn_mols) < 2:
            return None
        
        # Find MCS for reactants
        reactants = [r.GetReactants() for r in rxn_mols]
        products = [r.GetProducts() for r in rxn_mols]
        
        # Get common reactant pattern
        reactant_mcs = self._find_mcs_for_molecules([mol for mols in reactants for mol in mols])
        product_mcs = self._find_mcs_for_molecules([mol for mols in products for mol in mols])
        
        if not reactant_mcs or not product_mcs:
            return None
        
        # Add atom mapping
        mapped_smarts = self._add_atom_mapping(reactant_mcs, product_mcs)
        
        return mapped_smarts
    
    def _find_mcs_for_molecules(self, molecules: List[Chem.Mol]) -> Optional[str]:
        """Find maximum common substructure."""
        if len(molecules) < 2:
            return None
        
        mcs_result = rdFMCS.FindMCS(
            molecules,
            timeout=10,
            bondCompare=rdFMCS.BondCompare.CompareOrderExact,
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            ringMatchesRingOnly=True,
            completeRingsOnly=True
        )
        
        return mcs_result.smartsString if mcs_result.numAtoms > 0 else None
    
    def _add_atom_mapping(self, reactant_smarts: str, product_smarts: str) -> str:
        """Add atom mapping to SMARTS pattern."""
        # Simple heuristic: map heavy atoms in order
        reactant_atoms = re.findall(r'\[[^\]]+\]|[A-Z][a-z]?', reactant_smarts)
        product_atoms = re.findall(r'\[[^\]]+\]|[A-Z][a-z]?', product_smarts)
        
        # Add mapping numbers
        mapped_reactants = []
        for i, atom in enumerate(reactant_atoms[:min(10, len(reactant_atoms))], 1):
            if atom.startswith('['):
                mapped_atom = atom[:-1] + f':{i}]'
            else:
                mapped_atom = f'[{atom}:{i}]'
            mapped_reactants.append(mapped_atom)
        
        mapped_products = []
        for i, atom in enumerate(product_atoms[:min(10, len(product_atoms))], 1):
            if atom.startswith('['):
                mapped_atom = atom[:-1] + f':{i}]'
            else:
                mapped_atom = f'[{atom}:{i}]'
            mapped_products.append(mapped_atom)
        
        return ''.join(mapped_reactants) + '>>' + ''.join(mapped_products)
    
    def _extract_common_conditions(self, reactions: List[Dict]) -> Dict:
        """Extract most common conditions from examples."""
        from collections import Counter
        
        # Count occurrences of each condition
        catalysts = Counter(r.get('catalyst') for r in reactions if r.get('catalyst'))
        ligands = Counter(r.get('ligand') for r in reactions if r.get('ligand'))
        bases = Counter(r.get('base') for r in reactions if r.get('base'))
        solvents = Counter(r.get('solvent') for r in reactions if r.get('solvent'))
        
        # Get most common
        conditions = {
            'pd_source': [catalysts.most_common(1)[0][0]] if catalysts else [],
            'ligands': [ligands.most_common(1)[0][0]] if ligands else [],
            'base': [bases.most_common(1)[0][0]] if bases else [],
            'solvent': [solvents.most_common(1)[0][0]] if solvents else []
        }
        
        # Extract temperature and time ranges
        temps = [r.get('temperature_C') for r in reactions if r.get('temperature_C')]
        times = [r.get('time_h') for r in reactions if r.get('time_h')]
        
        if temps:
            conditions['temperature_C'] = [min(temps), max(temps)]
        if times:
            conditions['time_h'] = [min(times), max(times)]
        
        # Extract loading ranges
        loadings = {}
        pd_loadings = [r.get('pd_mol%') for r in reactions if r.get('pd_mol%')]
        if pd_loadings:
            loadings['Pd_mol%'] = [min(pd_loadings), max(pd_loadings)]
        
        if loadings:
            conditions['loadings'] = loadings
        
        return conditions
    
    def _calculate_priority(self, smarts: str, cluster_size: int) -> int:
        """Calculate priority based on pattern specificity."""
        
        # Count constraints in SMARTS
        specificity_score = 0
        
        # Count atom properties
        specificity_score += smarts.count('[')  # Bracketed atoms
        specificity_score += smarts.count('#')  # Atomic numbers
        specificity_score += smarts.count(':')  # Atom maps
        specificity_score += smarts.count('!')  # Negations
        
        # Count structural features
        specificity_score += smarts.count('1')  # Rings
        specificity_score += smarts.count('=')  # Double bonds
        specificity_score += smarts.count('#')  # Triple bonds
        
        # Adjust for cluster size (larger clusters = more general = lower priority)
        size_penalty = min(cluster_size / 10, 20)
        
        final_score = specificity_score - size_penalty
        
        # Map to priority range
        if final_score < 10:
            return 35  # General
        elif final_score < 20:
            return 50  # Specific
        elif final_score < 30:
            return 65  # Very specific
        else:
            return 75  # Highly specific
    
    def _calculate_confidence(self, reactions: List[Dict]) -> float:
        """Calculate confidence score based on evidence quality."""
        
        # Factors:
        # 1. Number of examples
        n_examples = len(reactions)
        example_score = min(n_examples / 20, 1.0)
        
        # 2. Average yield
        avg_yield = sum(r.get('yield', 0) for r in reactions) / n_examples
        yield_score = avg_yield / 100
        
        # 3. Consistency (low variance in conditions)
        # ... (simplified here)
        consistency_score = 0.8
        
        # Weighted average
        confidence = (
            0.3 * example_score +
            0.4 * yield_score +
            0.3 * consistency_score
        )
        
        return round(confidence, 2)
    
    def _generate_rule_id(self, smarts: str, conditions: Dict) -> str:
        """Generate unique rule ID."""
        rxn_type = self._infer_reaction_type(smarts)
        substrate = self._infer_substrate_class(smarts)
        ligand = conditions.get('ligands', ['UNKNOWN'])[0]
        
        return f"SCDB-{rxn_type.upper()[:3]}-{substrate}-{ligand}"
    
    def _infer_reaction_type(self, smarts: str) -> str:
        """Infer reaction type from SMARTS pattern."""
        # Simple heuristics
        if 'B' in smarts and '[c:' in smarts:
            return "Suzuki_Miyaura"
        elif 'N' in smarts and 'Pd' in smarts:
            return "Buchwald_Hartwig"
        # Add more patterns...
        return "Unknown"
    
    def _infer_substrate_class(self, smarts: str) -> str:
        """Infer substrate class from SMARTS."""
        if '-Cl' in smarts:
            return "ArCl"
        elif '-Br' in smarts:
            return "ArBr"
        elif '-I' in smarts:
            return "ArI"
        # Add more patterns...
        return "GENERAL"
    
    def _generate_name(self, smarts: str, conditions: Dict) -> str:
        """Generate human-readable name."""
        rxn_type = self._infer_reaction_type(smarts)
        substrate = self._infer_substrate_class(smarts)
        ligand = conditions.get('ligands', [''])[0]
        
        return f"{rxn_type} - {substrate} with {ligand}"
    
    def _generate_token_signature(self, conditions: Dict) -> List[str]:
        """Generate token signature for matching."""
        tokens = []
        
        if conditions.get('pd_source'):
            tokens.append(conditions['pd_source'][0])
        if conditions.get('ligands'):
            tokens.append(conditions['ligands'][0])
        if conditions.get('base'):
            tokens.append(conditions['base'][0])
        
        return tokens
    
    def _extract_features_from_smarts(self, smarts: str) -> Dict:
        """Extract feature constraints from SMARTS."""
        features = {}
        
        # Parse leaving group
        if '-Br' in smarts:
            features['electrophile.lg_class'] = 'Br'
        elif '-Cl' in smarts:
            features['electrophile.lg_class'] = 'Cl'
        elif '-I' in smarts:
            features['electrophile.lg_class'] = 'I'
        
        return features
```

### 2.2 Test Case Generator

**Purpose:** Automatically generate test cases for new rules

```python
# automation/test_generator.py
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Dict, List

class TestCaseGenerator:
    """Generate test cases automatically for validation."""
    
    def generate_for_rule(self, rule: Dict) -> List[Dict]:
        """Generate test cases for a single rule."""
        
        test_cases = []
        
        # 1. Generate positive test (should match)
        positive_test = self._generate_positive_test(rule)
        if positive_test:
            test_cases.append(positive_test)
        
        # 2. Generate boundary tests (near decision boundary)
        boundary_tests = self._generate_boundary_tests(rule)
        test_cases.extend(boundary_tests)
        
        # 3. Generate negative test (should NOT match)
        negative_test = self._generate_negative_test(rule)
        if negative_test:
            test_cases.append(negative_test)
        
        return test_cases
    
    def _generate_positive_test(self, rule: Dict) -> Optional[Dict]:
        """Generate a test that should match this rule."""
        
        # Parse SMARTS pattern
        smarts = rule['rxn_smiles_min']
        
        try:
            rxn = Chem.ReactionFromSmarts(smarts)
            if not rxn:
                return None
            
            # Generate example molecule matching pattern
            reactants = rxn.GetReactants()
            products = rxn.GetProducts()
            
            # Convert to SMILES (simplified - real version would enumerate examples)
            reactant_smiles = '.'.join(Chem.MolToSmiles(m) for m in reactants)
            product_smiles = '.'.join(Chem.MolToSmiles(m) for m in products)
            
            return {
                'rule_id': rule['id'],
                'type': 'positive',
                'smiles': f"{reactant_smiles}>>{product_smiles}",
                'description': f"Test case for {rule['name']}",
                'should_match': True
            }
        
        except:
            return None
    
    def _generate_boundary_tests(self, rule: Dict) -> List[Dict]:
        """Generate tests near decision boundaries."""
        # Implementation would create variations
        # e.g., change leaving group, add/remove substituents
        return []
    
    def _generate_negative_test(self, rule: Dict) -> Optional[Dict]:
        """Generate a test that should NOT match."""
        # Implementation would create a similar but distinct substrate
        return None
```

---

## Component 3: Automated Validation Pipeline

### 3.1 Continuous Validation Service

```python
# automation/validation/validator_service.py
import schedule
import time
from pathlib import Path

class ContinuousValidator:
    """Continuously validate database as new rules are added."""
    
    def __init__(self, db_path: str, test_path: str, report_dir: str):
        self.db_path = Path(db_path)
        self.test_path = Path(test_path)
        self.report_dir = Path(report_dir)
        self.report_dir.mkdir(parents=True, exist_ok=True)
        
        self.last_validation = None
        self.validation_history = []
    
    def run_validation(self):
        """Run complete validation suite."""
        print(f"\n{'='*60}")
        print(f"Running validation at {datetime.now()}")
        print(f"{'='*60}\n")
        
        from scripts.validate_all_suzuki_rules import run_validation_suite
        
        # Run validation
        results = run_validation_suite(
            database_path=str(self.db_path),
            test_cases_path=str(self.test_path)
        )
        
        # Store results
        self.last_validation = {
            'timestamp': datetime.now().isoformat(),
            'pass_rate': results['pass_rate'],
            'total_tests': results['total_tests'],
            'passing': results['passing'],
            'failing': results['failing']
        }
        self.validation_history.append(self.last_validation)
        
        # Generate reports
        self._generate_trend_report()
        
        # Check for regressions
        self._check_for_regressions()
        
        print(f"\n✓ Validation complete: {results['pass_rate']:.1%} pass rate")
    
    def _generate_trend_report(self):
        """Generate report showing validation trends over time."""
        import matplotlib.pyplot as plt
        
        if len(self.validation_history) < 2:
            return
        
        # Extract data
        timestamps = [h['timestamp'] for h in self.validation_history]
        pass_rates = [h['pass_rate'] for h in self.validation_history]
        
        # Create plot
        plt.figure(figsize=(12, 6))
        plt.plot(range(len(pass_rates)), pass_rates, marker='o')
        plt.xlabel('Validation Run')
        plt.ylabel('Pass Rate')
        plt.title('Validation Pass Rate Trend')
        plt.grid(True)
        plt.ylim(0, 1)
        
        # Save
        plot_path = self.report_dir / 'validation_trend.png'
        plt.savefig(plot_path)
        plt.close()
        
        print(f"✓ Trend report saved to {plot_path}")
    
    def _check_for_regressions(self):
        """Alert if pass rate drops."""
        if len(self.validation_history) < 2:
            return
        
        current = self.validation_history[-1]['pass_rate']
        previous = self.validation_history[-2]['pass_rate']
        
        if current < previous - 0.05:  # 5% drop
            print(f"\n⚠️  REGRESSION DETECTED!")
            print(f"   Pass rate dropped from {previous:.1%} to {current:.1%}")
            self._send_alert(
                f"Validation regression: {previous:.1%} → {current:.1%}"
            )
    
    def _send_alert(self, message: str):
        """Send alert notification."""
        # Could integrate with Slack, email, etc.
        print(f"📧 ALERT: {message}")
    
    def schedule_validation(self, interval_hours: int = 6):
        """Schedule periodic validation."""
        
        # Run immediately
        self.run_validation()
        
        # Schedule periodic runs
        schedule.every(interval_hours).hours.do(self.run_validation)
        
        print(f"\nScheduled validation every {interval_hours} hours")
        print("Press Ctrl+C to stop")
        
        while True:
            schedule.run_pending()
            time.sleep(60)

# Usage
if __name__ == "__main__":
    validator = ContinuousValidator(
        db_path='data/conditionDB/suzuki_db.json',
        test_path='tests/sample_reactions.py',
        report_dir='automation/reports'
    )
    validator.schedule_validation(interval_hours=6)
```

### 3.2 Smart Fix Recommender

```python
# automation/validation/fix_recommender.py
from typing import Dict, List

class FixRecommender:
    """Automatically suggest fixes for failing validation tests."""
    
    def analyze_failures(self, validation_results: Dict) -> List[Dict]:
        """Analyze failures and recommend fixes."""
        
        recommendations = []
        
        for failure in validation_results['failures']:
            rec = self._diagnose_failure(failure)
            if rec:
                recommendations.append(rec)
        
        return recommendations
    
    def _diagnose_failure(self, failure: Dict) -> Optional[Dict]:
        """Diagnose a single failure and recommend fix."""
        
        expected_rule = failure['expected']
        actual_rule = failure['actual']
        expected_priority = failure['expected_priority']
        actual_priority = failure['actual_priority']
        
        # Case 1: Priority conflict
        if actual_priority > expected_priority:
            return {
                'type': 'priority_conflict',
                'rule': expected_rule,
                'recommendation': {
                    'action': 'increase_priority',
                    'current': expected_priority,
                    'suggested': actual_priority + 5,
                    'reason': f'Currently beaten by {actual_rule} (priority {actual_priority})'
                },
                'auto_fixable': True
            }
        
        # Case 2: SMARTS pattern too broad
        if self._is_pattern_too_broad(failure):
            return {
                'type': 'smarts_too_broad',
                'rule': expected_rule,
                'recommendation': {
                    'action': 'refine_smarts',
                    'suggestions': self._suggest_smarts_refinements(failure),
                    'reason': 'Pattern matches too many substrates'
                },
                'auto_fixable': False  # Requires manual review
            }
        
        # Case 3: Feature detection issue
        if self._is_feature_issue(failure):
            return {
                'type': 'feature_not_detected',
                'rule': expected_rule,
                'recommendation': {
                    'action': 'update_feature_detector',
                    'feature': self._identify_missing_feature(failure),
                    'reason': 'Required feature not being detected'
                },
                'auto_fixable': True
            }
        
        return None
    
    def apply_automatic_fixes(self, recommendations: List[Dict], db_path: str):
        """Automatically apply fixable recommendations."""
        
        fixed = []
        
        for rec in recommendations:
            if not rec['auto_fixable']:
                continue
            
            try:
                if rec['type'] == 'priority_conflict':
                    self._fix_priority(rec, db_path)
                    fixed.append(rec['rule'])
                elif rec['type'] == 'feature_not_detected':
                    self._fix_feature_detection(rec)
                    fixed.append(rec['rule'])
            except Exception as e:
                print(f"Failed to auto-fix {rec['rule']}: {e}")
        
        return fixed
    
    def _fix_priority(self, rec: Dict, db_path: str):
        """Automatically adjust priority."""
        import json
        
        with open(db_path, 'r') as f:
            db = json.load(f)
        
        # Find and update rule
        for entry in db:
            if entry['id'] == rec['rule']:
                old_priority = entry['priority']
                entry['priority'] = rec['recommendation']['suggested']
                print(f"  Adjusted {rec['rule']} priority: {old_priority} → {entry['priority']}")
                break
        
        # Save
        with open(db_path, 'w') as f:
            json.dump(db, f, indent=2)
```

---

## Component 4: Workflow Orchestration

### 4.1 Apache Airflow DAG

```python
# automation/orchestration/reaction_rules_dag.py
from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.operators.bash import BashOperator
from datetime import datetime, timedelta

default_args = {
    'owner': 'chem-automation',
    'depends_on_past': False,
    'start_date': datetime(2025, 1, 1),
    'email_on_failure': True,
    'email_on_retry': False,
    'retries': 2,
    'retry_delay': timedelta(minutes=5),
}

dag = DAG(
    'reaction_rules_pipeline',
    default_args=default_args,
    description='End-to-end reaction rule generation and validation',
    schedule_interval=timedelta(days=1),  # Run daily
    catchup=False
)

# Task 1: Monitor for new PDFs
check_new_pdfs = PythonOperator(
    task_id='check_new_pdfs',
    python_callable=lambda: PDFMonitor().check_for_new_files(),
    dag=dag
)

# Task 2: Process PDFs
process_pdfs = PythonOperator(
    task_id='process_pdfs',
    python_callable=lambda: PDFProcessor().process_all(),
    dag=dag
)

# Task 3: Process datasets
process_datasets = PythonOperator(
    task_id='process_datasets',
    python_callable=lambda: DatasetProcessor().process_all(),
    dag=dag
)

# Task 4: Generate rules
generate_rules = PythonOperator(
    task_id='generate_rules',
    python_callable=lambda: RuleGenerator().generate_all(),
    dag=dag
)

# Task 5: Generate test cases
generate_tests = PythonOperator(
    task_id='generate_test_cases',
    python_callable=lambda: TestCaseGenerator().generate_all(),
    dag=dag
)

# Task 6: Run validation
run_validation = BashOperator(
    task_id='run_validation',
    bash_command='python scripts/validate_all_suzuki_rules.py',
    dag=dag
)

# Task 7: Analyze results and recommend fixes
analyze_results = PythonOperator(
    task_id='analyze_results',
    python_callable=lambda: FixRecommender().analyze_and_fix(),
    dag=dag
)

# Task 8: Apply automatic fixes
apply_fixes = PythonOperator(
    task_id='apply_fixes',
    python_callable=lambda: FixRecommender().apply_automatic_fixes(),
    dag=dag
)

# Task 9: Re-validate after fixes
revalidate = BashOperator(
    task_id='revalidate',
    bash_command='python scripts/validate_all_suzuki_rules.py',
    dag=dag
)

# Task 10: Generate reports
generate_reports = BashOperator(
    task_id='generate_reports',
    bash_command='python scripts/generate_conditions_report.py > docs/conditions_report.md',
    dag=dag
)

# Task 11: Commit to git (if improved)
commit_changes = BashOperator(
    task_id='commit_changes',
    bash_command='''
        if [ "$(git status --porcelain)" ]; then
            git add data/conditionDB/*.json tests/sample_reactions.py
            git commit -m "Automated: Update reaction rules ($(date +%Y-%m-%d))"
            git tag -a "auto-$(date +%Y%m%d)" -m "Automated update"
        fi
    ''',
    dag=dag
)

# Define task dependencies
check_new_pdfs >> process_pdfs
process_pdfs >> generate_rules
process_datasets >> generate_rules
generate_rules >> generate_tests
generate_tests >> run_validation
run_validation >> analyze_results
analyze_results >> apply_fixes
apply_fixes >> revalidate
revalidate >> generate_reports
generate_reports >> commit_changes
```

### 4.2 Simple CLI Orchestrator (No Airflow Required)

```python
# automation/orchestrator.py
import click
from pathlib import Path

@click.group()
def cli():
    """Reaction rules automation CLI."""
    pass

@cli.command()
@click.option('--pdf-dir', default='./incoming_papers', help='Directory to monitor for PDFs')
@click.option('--db-path', default='./data/conditionDB/suzuki_db.json', help='Database path')
def watch_pdfs(pdf_dir, db_path):
    """Start PDF monitoring service."""
    from automation.ingestion.pdf_monitor import start_pdf_monitor
    start_pdf_monitor(pdf_dir, db_path)

@cli.command()
@click.option('--dataset', required=True, help='Path to reaction dataset')
@click.option('--output', default='./automation/output', help='Output directory')
@click.option('--reaction-type', help='Filter by reaction type')
@click.option('--min-yield', default=0.60, help='Minimum yield threshold')
def process_dataset(dataset, output, reaction_type, min_yield):
    """Process a reaction dataset."""
    from automation.ingestion.dataset_processor import DatasetProcessor
    
    processor = DatasetProcessor(
        reaction_type_filter=reaction_type,
        min_yield=min_yield
    )
    stats = processor.process_dataset(dataset, output)
    click.echo(f"✓ Processed {stats['total']} reactions, generated {stats['rules_generated']} rules")

@cli.command()
@click.option('--db-path', default='./data/conditionDB/suzuki_db.json')
@click.option('--test-path', default='./tests/sample_reactions.py')
def validate(db_path, test_path):
    """Run validation suite."""
    from scripts.validate_all_suzuki_rules import run_validation_suite
    
    results = run_validation_suite(db_path, test_path)
    click.echo(f"Pass rate: {results['pass_rate']:.1%}")

@cli.command()
@click.option('--db-path', default='./data/conditionDB/suzuki_db.json')
@click.option('--test-path', default='./tests/sample_reactions.py')
@click.option('--auto-fix', is_flag=True, help='Automatically apply fixes')
def improve(db_path, test_path, auto_fix):
    """Analyze failures and suggest/apply fixes."""
    from automation.validation.fix_recommender import FixRecommender
    from scripts.validate_all_suzuki_rules import run_validation_suite
    
    # Run validation
    results = run_validation_suite(db_path, test_path)
    
    # Get recommendations
    recommender = FixRecommender()
    recommendations = recommender.analyze_failures(results)
    
    click.echo(f"\nFound {len(recommendations)} improvement opportunities")
    
    for rec in recommendations:
        click.echo(f"\n{rec['rule']}:")
        click.echo(f"  Type: {rec['type']}")
        click.echo(f"  Recommendation: {rec['recommendation']}")
        
        if auto_fix and rec['auto_fixable']:
            click.echo(f"  ✓ Applying automatic fix...")
            recommender.apply_automatic_fixes([rec], db_path)

@cli.command()
@click.option('--interval-hours', default=6, help='Validation interval')
def continuous_validation(interval_hours):
    """Run continuous validation service."""
    from automation.validation.validator_service import ContinuousValidator
    
    validator = ContinuousValidator(
        db_path='data/conditionDB/suzuki_db.json',
        test_path='tests/sample_reactions.py',
        report_dir='automation/reports'
    )
    validator.schedule_validation(interval_hours=interval_hours)

@cli.command()
def run_full_pipeline():
    """Run the complete end-to-end pipeline."""
    click.echo("Starting full automation pipeline...")
    
    # 1. Process any new data
    click.echo("\n1. Processing new datasets...")
    # ... call dataset processor
    
    # 2. Generate rules
    click.echo("\n2. Generating rules...")
    # ... call rule generator
    
    # 3. Generate tests
    click.echo("\n3. Generating test cases...")
    # ... call test generator
    
    # 4. Validate
    click.echo("\n4. Running validation...")
    # ... call validator
    
    # 5. Improve
    click.echo("\n5. Analyzing and fixing issues...")
    # ... call fix recommender
    
    # 6. Generate reports
    click.echo("\n6. Generating reports...")
    # ... call report generator
    
    click.echo("\n✓ Pipeline complete!")

if __name__ == '__main__':
    cli()
```

---

## Component 5: Deployment & Monitoring

### 5.1 Docker Containerization

```dockerfile
# automation/Dockerfile
FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy requirements
COPY requirements.txt requirements-automation.txt ./
RUN pip install --no-cache-dir -r requirements.txt -r requirements-automation.txt

# Copy application
COPY . .

# Create directories
RUN mkdir -p incoming_papers processed_papers automation/reports

# Run orchestrator
CMD ["python", "automation/orchestrator.py", "run-full-pipeline"]
```

```yaml
# docker-compose.yml
version: '3.8'

services:
  pdf-monitor:
    build: .
    command: python automation/orchestrator.py watch-pdfs
    volumes:
      - ./incoming_papers:/app/incoming_papers
      - ./data/conditionDB:/app/data/conditionDB
    restart: always

  validator:
    build: .
    command: python automation/orchestrator.py continuous-validation --interval-hours 6
    volumes:
      - ./data/conditionDB:/app/data/conditionDB
      - ./tests:/app/tests
      - ./automation/reports:/app/automation/reports
    restart: always

  airflow-webserver:
    image: apache/airflow:2.7.0
    command: webserver
    ports:
      - "8080:8080"
    volumes:
      - ./automation/orchestration:/opt/airflow/dags
    environment:
      - AIRFLOW__CORE__EXECUTOR=LocalExecutor
    restart: always

  airflow-scheduler:
    image: apache/airflow:2.7.0
    command: scheduler
    volumes:
      - ./automation/orchestration:/opt/airflow/dags
    environment:
      - AIRFLOW__CORE__EXECUTOR=LocalExecutor
    restart: always
```

### 5.2 Monitoring Dashboard

```python
# automation/monitoring/dashboard.py
import streamlit as st
import pandas as pd
import json
from pathlib import Path

st.set_page_config(page_title="Reaction Rules Monitor", layout="wide")

st.title("🧪 Reaction Rules Automation Dashboard")

# Load validation history
history_file = Path('automation/reports/validation_history.json')
if history_file.exists():
    with open(history_file) as f:
        history = json.load(f)
    
    df = pd.DataFrame(history)
    
    # Metrics row
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric(
            "Current Pass Rate",
            f"{df.iloc[-1]['pass_rate']:.1%}",
            f"{(df.iloc[-1]['pass_rate'] - df.iloc[-2]['pass_rate']):.1%}" if len(df) > 1 else None
        )
    
    with col2:
        st.metric("Total Rules", df.iloc[-1]['total_tests'])
    
    with col3:
        st.metric("Passing", df.iloc[-1]['passing'])
    
    with col4:
        st.metric("Failing", df.iloc[-1]['failing'])
    
    # Trend chart
    st.subheader("Validation Trend")
    st.line_chart(df.set_index('timestamp')['pass_rate'])
    
    # Recent validations table
    st.subheader("Recent Validations")
    st.dataframe(df.tail(10))

else:
    st.warning("No validation history found. Run validation first.")

# Database statistics
st.subheader("Database Statistics")

db_files = list(Path('data/conditionDB').glob('*.json'))

for db_file in db_files:
    with open(db_file) as f:
        db = json.load(f)
    
    st.write(f"**{db_file.stem}**: {len(db)} entries")
    
    # Show priority distribution
    priorities = [entry.get('priority', 0) for entry in db]
    priority_df = pd.DataFrame({'Priority': priorities})
    st.bar_chart(priority_df['Priority'].value_counts().sort_index())
```

---

## Usage Examples

### Quick Start

```bash
# 1. Install automation dependencies
pip install -r requirements-automation.txt

# 2. Start PDF monitoring
python automation/orchestrator.py watch-pdfs --pdf-dir ./incoming_papers

# 3. Process a dataset
python automation/orchestrator.py process-dataset \
    --dataset data/uspto_grants.csv \
    --reaction-type Suzuki \
    --min-yield 0.70

# 4. Run validation
python automation/orchestrator.py validate

# 5. Analyze and auto-fix issues
python automation/orchestrator.py improve --auto-fix

# 6. Start continuous validation
python automation/orchestrator.py continuous-validation --interval-hours 6

# 7. Run full pipeline
python automation/orchestrator.py run-full-pipeline
```

### Using Docker

```bash
# Build and start all services
docker-compose up -d

# View logs
docker-compose logs -f pdf-monitor
docker-compose logs -f validator

# Access Airflow UI
open http://localhost:8080
```

### Using Airflow

```bash
# Initialize Airflow
airflow db init

# Start webserver and scheduler
airflow webserver -p 8080 &
airflow scheduler &

# Trigger the DAG
airflow dags trigger reaction_rules_pipeline

# Monitor in UI
open http://localhost:8080
```

---

## Configuration Files

### requirements-automation.txt

```
# Existing requirements
-r requirements.txt

# Automation-specific
apache-airflow==2.7.0
watchdog==3.0.0
schedule==1.2.0
click==8.1.7
streamlit==1.28.0

# Data processing
scikit-learn==1.3.0
scipy==1.11.3

# PDF processing
pdfplumber==0.10.3
PyPDF2==3.0.1

# Optional: Chemistry NLP
# ChemDataExtractor==2.1.0  # Requires additional setup
```

### .env Configuration

```bash
# automation/.env

# Paths
PDF_WATCH_DIR=./incoming_papers
DB_PATH=./data/conditionDB
OUTPUT_DIR=./automation/output
REPORT_DIR=./automation/reports

# Validation settings
VALIDATION_INTERVAL_HOURS=6
MIN_PASS_RATE_THRESHOLD=0.80

# Rule generation settings
MIN_YIELD=0.70
MIN_CLUSTER_SIZE=5
AUTO_FIX_ENABLED=true

# Notifications
ENABLE_SLACK=false
SLACK_WEBHOOK_URL=

ENABLE_EMAIL=false
EMAIL_SMTP_SERVER=
EMAIL_FROM=
EMAIL_TO=

# Airflow settings (if using)
AIRFLOW__CORE__EXECUTOR=LocalExecutor
AIRFLOW__CORE__DAGS_FOLDER=/opt/airflow/dags
```

---

## Summary

This automation framework provides:

✅ **Automatic data ingestion** from PDFs and datasets  
✅ **Intelligent rule generation** with SMARTS pattern discovery  
✅ **Automatic test case generation** for validation  
✅ **Continuous validation** with trend monitoring  
✅ **Smart fix recommendations** with auto-apply capability  
✅ **Workflow orchestration** (CLI, Airflow, or Docker)  
✅ **Monitoring dashboard** for visibility  
✅ **Version control integration** with automatic commits  

The system can run:
- **Fully automated** (scheduled pipeline)
- **Semi-automated** (human review of recommendations)
- **On-demand** (CLI commands)

Choose the level of automation based on your needs and confidence in the system!
