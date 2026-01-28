"""
ARCHITECTURAL ISSUE: Detection System Not Using reacted_motifs

USER OBSERVATION:
-----------------
The system correctly identifies:
  reacted_motifs: ['Ar-B(OH)2', 'Ar-Br']
  formed_motifs: ['Ar-Ar']

This is clearly a Suzuki pattern, so why isn't it used for detection?

ROOT CAUSE:
-----------
Order of operations in featurize_reaction():

1. detect_reaction_types() called FIRST (line ~212)
   - Uses ALL motifs in reactants
   - Sees Ar-B(OH)2, Ar-Br, AND Ar-H (from boronic acid)
   - Matches Arylation_Ar_H rule (electrophile + Ar-H substrate)

2. Featurize reactants and products (line ~233-250)
   
3. Calculate aggregates LATER (line ~260)
   - Compares reactant vs product motifs
   - Correctly identifies reacted_motifs: ['Ar-B(OH)2', 'Ar-Br']
   
4. Generate reaction_key (line ~303)
   - Uses reacted_motifs for key generation

The detection system CAN'T use reacted_motifs because they haven't been
calculated yet!

WHY THIS DESIGN:
----------------
Historical reasons:
- Detection was designed to work on raw reaction SMILES only
- Product confirmation is optional (some reactions don't have products)
- Early detection helps filter/route reactions before expensive featurization

PROPOSED SOLUTIONS:
-------------------

Option 1: TWO-PASS DETECTION (Recommended)
- First pass: Current slot-based detection
- Second pass: Use reacted_motifs to validate/refine
- If reacted_motifs contradicts first pass, prefer second pass result
- Example: If reacted=['Ar-B(OH)2', 'Ar-Br'], prioritize Suzuki over Arylation

Option 2: REORDER OPERATIONS
- Calculate aggregates BEFORE detection
- Requires products (not always available)
- Detection uses reacted_motifs as primary evidence
- Fallback to slot-based if products unavailable

Option 3: IMPROVED SLOT MATCHING
- Add exclusion rules to reaction definitions
- Example: Arylation_Ar_H excludes molecules with @organoboron
- Requires updating all reaction rules in taxonomy
- More maintainable than detection logic changes

Option 4: USE ROLE CLASSIFICATION (Current Workaround)
- Role classification uses motif analysis + slot semantics
- Already correctly identifies Suzuki_miyaura
- CLI now shows "Role-Based Reaction Type" prominently
- Recommendation: Use role classification for conditions matching

RECOMMENDED APPROACH:
---------------------
Immediate: Use role classification result (already implemented in CLI)
Short-term: Add two-pass detection with reacted_motifs validation
Long-term: Redesign detection to use product analysis by default

CODE CHANGES NEEDED:
--------------------
1. In chemtools/featurizers/formatters/reaction.py:
   - Add validate_detection_with_aggregates() function
   - Call after aggregates calculated
   - Update reaction_type if reacted_motifs suggest different family

2. In chemtools/reaction_type_detection.py:
   - Add reacted_motifs_based_detection() function
   - Priority matching: Ar-B(OH)2 + Ar-Br = Suzuki, etc.

3. Update CLI to show both:
   - Slot-based detection (current)
   - Reacted-motifs detection (validated)
   - Role classification (most accurate)

IMPACT ASSESSMENT:
------------------
This affects ALL cross-coupling reactions where one reactant contains
spectator Ar-H bonds:
- Suzuki (boronic acids have Ar-H)
- Stille (organostannanes may have Ar-H)
- Negishi (organozincs may have Ar-H)
- Heck (alkenes with Ar may have Ar-H)

Estimated false positives: 10-20% of cross-coupling reactions
