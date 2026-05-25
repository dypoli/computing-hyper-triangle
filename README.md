Source code for the paper Efficient Computation of Hyper-triangles on Hypergraphs

The real-world datasets used in the paper are available: https://www.cs.cornell.edu/~arb/data/

ApproOptimal.cpp corresponding to the algorithm: Appro-adv
CountDenseTTT.cpp corresponding to the algorithm: CountDenseTTT
CountSpareTTT.cpp corresponding to the algorithm: CountSpareTTT
ExactAlgorithm.cpp corresponding to the algorithm: CountTTT   CountTTC   CountTCC   CountCCC  Exact-bs
ParOptimal.cpp corresponding to the algorithm: Par-adv

For the algorithms: Exact-bs, Par-bs, Appro-bs, we simply modify the algorithms: MoCHy-E,  parallelized MoCHy-E, MoCHy-A+, of the paper "Hypergraph Motifs: Concepts, Algorithms, and Discoveries"


## Dynamic Algorithms

The `Dynamic` folder contains four main directories, each implementing different update operations for hyperedge triangle counting:

### 1. AddHyperedge

**Location:** `Dynamic/AddHyperedge/`

**Description:**
This module implements algorithms for computing hyper-triangles when new hyperedges are dynamically added to the hypergraph.

**Key Files:**
- `HIU.cpp` - Basic Hyperedge Insertion Update algorithm
- `HIU+.cpp` - Optimized version of the Hyperedge Insertion Update algorithm
- `HelperFunc.h` - Utility functions and data structure definitions
- `read_data.cpp` - Data input/parsing functions
- `NewHyperedge.txt` - Test data file containing new hyperedges to be inserted
- `TriangleNumber.txt` - input file for hyper-triangle counts before inserting
- `unique-email-Enron.txt` - Enron email dataset used for testing

**Algorithm Overview:**
- Reads the initial hypergraph structure from the Enron dataset
- Dynamically inserts new hyperedges from `NewHyperedge.txt`
- For each new hyperedge, identifies neighboring hyperedges (those sharing vertices)
- Finds newly formed hyper-triangles (three pairwise adjacent hyperedges)
- Maintains and updates the hyper-triangle count by pattern type
- Enforces ordering constraints to avoid duplicate counting

**Main Operations:**
1. Load original hypergraph and build adjacency structures
2. For each new hyperedge:
   - Update the hyperedge-to-node mapping
   - Find neighboring hyperedges
   - Update reverse adjacency links
   - Identify and count newly formed triangles

---

### 2. AddVertex

**Location:** `Dynamic/AddVertex/`

**Description:**
This module implements algorithms for computing hyper-triangles when new vertices are dynamically added to existing hyperedges in the hypergraph.

**Key Files:**
- `VIU.cpp` - Basic Vertex Insertion Update algorithm
- `VIU+.cpp` - Optimized version of the Vertex Insertion Update algorithm
- `HelperFunc.h` - Utility functions and pattern detection
- `read_data.cpp` - Data loading functions
- `NewVertices.txt` - File specifying which vertices to add to which hyperedges 
- `MarkTable.txt` - Table marking affected hyperedges as decribed in the paper
- `PreRecordedTable.txt` - Pre-recorded table as decribed in the paper
- `TriangleNumber.txt` - input file for hyper-triangle counts before inserting
- `unique-email-Enron.txt` - Enron email dataset

**Algorithm Overview:**
- Loads the initial hypergraph structure
- Reads vertex additions from `NewVertices.txt` (specifying which vertices to add to which hyperedges)
- For each affected hyperedge, tracks structural changes
- Classifies affected hyperedges into three categories:
  - Old adjacencies (unchanged neighbors)
  - Adjusted adjacencies (neighbors with modified vertex sets)
  - New adjacencies (freshly formed connections)
- Updates hyper-triangle counts by analyzing:
  - Triangles with old neighbors
  - Triangles with adjusted neighbors
  - Triangles involving newly formed connections

**Main Operations:**
1. Load original hypergraph and identify affected hyperedges
2. Build three adjacency lists (old, adjust, new)
3. Update triangles by:
   - Decrementing counts for old triangle patterns
   - Incrementing counts for new triangle patterns after vertex addition

---

### 3. DeleteHyperedge

**Location:** `Dynamic/DeleteHyperedge/`

**Description:**
This module implements algorithms for computing hyper-triangles when hyperedges are dynamically deleted from the hypergraph.

**Key Files:**
- `HDU.cpp` - Basic Hyperedge Deletion Update algorithm
- `HDU+.cpp` - Optimized version of the Hyperedge Deletion Update algorithm
- `HelperFunc.h` - Utility functions for deletion operations
- `read_data.cpp` - Data input functions
- `DeleteHyperedge.txt` - File containing IDs of hyperedges to be deleted
- `TriangleNumber.txt` - input file for hyper-triangle counts before deleting
- `unique-email-Enron.txt` - Enron email dataset

**Algorithm Overview:**
- Reads the initial hypergraph structure
- Loads the list of hyperedges to delete from `DeleteHyperedge.txt`
- Processes each deletion sequentially
- For each deletion:
  - Identifies all hyper-triangles containing the hyperedge to be deleted
  - Decrements the count for each affected triangle pattern
  - Removes the hyperedge from all adjacency structures
  - Cleans up node-to-hyperedge mappings

**Main Operations:**
1. Load initial hypergraph and adjacency lists
2. For each hyperedge to delete:
   - Find all triangles including this hyperedge (before deletion)
   - Decrement triangle counts
   - Remove from neighbors' adjacency lists
   - Remove from node-to-hyperedge mappings
   - Mark hyperedge as deleted

---

### 4. DeleteVertex

**Location:** `Dynamic/DeleteVertex/`

**Description:**
This module implements algorithms for computing hyper-triangles when vertices are dynamically deleted from hyperedges in the hypergraph.

**Key Files:**
- `VDU.cpp` - Vertex Deletion Update algorithm
- `HelperFunc.h` - Utility functions for vertex deletion
- `read_data.cpp` - Data loading functions
- `DeleteVertices.txt` - File specifying which vertices to delete
- `TriangleNumber.txt` - input file for hyper-triangle counts before deleting
- `unique-email-Enron.txt` - Enron email dataset

**Algorithm Overview:**
- Loads the initial hypergraph structure
- Reads vertices to be deleted from `DeleteVertices.txt`
- Identifies all hyperedges affected by vertex deletions
- Classifies affected hyperedges based on adjacency changes:
  - Old adjacencies (unaffected neighbors)
  - Adjusted adjacencies (neighbors whose vertex sets changed)
  - Deleted adjacencies (lost connections)
- Updates hyper-triangle counts by:
  - Finding and decrementing counts for triangles affected by the deletion
  - Properly handling cases where both hyperedges in a triangle are partially modified

**Main Operations:**
1. Load initial hypergraph
2. Remove deleted vertices from all hyperedges
3. Identify completely deleted hyperedges
4. Build updated adjacency structures
5. Update triangle counts by analyzing:
   - Triangles with unmodified neighbors (old)
   - Triangles with modified neighbors (adjust)
   - Triangles involving deleted connections (delete)

---

## Data Format

### Input Files

**Hypergraph Format (unique-email-Enron.txt):**
```
vertex1 vertex2 vertex3 ... vertexN
```
Each line represents one hyperedge containing the listed vertices (space-separated integers).

**New Hyperedges/Vertices (NewHyperedge.txt, NewVertices.txt):**
Similar format, specifying which vertices belong to new/modified hyperedges.

**Deletions (DeleteHyperedge.txt, DeleteVertices.txt):**
List of vertex or hyperedge IDs to be deleted, one per line.

### Output Files

**TriangleNumber.txt:**
Contains the count of each hyper-triangle pattern type. The file has multiple integer values representing different patterns of 3-vertex subsets within the three hyperedges.

---

## Hyper-Triangle Definition

A hyper-triangle consists of three pairwise adjacent hyperedges:
- Hyperedges `e1`, `e2`, `e3` are adjacent if they share at least one vertex
- A hyper-triangle is defined by the pattern of shared vertices between the three hyperedges
- Different patterns are counted separately (e.g., patterns distinguishing different vertex overlap configurations)

---

## Building and Running

Each folder contains compiled C++ code designed to:
1. Read the initial hypergraph
2. Load update specifications (additions or deletions)
3. Compute the impact on hyper-triangle counts
4. Output results to `TriangleNumber.txt`

The algorithms use:
- Vector-based adjacency lists
- Hash maps for fast lookup
- Binary search for membership testing
- Two-pointer technique for set intersection operations

---

## Performance Characteristics

The algorithms are designed for incremental updates with complexity proportional to:
- The number of neighbors of affected hyperedges
- The overlap between adjacency lists
- The complexity of pattern detection (constant for each potential triangle)

The optimized versions (HIU+, VIU+, HDU+) implement additional data structures and caching mechanisms to improve practical performance on large dynamic hypergraphs.

---

## Dataset

All modules use the **Enron email dataset** as their primary test case:
- **Source:** Enron email communication network
- **Description:** Each hyperedge represents an email with multiple recipients
- **File:** `unique-email-Enron.txt`
- **Purpose:** Real-world benchmark for dynamic hyper-triangle counting algorithms
