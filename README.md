# DVOUG: Dynamic Variable-Order Unitig Genome Graph

## Overview
**DVOUG** (Dynamic Variable-Order Unitig Genome Graph) is a novel graph construction framework designed for robust genome assembly and DNA storage data reconstruction under low sequencing coverage or high error rates. By dynamically adjusting the k-mer size during graph extension, DVOUG enhances graph connectivity while maintaining structural accuracy, effectively addressing fragmentation and noise issues in traditional de Bruijn graph (DBG) approaches.

### Key Features
- **Dynamic k-mer Adjustment**: Utilizes variable-order k-mers to repair fragmented regions in low-coverage or noisy data.
- **Source-Aware Indexing**: Integrates FM-index and positional context to filter redundant k-mers and guide precise graph extension.
- **Two-Stage Extension**: Combines **Mid-Max** (high k-value) and **Min-Mid** (low k-value) extension phases for optimal balance between specificity and connectivity.

---

## Installation

### Prerequisites
- **C++ Compiler** (supporting C++11 or higher)
- **Python 3.8+** (for auxiliary scripts)
- **BCALM2** ([GitHub](https://github.com/GATB/bcalm)): For initial unitig-level graph construction.
- **kISS**  ([GitHub](https://github.com/jhhung/kISS)): For efficient k-mer indexing.
- **WFA2-lib**  ([GitHub](https://github.com/smarco/WFA)): the wavefront alignment algorithm.

### Steps
1. Clone the repository:
   ```bash
   git clone https://github.com/lzqzz/DVOUG.git
   cd DVOUG