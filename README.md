# Fast-decoding-of-group-testing-results-from-Reed-Solomon-codes

This repository contains the Java code used in the research paper titled **"Fast Decoding of Group Testing Results from Reed-Solomon d-Disjunct Matrices"** by **Dongxia (Mico) Luo** and **Lucia Moura**. 

The paper has been accepted for presentation at the **International Workshop on the Arithmetic of Finite Fields (WAIFI 2024)** and will be published in **Lecture Notes in Computer Science (LNCS)**.

## Paper Abstract

Non-adaptive combinatorial group testing has applications in disease screening as well as in many problems in digital security and communications. Matrices that are d-disjunct (also called d-cover-free) can be built using codes and allow for the detection of d defective items using group testing. In this paper, we study d-disjunct matrices built from Reed-Solomon codes, and design a specialized algorithm for decod- ing the results of group testing using these matrices. We do an exper- imental comparison between our method and the naive one that only uses the d-disjunct property of the matrix, and show that the former outperforms the latter as the size of the problem grows.

## Code Overview

This repository includes two Java files that implement the algorithms described in the paper:

- **`FindDefectives.java`**: Implements the fast decoding algorithm for identifying defective items using the Reed-Solomon d-disjunct matrices.
- **`ModInverse.java`**: Implements the modular inverse matrix algorithm used in the construction and decoding process.

Both files are designed to work together and are optimized for the scenarios described in the paper.

## Requirements

No additional libraries are required to run the provided Java code.

## Usage

To compile and run the code:

1. Download the repository.
2. Compile the Java files using the following command:
   ```bash
   javac FindDefectives.java ModInverse.java
