Mutation Picker Tool
Introduction:

The Mutation Picker Tool is designed to facilitate identifying and analyzing protein mutations across different generations and conditions. This tool automates the process of comparing protein sequences, generating mutation data, and aiding in the analysis of protein variations.

Usage:

Running the Script:
The core functionality of the tool is encapsulated within the main.py script. Upon execution, this script opens a dialogue to input essential parameters and settings. Users are required to provide input values in the designated fields before initiating the analysis.

Input Fields:

Input File: The tool requires an aligned protein FASTA file to perform a thorough analysis. This file should contain the protein sequences from different generations consolidated together. This alignment process ensures that similar proteins from distinct strains are positioned together, enabling accurate mutation tracking within specific proteins across organisms under varying conditions.

Path: Path of the input file
Wild ID: This field corresponds to the identification of the wild-type protein.

Generations: Specify the number of protein generations for comparison and analysis.
Example Scenario:

Consider the case of studying the proteome of an E. coli strain ("ecoliT") in various conditions. The goal is to compare this wild-type protein against six treated E. coli strains under different conditions, spanning seven generations. The use of an aligned FASTA file is imperative for this analysis. Aligning protein sequences enables the consolidation of similar proteins from distinct strains, facilitating mutation tracking within specific proteins across organisms under differing conditions.

In this example, the "Wild ID" would be set as "ecoliT," this should be consistently reflected in the protein names provided within the aligned FASTA file.

Output:

Upon input validation, the script imports the task.py module. This module then performs a comprehensive analysis, generating three distinct CSV files. These files contain the following information:

Mutations between different generations of the same protein.
Truncations in protein sequences.
Additional truncations are specific to the analyzed protein.
All generated files are saved within the same directory.

Conclusion:

The Mutation Picker Tool is invaluable for researchers aiming to analyze and compare protein mutations across generations and conditions. This tool accelerates mutation tracking and facilitates informed biological insights by automating the analysis process.

For further information or troubleshooting, please refer to the tool's documentation or contact our support team.
