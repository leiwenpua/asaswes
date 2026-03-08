# ASASWES: Automated Sequence Assembly System for WES

[![Python Version](https://img.shields.io/badge/python-3.11.7-blue.svg)](https://www.python.org/)
[![OS](https://img.shields.io/badge/OS-Ubuntu_22.04-orange.svg)](https://ubuntu.com/)
[![Pipeline](https://img.shields.io/badge/Workflow-GATK%20|%20BWA%20|%20Samtools-green.svg)]()

## Overview
**ASASWES** is an end-to-end, automated pipeline designed to streamline the **Whole Exome Sequencing (WES)** analysis workflow. It transforms raw sequencing data into annotated variant reports with minimal manual intervention.

Developed with a focus on **efficiency and reproducibility**, this system replaces labor-intensive manual command-line execution with a cohesive, modular framework, featuring a real-time monitoring GUI for enhanced process management.

---

## Key Features
* **End-to-End Automation:** Integrated QC, Alignment, Post-processing, Variant Calling, and Annotation.
* **Customizable ETL Logic:** Granular control over parameters for GATK, BWA, and Picard tools via a dedicated Configuration UI.
* **Real-time Monitoring:** Graphical widget to track command execution and logs in real-time.
* **Standardized Output:** Ensures consistent data structures and clinical-grade report readiness.

---

## Pipeline Workflow
The system orchestrates the following bioinformatics tools:
1.  **Quality Control:** `FastQC` (Raw data assessment)
2.  **Data Transformation (Picard):** `FastqToSam` -> `MarkIlluminaAdapters` -> `SamToFastq`
3.  **Alignment:** `BWA MEM` (Burrows-Wheeler Aligner)
4.  **Post-Alignment (Picard):** `MergeBamAlignment` -> `SortSam` -> `SetNmMdAndUqTags`
5.  **Variant Calling:** `GATK 4.6.0.0`
6.  **Annotation:** `ANNOVAR` 

---

## System Requirements & Dependencies
### **Environment**
* **OS:** Ubuntu 22.04.4 LTS
* **Hardware:** 32GB RAM | 1TB Disk Space | 64-bit Architecture
* **Runtime:** Python 3.11.7 | OpenJDK 17.0.12

### **Core Tools Integrated**
| Tool | Version | Purpose |
| :--- | :--- | :--- |
| **GATK** | 4.6.0.0 | Variant Discovery |
| **BWA** | 0.7.17 | Genome Mapping |
| **Samtools** | 1.13 | High-throughput sequencing data manipulation |
| **FastQC** | 0.12.1 | Sequence Quality Control |
| **ANNOVAR** | 2020-06-07 | Functional Annotation |

---

## Getting Started
### **Installation**
1. Clone the repository:
   ```bash
   git clone [https://github.com/leiwenpua/asaswes.git](https://github.com/leiwenpua/asaswes.git)
   cd asaswes

2. Ensure all dependencies (GATK, BWA, etc.) are installed and paths are configured in the ProgramPage within the GUI.

### **Usage**
1. Launch the application:
   ```bash
   python main.py

2. MainPage: Set output directory and load FASTQ/Reference files.
   <img width="728" height="508" alt="image" src="https://github.com/user-attachments/assets/81f9b76c-8f54-457f-9a44-32bc5d8fc350" />

3. ConfigurationPage: Fine-tune tool arguments (e.g., Sample Name, Platform).
   <img width="745" height="484" alt="image" src="https://github.com/user-attachments/assets/8c7a82b0-70dd-4d87-b40b-29cc611e85fd" />

4. ProgramPage: Map tool executable paths.
   <img width="784" height="500" alt="image" src="https://github.com/user-attachments/assets/869995d5-2893-437c-a7ca-bc7b7cd09cd8" />

5. Run: Execute and monitor the progress.

### **Documentation**
A comprehensive User Manual is available in this repository, providing detailed instructions on:

Advanced parameter tuning for Picard & GATK.

Troubleshooting and error handling.

GUI navigation and interface overview.
