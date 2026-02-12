# EHT Pacing-Following Processing Pipeline (Script 01/03)

This repository contains an R script for preprocessing engineered heart tissue (EHT) measurements generated using the **NovoHeart CTScreen** platform.

This script represents the first step (01/03) of a modular EHT analysis pipeline.

The pipeline combines CTScreen output files, extracts pacing step information, computes expected pacing frequencies, and determines whether tissues are following external pacing 
stimulation based on a defined tolerance threshold.

---

## What the pipeline does

Starting from raw CTScreen **Summary_Data.csv** files, the script:

  1. Loads and merges multiple Summary_Data CSV files (recursively)
  2. Automatically detects file delimiter and header structure
  3. Extracts pacing step information from file names (e.g. step_001, step_002)
  4. Maps expected pacing frequency for each step
  5. Adds experimental metadata:
        - Batch
        - Donor
        - Week_Time
  6. Computes:
        - Deviation from expected pacing frequency
        - Relative error
  7. Classifies each recording as Following or Not Following based on a user-defined tolerance (e.g., 5%, 10%, 15%)
  8. Writes a timestamped combined dataset to a structured output folder

All file paths and parameters are selected interactively via GUI dialogs.

---

## Required inputs

1. NovoHeart CTScreen Summary files

  - One or multiple Summary_Data.csv files
  - Files may be stored in nested subfolders (recursive search supported)
  - Each file must contain:
          - File_Name column
          - Frequency_Hz column
          - The file name must contain a Pacing Step

2. Experimental metadata (user-defined via GUI)
   
  - The user provides:
          - Patient_ID (Donor)
          - Week_Time
          - Batch (yymmdd format)
          - Following tolerance threshold

---

## Pacing Step Mapping

The script automatically maps pacing step codes to expected frequencies:

  Step	      Expected_Frequency
  step_001	  0.0 Hz (spontaneous)
  step_002	  1.0 Hz
  step_003	  1.5 Hz
  step_004	  2.0 Hz
  step_005	  2.5 Hz
  step_006	  3.0 Hz

For spontaneous recordings (step_001), Following is defined as:
- Yes → if measured frequency = 0 Hz
- No → otherwise

For paced recordings, Following is defined as:
- abs(Deviation) ≤ (Tolerance × Expected_Hz)

---

## Cleaned summary output

The script generates a timestamped output folder:

_EHT_Analysis_YYYY-MM-DD_HH-MM/
    Combined_Summary_with_Following_<tag>_tolXXpct.csv

This dataset contains:
- Original measurement values
- Metadata (Batch, Donor, Week_Time)
- Expected_Hz
- Deviation_Hz
- Rel_Error
- Following classification

The output is ready for:
- Statistical analysis
- Donor-level aggregation
- Visualization
- Downstream merging (Script 02)

---

## Typical use cases

1. EHT pacing-following validation
2. Drug response studies (pre vs post treatment)
3. Patient-specific hiPSC-derived EHT comparison
4. Batch-controlled contractility experiments
5. Standardized preprocessing in collaborative projects

---

## Position in the EHT pipeline

This script is Script 01 of a 3-step EHT analysis workflow:

Script 01 (this repository) – Combined dataset generation and Following classification
Script 02 – Donor merging and filtering (baseline + Following-only selection)
Script 03 - 

---

## Methods Description

EHT contractility measurements generated using the NovoHeart CTScreen platform were combined and annotated using a custom R-based GUI pipeline. Pacing step information was extracted from file names and mapped to expected frequencies. Deviation from expected pacing was calculated, and recordings were classified as Following based on a predefined tolerance threshold. The resulting structured dataset was exported for downstream statistical analysis.

---

## Authorship

This script was developed by Michele Buono and can be used freely for research purposes, provided appropriate citation of the author.
The overall workflow, structure, and clarity of the pipeline were iteratively refined with assistance from OpenAI - ChatGPT 5.2, which was used as a tool to improve code 
organization, documentation, and usability.
