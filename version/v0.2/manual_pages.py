# This file defines the various informational pages for the Automated Sequence Assembly Software's GUI.
# Each class in this file represents a different page within the application's user manual or help system.


import tkinter as tk


class AboutPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        # Define the content for each section
        overview_content = ("This user manual describes an automated sequence assembly pipeline designed specifically "
                            "for Whole Exome Sequencing (WES). The pipeline aims to streamline the sequence assembly "
                            "process, providing a comprehensive, user-friendly solution that manages various stages of "
                            "data processing and analysis. By automating sequence assembly, the pipeline significantly "
                            "reduces the time and effort required compared to traditional manual methods, enabling "
                            "users to focus on other critical aspects of their work.")

        purpose_content = (
            "The primary purpose of this pipeline is to efficiently manage the sequence assembly process "
            "for WES data. Traditional manual sequence assembly methods are often time-consuming and "
            "labor-intensive, requiring extensive manual intervention and attention to detail. This "
            "automated pipeline addresses these challenges by providing a cohesive system that integrates "
            "multiple bioinformatics tools and processes, facilitating a smooth and accurate assembly of "
            "sequences. It allows users to input raw sequencing data and proceed through quality control, "
            "alignment, variant calling, and annotation with minimal manual intervention, ultimately saving "
            "valuable time and resources.")

        audience_content = (
            "This manual is intended for bioinformaticians, researchers, and lab technicians involved in "
            "genomic data analysis, particularly those working with WES data. The pipeline is designed to "
            "be accessible to users with a basic familiarity with genomic data and fundamental command-line "
            "usage. It provides a straightforward approach to handling the complexities of sequence assembly, "
            "making it suitable for a broad range of users within the genomic research community.")

        # Create frames and text widgets for each section
        self.create_section("Overview", overview_content)
        self.create_section("Purpose", purpose_content)
        self.create_section("Audience", audience_content)

    def create_section(self, title, content):
        frame = tk.Frame(self, borderwidth=2, relief="groove")
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        label = tk.Label(frame, text=title, font=("Arial", 16, "bold"))
        label.pack(pady=(5, 10))

        text_widget = tk.Text(frame, wrap="word", height=10, width=80, state="disabled", bg="#f0f0f0", bd=0)
        text_widget.pack(padx=10, pady=10, fill="both", expand=True)
        text_widget.configure(state="normal")
        text_widget.insert("1.0", content)
        text_widget.configure(state="disabled")


class InstallationPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        # Software Requirement content
        software_requirements = """
• Python 3.11.7:
  - Source: https://www.python.org/downloads/
  - Packages required:
    • os
    • json
    • threading
    • subprocess
    • tkinter (including ttk, filedialog, and scrolledtext)
    These packages are essential for ensuring the proper functioning of the pipeline.

• OpenJDK 17.0.12:
  - Source: https://openjdk.org/install/

• FastQC 0.12.1 (Andrews, 2010):
  - Source: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

• GATK 4.6.0.0 (McKenna et al., 2010):
  - Source: https://github.com/broadinstitute/gatk/releases

• Samtools 1.13 (Danecek et al., 2021):
  - Source: https://www.htslib.org/

• Burrows-Wheeler Aligner 0.7.17 (Li & Durbin, 2009):
  - Source: https://bio-bwa.sourceforge.net/

• AVVOVAR (Wang et al., 2010):
  - Version Date: 2020-06-07
  - Source: https://annovar.openbioinformatics.org/en/latest/user-guide/download/#annovar-main-package
"""

        # Hardware Requirement content
        hardware_requirements = """
• Operating system:
  - Ubuntu 22.04.4 LTS

• RAM:
  - 32.0 GiB 

• Bits:
  - 64-bit

• Hard Disk:
  - 1.0 TB

• Swap Space:
  - 50 GiB 

Notice: The hardware requirements shown are for efficiently operating the automated sequence assembly pipeline when processing sequence reads with a size of around 2 GiB. If  sequence files larger than 2 GiB are processed, it may be necessary to increase the RAM, Hard Disk, and Swap Space to ensure that the pipeline runs efficiently and without interruptions. Larger datasets can significantly increase the computational and storage demands of the pipeline, necessitating more robust hardware specifications.

"""

        # Create sections for each category
        self.create_section("Software Requirements", software_requirements)
        self.create_section("Hardware Requirements", hardware_requirements)

    def create_section(self, title, content):
        frame = tk.Frame(self, borderwidth=2, relief="groove")
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        label = tk.Label(frame, text=title, font=("Arial", 16, "bold"))
        label.pack(pady=(5, 10))

        text_widget = tk.Text(frame, wrap="word", height=10, width=80, state="disabled", bg="#f0f0f0", bd=0)
        text_widget.pack(padx=10, pady=10, fill="both", expand=True)
        text_widget.configure(state="normal")
        text_widget.insert("1.0", content)
        text_widget.configure(state="disabled")


class GettingStartedPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        # Launching the Application content
        launching_content = """
• To start the automated sequence assembly pipeline, you can use either a terminal or an integrated development environment (IDE).
  - Open the main.py in the IDE and run it.
  - Or, navigate to the directory containing the main.py file and run the following command in a terminal:
    \n    python main.py
  \nThis command launches the application, initiating the graphical user interface (GUI) for the pipeline. Ensure that all necessary dependencies and environment variables are set up correctly before launching.
"""

        # User Interface Overview content
        ui_overview_content = """
• MainPage:
  - This is the central hub where users can input their sequencing data. You can provide the paths to your sequencing output files in FASTQ format and the reference sequence. The MainPage also includes a "Run" button to initiate the sequence assembly process. This page is designed for quick access and easy initiation of the pipeline.

• ConfigurationPage:
  - The ConfigurationPage allows users to fine-tune the settings for various tools used in the pipeline. Here, you can adjust the arguments and parameters for tools such as GATK, BWA, and ANNOVAR. This customization enables more precise control over the analysis and ensures that the pipeline can be tailored to specific research needs.

• ProgramPage:
  - On the ProgramPage, users can provide the paths to the executable files of the external tools required by the pipeline. This includes specifying the locations for FASTQC, GATK, BWA, and ANNOVAR. Correctly setting these paths ensures that the pipeline can correctly call and utilize these tools during the analysis process.

• HelpPage:
  - The HelpPage provides guidance and support resources. It contains information about the application and troubleshooting. This page is designed to assist users in navigating the application and resolving common issues.
"""

        # Create sections for each category
        self.create_section("Launching the Application", launching_content)
        self.create_section("User Interface Overview", ui_overview_content)

    def create_section(self, title, content):
        frame = tk.Frame(self, borderwidth=2, relief="groove")
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        label = tk.Label(frame, text=title, font=("Arial", 16, "bold"))
        label.pack(pady=(5, 10))

        text_widget = tk.Text(frame, wrap="word", height=10, width=80, state="disabled", bg="#f0f0f0", bd=0)
        text_widget.pack(padx=10, pady=10, fill="both", expand=True)
        text_widget.configure(state="normal")
        text_widget.insert("1.0", content)
        text_widget.configure(state="disabled")


class WorkflowPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        workflow_content = """
1. Initial Setup and Configuration:
   - Install necessary toolkits and programs including Python, OpenJDK, FastQC, GATK, Samtools, BWA, and ANNOVAR.
   - Set up a directory for raw sequence files (FASTQ format) and reference genome (FASTA format).
   - Configure the pipeline with parameters and paths for the tools and input data.

2. Quality Control:
   - Tool: FastQC
   - Process: Assess the quality of raw sequencing data.
   - Outcome: Generates metrics such as quality scores, GC content distribution, and sequence length distribution to identify any irregularities that could impact downstream analyses.

3. Preprocessing:
   - Tool: Picard (FastqToSam, MarkIlluminaAdapters, SamToFastq)
   - Process:
     - Convert raw FASTQ sequences to BAM format, adding read group information.
     - Mark adapter sequences in BAM files to avoid misinterpretation.
     - Convert BAM files back to FASTQ format for alignment.
   - Outcome: Generates a BAM file with marked adapters, ready for alignment.

4. Sequence Alignment:
   - Tool: BWA (BWA-MEM), Picard (MergeBamAlignment, CreateSequenceDictionary, SortSam), Samtools
   - Process:
     - Index the reference genome.
     - Align sequencing reads to the reference genome using BWA-MEM.
     - Merge alignment data with unmapped BAM files.
     - Sort and prepare BAM files for downstream analysis.
   - Outcome: A merged and sorted BAM file containing aligned reads.

5. Post-Alignment Processing:
   - Tool: Picard (SetNmMdAndUqTags, MarkDuplicates, SortSam), GATK (BaseRecalibrator, ApplyBQSR)
   - Process:
     - Calculate and add crucial tags (NM, MD, UQ) to BAM files.
     - Mark duplicate reads to reduce redundancy.
     - Recalibrate base quality scores using BaseRecalibrator and ApplyBQSR.
   - Outcome: A high-quality, duplicate-marked, and recalibrated BAM file.

6. Variant Calling:
   - Tool: GATK (HaplotypeCaller, SelectVariants)
   - Process:
     - Call variants (SNVs and Indels) from the recalibrated BAM file.
     - Extract specific variant types (SNVs and Indels) for further analysis.
   - Outcome: VCF files containing detailed information about detected variants.

7. Quality Evaluation:
   - Tool: Picard (CollectVariantCallingMetrics)
   - Process: Collect metrics to evaluate the quality of variant calls, including comparison with known variant databases.
   - Outcome: Detailed metrics and statistics that assess the accuracy and quality of the variant calls.

8. Variant Annotation:
   - Tool: ANNOVAR
   - Process: Annotate variants with additional information such as gene names, predicted effects, allele frequencies, and clinical significance.
   - Outcome: Annotated VCF files that provide comprehensive insights into the biological and clinical relevance of the identified variants.
"""

        self.create_section("Workflow Overview", workflow_content)

    def create_section(self, title, content):
        frame = tk.Frame(self, borderwidth=2, relief="groove")
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        label = tk.Label(frame, text=title, font=("Arial", 16, "bold"))
        label.pack(pady=(5, 10))

        text_widget = tk.Text(frame, wrap="word", height=10, width=80, state="disabled", bg="#f0f0f0", bd=0)
        text_widget.pack(padx=10, pady=10, fill="both", expand=True)
        text_widget.configure(state="normal")
        text_widget.insert("1.0", content)
        text_widget.configure(state="disabled")


class IssuesPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        # Common Issues content
        common_issues_content = """
• Input File Format Errors:
  - Issue: The pipeline does not recognize the input files.
  - Solution: Ensure that the input files are in the correct format (FASTQ for sequencing reads and FASTA for the reference genome). Verify that the files are properly named and located in the specified input directory.

• Insufficient Resources:
  - Issue: The pipeline crashes or runs very slowly.
  - Solution: Check that your system meets the minimum hardware requirements, including RAM and disk space. Consider increasing swap space if memory is insufficient. Ensure that your system has enough free disk space for temporary files and output data.

• Tool Path Configuration:
  - Issue: The pipeline cannot find external tools (e.g., GATK, BWA).
  - Solution: Verify that the paths to external tools are correctly set in the ProgramPage. Ensure that the specified paths are accurate and point to the correct executable files.

• Missing Required Arguments:
  - Issue: The pipeline does not start or produces errors due to missing required arguments.
  - Solution: Ensure all required fields are filled in the ConfigurationPage, such as input files, output directories, and reference files.

• Known Site and Reference Sequence Mismatch:
  - Issue: The known sites file does not match the reference sequence being used.
  - Solution: Ensure that the known sites file (e.g., dbSNP, 1000 Genomes) corresponds to the same genome build as the reference sequence. Using known sites from a different build can result in errors or incorrect variant calling.

• Invalid Directory Names:
  - Issue: The directory name contains spaces or special characters, causing issues with file retrieval and output.
  - Solution: Ensure that directory names used in the pipeline do not contain spaces or special characters. Use underscores (_) or hyphens (-) instead of spaces, and avoid characters like &, %, $, etc.
"""

        # Logs and Debugging content
        logs_debugging_content = """
• Accessing Logs: All logs and error messages are recorded in the "log" folder, located within the same directory as the pipeline. This folder contains detailed logs of each step, including command executions, outputs, and error messages.

• Interpreting Logs: Logs provide information on the progress and status of each step in the pipeline. They can be used to identify where the pipeline encountered issues and provide details about what went wrong. Look for lines marked as "ERROR" or "WARNING" to quickly find potential problems.

• Important Note: The log files in the "log" folder will be overwritten each time a new sequence assembly is run. If you need to retain logs from a specific analysis, please save them to another directory before starting a new run. This ensures that you have a permanent record of the logs for future reference or troubleshooting. You can copy the log files to a designated backup location or rename them to include relevant details about the analysis they correspond to.
"""

        # Create sections for each category
        self.create_section("Common Issues", common_issues_content)
        self.create_section("Logs and Debugging", logs_debugging_content)

    def create_section(self, title, content):
        frame = tk.Frame(self, borderwidth=2, relief="groove")
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        label = tk.Label(frame, text=title, font=("Arial", 16, "bold"))
        label.pack(pady=(5, 10))

        text_widget = tk.Text(frame, wrap="word", height=10, width=80, state="disabled", bg="#f0f0f0", bd=0)
        text_widget.pack(padx=10, pady=10, fill="both", expand=True)
        text_widget.configure(state="normal")
        text_widget.insert("1.0", content)
        text_widget.configure(state="disabled")
