import tkinter as tk
from tkinter import ttk, filedialog


class ConfigurationPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        # Create a canvas and a vertical scrollbar for scrolling it
        canvas = tk.Canvas(self)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)

        # Configure the scrollable frame to resize with the canvas
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        # Place the scrollable frame inside the canvas
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        # Pack the canvas and scrollbar into the frame
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        # Add a label to the scrollable frame
        label = tk.Label(scrollable_frame, text="Configuration Page", font=("Arial", 24))
        label.pack(pady=10)

        # FastqToSam (Picard) section
        self.create_picard_section(scrollable_frame, "FastqToSam (Picard)", [
            ("Required Arguments for FastqToSam", [
                ("FASTQ to SAM Output", "FastqToSam.bam", False, tk.StringVar),
                ("Sample Name", "sample001", False, tk.StringVar),
                ("Platform", "ILLUMINA", False, tk.StringVar),
            ]),
            ("Optional Arguments for FastqToSam", [
                ("Allow and Ignore Empty Lines", False, True, tk.BooleanVar),
                ("Use Sequential FASTQs", False, True, tk.BooleanVar),
                ("Comment", "", False, tk.StringVar),
                ("Description", "", False, tk.StringVar),
                ("Library Name", "", False, tk.StringVar),
                ("Max Q", 93, False, tk.IntVar),
                ("Min Q", 0, False, tk.IntVar),
                ("Platform Model", "", False, tk.StringVar),
                ("Platform Unit", "", False, tk.StringVar),
                ("Predicted Insert Size", 0, False, tk.IntVar),
                ("Program Group", "", False, tk.StringVar),
                ("Quality Format", "", False, tk.StringVar),
                ("Read Group Name", "A", False, tk.StringVar),
                ("Run Date", "", False, tk.StringVar),
                ("Sequencing Center", "", False, tk.StringVar),
                ("FastqToSam Sort Order", "queryname", False, tk.StringVar)
            ])
        ])

        # MarkIlluminaAdapters (Picard) section
        self.create_picard_section(scrollable_frame, "MarkIlluminaAdapters (Picard)", [
            ("Required Arguments for MarkIlluminaAdapters", [
                ("Mark Illumina Metrics", "metrics.txt", False, tk.StringVar),
                ("Mark Illumina Output", "MarkIlluminaAdapters_read_pairs.bam", False, tk.StringVar)
            ]),
            ("Optional Arguments for MarkIlluminaAdapters", [
                ("Adapter Truncation Length", 30, False, tk.IntVar),
                ("Adapters", "PAIRED_END", False, tk.StringVar),
                ("Five Prime Adapter", "", False, tk.StringVar),
                ("Max Error Rate PE", 0.1, False, tk.DoubleVar),
                ("Max Error Rate SE", 0.1, False, tk.DoubleVar),
                ("Min Match Bases PE", 6, False, tk.IntVar),
                ("Min Match Bases SE", 12, False, tk.IntVar),
                ("Num Adapters to Keep", 1, False, tk.IntVar),
                ("Three Prime Adapter", "", False, tk.StringVar)
            ])
        ])

        # SamToFastq section
        self.create_picard_section(scrollable_frame, "SamToFastq (Picard)", [
            ("Required Arguments for SamToFastq", [
                ("SAM to FASTQ Output", "SamToFastq.fastq", False, tk.StringVar)
            ]),
            ("Optional Arguments for SamToFastq", [
                ("Clipping Action", "2", False, tk.StringVar),
                ("Clipping Attribute", "XT", False, tk.StringVar),
                ("Clipping Min Length", 0, False, tk.IntVar),
                ("Include Non PF Reads", True, True, tk.BooleanVar),
                ("Include Non Primary Alignments", False, True, tk.BooleanVar),
                ("Interleave", True, True, tk.BooleanVar),
                ("Quality", 0, False, tk.IntVar),
                ("Read1 Max Bases to Write", 0, False, tk.IntVar),
                ("Read1 Trim", 0, False, tk.IntVar),
                ("Read2 Max Bases to Write", 0, False, tk.IntVar),
                ("Read2 Trim", 0, False, tk.IntVar),
                ("RG Tag", "PU", False, tk.StringVar)
            ])
        ])

        # BWA alignment section
        self.create_picard_section(scrollable_frame, "BWA Alignment", [
            ("Required Arguments for BWA Alignment", [
                ("Threads", 10, False, tk.IntVar),
                ("Index Algorithm", "bwtsw", False, tk.StringVar)
            ])
        ])

        # CreateSequenceDictionary section
        self.create_picard_section(scrollable_frame, "CreateSequenceDictionary (Picard)", [
            ("Optional Arguments for CreateSequenceDictionary", [
                ("Genome Assembly", "", False, tk.StringVar),
                ("Number of Sequences", 2147483647, False, tk.IntVar),  # Ensure this line exists
                ("Species", "", False, tk.StringVar),
                ("Truncate Names at Whitespace", True, True, tk.BooleanVar),
                ("URI", "", False, tk.StringVar)
            ])
        ])

        # MergeBamAlignment (Picard) section
        self.create_picard_section(scrollable_frame, "MergeBamAlignment (Picard)", [
            ("Required Arguments for MergeBamAlignment", [
                ("Merge Align BAM Output", "merge_align.bam", False, tk.StringVar)
            ]),
            ("Optional Arguments for MergeBamAlignment", [
                ("Add Mate Cigar", True, True, tk.BooleanVar),
                ("Aligned Reads Only", False, True, tk.BooleanVar),
                ("Aligner Proper Pair Flags", True, True, tk.BooleanVar),
                ("Attributes To Remove", "", False, tk.StringVar),
                ("Attributes To Retain", "X0,XS", False, tk.StringVar),
                ("Attributes To Reverse", "", False, tk.StringVar),
                ("Attributes To Reverse Complement", "", False, tk.StringVar),
                ("Clip Adapters", False, True, tk.BooleanVar),
                ("Clip Overlapping Reads", True, True, tk.BooleanVar),
                ("Expected Orientations", "FR", False, tk.StringVar),
                ("Include Secondary Alignments", True, True, tk.BooleanVar),
                ("MergeBamAlignment Is Bisulfite Sequence", False, True, tk.BooleanVar),
                ("Matching Dictionary Tags", "", False, tk.StringVar),
                ("Max Insertions Or Deletions", -1, False, tk.IntVar),
                ("Min Unclipped Bases", 32, False, tk.IntVar),
                ("Primary Alignment Strategy", "MostDistant", False, tk.StringVar),
                ("Program Group Command Line", "", False, tk.StringVar),
                ("Program Group Name", "", False, tk.StringVar),
                ("Program Group Version", "", False, tk.StringVar),
                ("Program Record Id", "", False, tk.StringVar),
                ("Read1 Aligned Bam", "", False, tk.StringVar),
                ("Read1 Trim", 0, False, tk.IntVar),
                ("Read2 Aligned Bam", "", False, tk.StringVar),
                ("Read2 Trim", 0, False, tk.IntVar),
                ("MergeBamAlignment Sort Order", "unsorted", False, tk.StringVar),
                ("Unmap Contaminant Reads", True, True, tk.BooleanVar),
                ("Unmapped Read Strategy", "COPY_TO_TAG", False, tk.StringVar)
            ])

        ])

        # SortSam (Picard) For Merged Bam section
        self.create_picard_section(scrollable_frame, "SortSam (Picard) For Merged BAM", [
            ("Required Arguments for SortSam", [
                ("SortSam For Merged BAM Output", "sorted_merge.bam", False, tk.StringVar),
                ("SortSam For Merged BAM Sort Order", "coordinate", False, tk.StringVar)
            ]),
            ("Optional Arguments for SortSam", [
                ("SortSam For Merged BAM Create Index", False, True),
                ("SortSam For Merged BAM Create MD5 File", False, True)
            ])
        ])

        # SetNmMdAndUqTags (Picard) section
        self.create_picard_section(scrollable_frame, "SetNmMdAndUqTags (Picard)", [
            ("Required Arguments for SetNmMdAndUqTags", [
                ("SetnmmdUqtags Output", "SetnmmdUqtags.bam")
            ]),
            ("Optional Arguments for SetNmMdAndUqTags", [
                ("SetnmmdUqtags Is Bisulfite Sequence", False, True),
                ("Set Only Uq", False, True),
                ("SetnmmdUqtags Create Index", True, True),
                ("SetnmmdUqtags Create MD5 File", True, True)
            ])
        ])

        # MarkDuplicates (Picard) section
        self.create_picard_section(scrollable_frame, "MarkDuplicates (Picard)", [
            ("Required Arguments for MarkDuplicates", [
                ("Metrics File", "marked_dup_metrics.txt"),
                ("Mark Duplicates Output", "marked_duplicates.bam")
            ]),
            ("Optional Arguments for MarkDuplicates", [
                ("Assume Sort Order", "queryname", False, tk.StringVar),
                ("Barcode Tag", "", False, tk.StringVar),
                ("Clear Dt", True, True),
                ("MarkDuplicates Comment", "", False, tk.StringVar),
                ("Duplicate Scoring Strategy", "SUM_OF_BASE_QUALITIES", False, tk.StringVar),
                ("Max File Handles For Read Ends Map", 8000, False, tk.IntVar),
                ("Max Optical Duplicate Set Size", 300000, False, tk.IntVar),
                ("Max Sequences For Disk Read Ends Map", 50000, False, tk.IntVar),
                ("Optical Duplicate Pixel Distance", 2500, False, tk.IntVar),
                ("Program Group Command Line", "", False, tk.StringVar),
                ("Program Group Name", "", False, tk.StringVar),
                ("Program Group Version", "", False, tk.StringVar),
                ("Program Record Id", "", False, tk.StringVar),
                ("Read Name Regex", "", False, tk.StringVar),
                ("Read One Barcode Tag", "", False, tk.StringVar),
                ("Read Two Barcode Tag", "", False, tk.StringVar),
                ("Sorting Collection Size Ratio", 0.25, False, tk.DoubleVar),
                ("Tag Duplicate Set Members", False, True),
                ("Tagging Policy", "DontTag", False, tk.StringVar),
                ("MarkDuplicates Create Index", False, True),
                ("MarkDuplicates Create MD5 File", True, True)
            ])
        ])

        # SortSam (Picard) For Marked Duplicate BAM section
        self.create_picard_section(scrollable_frame, "SortSam (Picard) For Marked Duplicate BAM", [
            ("Required Arguments for SortSam", [
                ("SortSam For Marked Duplicate Output", "sorted_marked.bam", False, tk.StringVar),
                ("SortSam For Marked Duplicate Sort Order", "coordinate", False, tk.StringVar)
            ]),
            ("Optional Arguments for SortSam", [
                ("SortSam For Marked Duplicate Create Index", True, True),
                ("SortSam For Marked Duplicate Create MD5 File", False, True)
            ])
        ])

        # BaseRecalibrator section
        self.create_picard_section(scrollable_frame, "BaseRecalibrator (Picard)", [
            ("Required Arguments for BaseRecalibrator", [
                ("Base Recal Output", "recal_data.table", False, tk.StringVar),
                ("Base Recal Known Sites", "/media/STORAGE1/base_directory/1000G_phase1.snps.high_confidence.hg38.vcf.gz, /media/STORAGE1/base_directory/All_20180418.vcf.gz, /media/STORAGE1/base_directory/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz", False, tk.StringVar)
            ]),
            ("Optional Arguments for BaseRecalibrator", [
                ("BQSR BAQ Gap Open Penalty", 40.0, False, tk.DoubleVar),
                ("Cloud Index Prefetch Buffer", -1, False, tk.IntVar),
                ("Cloud Prefetch Buffer", 40, False, tk.IntVar),
                ("Default Base Qualities", -1, False, tk.IntVar),
                ("Deletions Default Quality", 45, False, tk.IntVar),
                ("Disable BAM Index Caching", False, True, tk.BooleanVar),
                ("GCS Max Retries", 20, False, tk.IntVar),
                ("Indels Context Size", 3, False, tk.IntVar),
                ("Insertions Default Quality", 45, False, tk.IntVar),
                ("Interval Merging Rule", "ALL", False, tk.StringVar),
                ("Intervals", "", False, tk.StringVar),
                ("Low Quality Tail", 2, False, tk.IntVar),
                ("Maximum Cycle Value", 500, False, tk.IntVar),
                ("Mismatches Context Size", 2, False, tk.IntVar),
                ("Mismatches Default Quality", -1, False, tk.IntVar),
                ("Preserve QScores Less Than", 6, False, tk.IntVar),
                ("Quantizing Levels", 16, False, tk.IntVar),
                ("Use Original Qualities", False, True, tk.BooleanVar)
            ])
        ])

        # ApplyBQSR (Picard) section
        self.create_picard_section(scrollable_frame, "ApplyBQSR (Picard)", [
            ("Required Arguments for ApplyBQSR", [
                ("ApplyBQSR Output", "bqsr_output.bam", False, tk.StringVar)
            ]),
            ("Optional Arguments for ApplyBQSR", [
                ("ApplyBQSR Cloud Index Prefetch Buffer", -1, False, tk.IntVar),
                ("ApplyBQSR Cloud Prefetch Buffer", 40, False, tk.IntVar),
                ("ApplyBQSR Disable BAM Index Caching", False, True),
                ("ApplyBQSR Emit Original Quals", False, True),
                ("ApplyBQSR GCS Max Retries", 20, False, tk.IntVar),
                ("ApplyBQSR Global Qscore Prior", -1.0, False, tk.DoubleVar),
                ("ApplyBQSR Interval Merging Rule", "ALL", False, tk.StringVar),
                ("ApplyBQSR Intervals", "", False, tk.StringVar),
                ("ApplyBQSR Preserve Qscores Less Than", 6, False, tk.IntVar),
                ("ApplyBQSR Quantize Quals", 0, False, tk.IntVar),
                ("ApplyBQSR Use Original Qualities", True, True)
            ])
        ])

        # HaplotypeCaller (Picard) section
        self.create_picard_section(scrollable_frame, "HaplotypeCaller (Picard)", [
            ("Required Arguments for HaplotypeCaller", [
                ("HaplotypeCaller Output", "haplotypecaller.vcf", False, tk.StringVar)
            ]),
            ("Optional Arguments for HaplotypeCaller", [
                ("HaplotypeCaller Activity Profile Out", "", False, tk.StringVar),
                ("HaplotypeCaller Alleles", "", False, tk.StringVar),
                ("HaplotypeCaller Annotate With Num Discovered Alleles", False, True),
                ("HaplotypeCaller Annotation", "", False, tk.StringVar),
                ("HaplotypeCaller Annotation Group", "", False, tk.StringVar),
                ("HaplotypeCaller Annotations To Exclude", "", False, tk.StringVar),
                ("HaplotypeCaller Assembly Region Out", "", False, tk.StringVar),
                ("HaplotypeCaller Base Quality Score Threshold", 18, False, tk.IntVar),
                ("HaplotypeCaller Cloud Index Prefetch Buffer", -1, False, tk.IntVar),
                ("HaplotypeCaller Cloud Prefetch Buffer", 40, False, tk.IntVar),
                ("HaplotypeCaller Contamination Fraction To Filter", 0.0, False, tk.DoubleVar),
                ("HaplotypeCaller Correct Overlapping Quality", False, True),
                ("HaplotypeCaller Dbsnp", "", False, tk.StringVar),
                ("HaplotypeCaller Disable Bam Index Caching", False, True),
                ("HaplotypeCaller Disable Sequence Dictionary Validation", False, True),
                ("HaplotypeCaller Founder Id", "", False, tk.StringVar),
                ("HaplotypeCaller Gcs Max Retries", 20, False, tk.IntVar),
                ("HaplotypeCaller Gcs Project For Requester Pays", "", False, tk.StringVar),
                ("HaplotypeCaller Graph Output", "", False, tk.StringVar),
                ("HaplotypeCaller Heterozygosity", 0.001, False, tk.DoubleVar),
                ("HaplotypeCaller Heterozygosity Stdev", 0.01, False, tk.DoubleVar),
                ("HaplotypeCaller Indel Heterozygosity", 1.25E-4, False, tk.DoubleVar),
                ("HaplotypeCaller Interval Merging Rule", "ALL", False, tk.StringVar),
                ("HaplotypeCaller Intervals", "", False, tk.StringVar),
                ("HaplotypeCaller Max Reads Per Alignment Start", 50, False, tk.IntVar),
                ("HaplotypeCaller Min Base Quality Score", 10, False, tk.IntVar),
                ("HaplotypeCaller Native Pair Hmm Threads", 4, False, tk.IntVar),
                ("HaplotypeCaller Native Pair Hmm Use Double Precision", False, True),
                ("HaplotypeCaller Num Reference Samples If No Call", 0, False, tk.IntVar),
                ("HaplotypeCaller Output Mode", "EMIT_VARIANTS_ONLY", False, tk.StringVar),
                ("HaplotypeCaller Pedigree", "", False, tk.StringVar),
                ("HaplotypeCaller Population Callset", "", False, tk.StringVar),
                ("HaplotypeCaller Sample Name", "", False, tk.StringVar),
                ("HaplotypeCaller Sample Ploidy", 2, False, tk.IntVar),
                ("HaplotypeCaller Sites Only Vcf Output", False, True),
                ("HaplotypeCaller Stand Call Conf", 30.0, False, tk.DoubleVar)
            ])
        ])

        # SelectVariants (GATK) SNVs Extraction section
        self.create_picard_section(scrollable_frame, "SelectVariants SNVs Extraction", [
            ("Required Arguments for SNV", [
                ("SNV Output", "snps.vcf", False, tk.StringVar)
            ]),
            ("Optional Arguments for SNV", [
                ("SNV Cloud Index Prefetch Buffer", -1, False, tk.IntVar),
                ("SNV Cloud Prefetch Buffer", 40, False, tk.IntVar),
                ("SNV Concordance", "", False, tk.StringVar),
                ("SNV Disable Bam Index Caching", False, True),
                ("SNV Discordance", "", False, tk.StringVar),
                ("SNV Exclude Filtered", False, True),
                ("SNV Exclude Ids", "", False, tk.StringVar),
                ("SNV Exclude Non Variants", False, True),
                ("SNV Exclude Sample Expressions", "", False, tk.StringVar),
                ("SNV Exclude Sample Name", "", False, tk.StringVar),
                ("SNV Gcs Max Retries", 20, False, tk.IntVar),
                ("SNV Interval Merging Rule", "ALL", False, tk.StringVar),
                ("SNV Intervals", "", False, tk.StringVar),
                ("SNV Invert Mendelian Violation", False, True),
                ("SNV Invert Select", False, True),
                ("SNV Keep Ids", "", False, tk.StringVar),
                ("SNV Keep Original Ac", False, True),
                ("SNV Keep Original Dp", False, True),
                ("SNV Max Filtered Genotypes", 2147483647, False, tk.IntVar),
                ("SNV Max Fraction Filtered Genotypes", 1.0, False, tk.DoubleVar),
                ("SNV Max Indel Size", 2147483647, False, tk.IntVar),
                ("SNV Max Nocall Fraction", 1.0, False, tk.DoubleVar),
                ("SNV Max Nocall Number", 2147483647, False, tk.IntVar),
                ("SNV Mendelian Violation", False, True),
                ("SNV Mendelian Violation Qual Threshold", 0.0, False, tk.DoubleVar),
                ("SNV Min Filtered Genotypes", 0, False, tk.IntVar),
                ("SNV Min Fraction Filtered Genotypes", 0.0, False, tk.DoubleVar),
                ("SNV Min Indel Size", 0, False, tk.IntVar),
                ("SNV Pedigree", "", False, tk.StringVar),
                ("SNV Preserve Alleles", False, True),
                ("SNV Remove Fraction Genotypes", 0.0, False, tk.DoubleVar),
                ("SNV Remove Unused Alternates", False, True),
                ("SNV Restrict Alleles To", "ALL", False, tk.StringVar),
                ("SNV Sample Expressions", "", False, tk.StringVar),
                ("SNV Sample Name", "", False, tk.StringVar),
                ("SNV Select Random Fraction", 0.0, False, tk.DoubleVar),
                ("SNV Select Type To Exclude", "", False, tk.StringVar),
                ("SNV Select Type To Include", "SNP", False, tk.StringVar),
                ("SNV Select Expressions", "", False, tk.StringVar),
                ("SNV Set Filtered Gt To Nocall", False, True)
            ])
        ])

        # SelectVariants (GATK) Indel Extraction section
        self.create_picard_section(scrollable_frame, "Indel (Picard)", [
            ("Required Arguments for Indel", [
                ("Indel Output", "indel.vcf", False, tk.StringVar)
            ]),
            ("Optional Arguments for Indel", [
                ("Indel Cloud Index Prefetch Buffer", -1, False, tk.IntVar),
                ("Indel Cloud Prefetch Buffer", 40, False, tk.IntVar),
                ("Indel Concordance", "", False, tk.StringVar),
                ("Indel Disable Bam Index Caching", False, True),
                ("Indel Discordance", "", False, tk.StringVar),
                ("Indel Exclude Filtered", False, True),
                ("Indel Exclude Ids", "", False, tk.StringVar),
                ("Indel Exclude Non Variants", False, True),
                ("Indel Exclude Sample Expressions", "", False, tk.StringVar),
                ("Indel Exclude Sample Name", "", False, tk.StringVar),
                ("Indel Gcs Max Retries", 20, False, tk.IntVar),
                ("Indel Interval Merging Rule", "ALL", False, tk.StringVar),
                ("Indel Intervals", "", False, tk.StringVar),
                ("Indel Invert Mendelian Violation", False, True),
                ("Indel Invert Select", False, True),
                ("Indel Keep Ids", "", False, tk.StringVar),
                ("Indel Keep Original Ac", False, True),
                ("Indel Keep Original Dp", False, True),
                ("Indel Max Filtered Genotypes", 2147483647, False, tk.IntVar),
                ("Indel Max Fraction Filtered Genotypes", 1.0, False, tk.DoubleVar),
                ("Indel Max Indel Size", 2147483647, False, tk.IntVar),
                ("Indel Max Nocall Fraction", 1.0, False, tk.DoubleVar),
                ("Indel Max Nocall Number", 2147483647, False, tk.IntVar),
                ("Indel Mendelian Violation", False, True),
                ("Indel Mendelian Violation Qual Threshold", 0.0, False, tk.DoubleVar),
                ("Indel Min Filtered Genotypes", 0, False, tk.IntVar),
                ("Indel Min Fraction Filtered Genotypes", 0.0, False, tk.DoubleVar),
                ("Indel Min Indel Size", 0, False, tk.IntVar),
                ("Indel Pedigree", "", False, tk.StringVar),
                ("Indel Preserve Alleles", False, True),
                ("Indel Remove Fraction Genotypes", 0.0, False, tk.DoubleVar),
                ("Indel Remove Unused Alternates", False, True),
                ("Indel Restrict Alleles To", "ALL", False, tk.StringVar),
                ("Indel Sample Expressions", "", False, tk.StringVar),
                ("Indel Sample Name", "", False, tk.StringVar),
                ("Indel Select Random Fraction", 0.0, False, tk.DoubleVar),
                ("Indel Select Type To Exclude", "", False, tk.StringVar),
                ("Indel Select Type To Include", "INDEL", False, tk.StringVar),
                ("Indel Select Expressions", "", False, tk.StringVar),
                ("Indel Set Filtered Gt To Nocall", False, True)
            ])
        ])

        # CollectVariantCallingMetrics (Picard) section
        self.create_picard_section(scrollable_frame, "CollectVariantCallingMetrics (Picard)", [
            ("Required Arguments for CollectVariantCallingMetrics", [
                ("Collectvariantcallingmetrics Output", "variant_metrics", False, tk.StringVar),
                ("Collectvariantcallingmetrics Dbsnp", "/media/STORAGE1/base_directory/All_20180418.vcf.gz", False, tk.StringVar)
            ]),
            ("Optional Arguments for CollectVariantCallingMetrics", [
                ("Collectvariantcallingmetrics Gvcf Input", False, True),
                ("Collectvariantcallingmetrics Sequence Dictionary", "", False, tk.StringVar),
                ("Collectvariantcallingmetrics Target Intervals", "", False, tk.StringVar),
                ("Collectvariantcallingmetrics Thread Count", 10, False, tk.IntVar)
            ])
        ])

        # Convert2ANNOVAR Section in the configuration page
        self.create_convert2annovar_section(scrollable_frame)

        # ANNOVAR Section
        self.create_picard_section(scrollable_frame, "ANNOVAR", [
            ("Required Arguments for ANNOVAR", [
                ("Database Directory", "/media/STORAGE1/programs/annovar/humandb", False, tk.StringVar),
                ("Genome Build Version", "hg38", False, tk.StringVar),
                ("Output Prefix", "myanno", False, tk.StringVar)
            ]),
            ("Optional Arguments for ANNOVAR", [
                ("Protocol", "refGene,cytoBand,exac03,avsnp147,dbnsfp30a", False, tk.StringVar),
                ("Operation", "g,r,f,f,f", False, tk.StringVar),
                ("Nastring", ".", False, tk.StringVar),
                ("No Polish", False, True),
                ("Remove Intermediate Files", False, True)
            ])
        ])

    # Browsing function for ANNOVAR Database Directory
    def browse_annovar_db_directory(self, entry):
        # Open a directory selection dialog
        directory = filedialog.askdirectory(title="Select ANNOVAR Database Directory")
        if directory:
            # If a directory was selected, update the entry widget with the directory path
            entry.delete(0, tk.END)
            entry.insert(0, directory)

    def create_picard_section(self, parent, title, arguments):
        # Create a LabelFrame with the given title for grouping related widgets
        frame = tk.LabelFrame(parent, text=title, font=("Arial", 16, "bold"))
        frame.pack(fill="x", padx=50, pady=10)  # Pack with padding and fill in x direction

        row = 0
        entry_width = 30

        # Iterate over each section title and its corresponding arguments
        for section_title, section_arguments in arguments:
            # Add a bold label for each section title
            tk.Label(frame, text=section_title, font=("Arial", 12, "bold")).grid(row=row, column=0, sticky="w", padx=10,
                                                                                 pady=5, columnspan=6)
            row += 1
            for i, arg in enumerate(section_arguments):
                label_text = arg[0]
                default = arg[1]
                is_checkbox = len(arg) > 2 and arg[2]
                entry_type = arg[3] if len(arg) > 3 else tk.StringVar
                col = (i % 3) * 2
                label_key = label_text.lower().replace(" ", "_")

                # Check for specific fields that require a file browser
                if label_key in {"collectvariantcallingmetrics_dbsnp", "base_recal_known_sites"}:
                    self.create_browse_file_section(frame, label_text, row, column=col // 2, default=default,
                                                    entry_type=entry_type, width=entry_width)
                elif is_checkbox:
                    self.create_checkbox_section(frame, label_text, row, column=col // 2, default=default)
                else:
                    self.create_entry_section(frame, label_text, row, column=col // 2, default=default,
                                              entry_type=entry_type, width=entry_width)
                if (i + 1) % 3 == 0:
                    row += 1
            row += 1

    def create_entry_section(self, parent, label_text, row, column=0, default="", entry_type=tk.StringVar, width=30):
        # Create a label for the entry
        label = tk.Label(parent, text=label_text, width=35, anchor="w")
        label.grid(row=row, column=column * 3, padx=5, pady=5, sticky="w")

        # Create a variable to hold the entry's value
        entry_var = entry_type()
        entry_var.set(default)  # Set the default value
        entry = tk.Entry(parent, textvariable=entry_var, width=width)
        entry.grid(row=row, column=column * 3 + 1, padx=5, pady=5, sticky="we")

        # Store the entry variable in the instance's attributes
        setattr(self, label_text.lower().replace(" ", "_"), entry_var)

        return entry_var

    def create_checkbox_section(self, parent, label_text, row, column=0, default=False):
        # Create a label for the checkbox
        label = tk.Label(parent, text=label_text, width=35, anchor="w")
        label.grid(row=row, column=column * 3, padx=5, pady=5, sticky="w")

        # Create a BooleanVar to hold the checkbox's value
        var = tk.BooleanVar(value=default)
        checkbox = tk.Checkbutton(parent, variable=var)
        checkbox.grid(row=row, column=column * 3 + 1, padx=5, pady=5, sticky="w")

        # Store the checkbox variable in the instance's attributes
        setattr(self, label_text.lower().replace(" ", "_"), var)

        return var

    def create_browse_file_section(self, parent, label_text, row, column=0, default="", entry_type=tk.StringVar,
                                   width=30):
        # Create a label for the entry and browse button
        label = tk.Label(parent, text=label_text, width=35, anchor="w")
        label.grid(row=row, column=column * 3, padx=5, pady=5, sticky="w")

        # Create a variable to hold the entry's value
        entry_var = entry_type()
        entry_var.set(default)  # Set the default value
        entry = tk.Entry(parent, textvariable=entry_var, width=width)
        entry.grid(row=row, column=column * 3 + 1, padx=5, pady=5, sticky="we")

        # Create a browse button that opens a file dialog for multiple file selections
        button = ttk.Button(parent, text="Browse", command=lambda: self.browse_files(entry))
        button.grid(row=row, column=column * 3 + 2, padx=5, pady=5)

        # Store the entry variable in the instance's attributes
        setattr(self, label_text.lower().replace(" ", "_"), entry_var)

        return entry_var

    def browse_files(self, entry):
        # Open a file selection dialog for multiple files
        file_paths = filedialog.askopenfilenames(title="Select files",
                                                 filetypes=[("VCF files", "*.vcf"), ("All files", "*.*")])
        if file_paths:
            # If files were selected, update the entry widget with the joined file paths
            entry.delete(0, tk.END)
            entry.insert(0, ', '.join(file_paths))  # Join paths with a semicolon and space

    def create_convert2annovar_section(self, parent):
        frame = tk.LabelFrame(parent, text="Convert2ANNOVAR", font=("Arial", 16, "bold"))
        frame.pack(fill="x", padx=50, pady=10)

        # Optional Arguments for Convert2ANNOVAR
        self.convert2annovar_with_zyg = tk.BooleanVar()
        self.withzyg_checkbox = tk.Checkbutton(frame, text="With Zyg", variable=self.convert2annovar_with_zyg,
                                               command=self.on_withzyg_or_withfreq_change)
        self.withzyg_checkbox.pack(side="left", padx=5, pady=5, anchor="w")

        self.convert2annovar_with_freq = tk.BooleanVar()
        self.withfreq_checkbox = tk.Checkbutton(frame, text="With Freq", variable=self.convert2annovar_with_freq,
                                                command=self.on_withzyg_or_withfreq_change)
        self.withfreq_checkbox.pack(side="left", padx=5, pady=5, anchor="w")

        self.convert2annovar_include_info = tk.BooleanVar()
        tk.Checkbutton(frame, text="Include Info", variable=self.convert2annovar_include_info).pack(side="left", padx=5,
                                                                                                    pady=5, anchor="w")

        self.convert2annovar_comment = tk.BooleanVar()
        tk.Checkbutton(frame, text="Comment", variable=self.convert2annovar_comment).pack(side="left", padx=5, pady=5,
                                                                                          anchor="w")

    def on_withzyg_or_withfreq_change(self):
        """Callback function to ensure that 'With Zyg' and 'With Freq' are mutually exclusive."""
        if self.convert2annovar_with_zyg.get():
            self.withfreq_checkbox.config(state=tk.DISABLED)
        else:
            self.withfreq_checkbox.config(state=tk.NORMAL)

        if self.convert2annovar_with_freq.get():
            self.withzyg_checkbox.config(state=tk.DISABLED)
        else:
            self.withzyg_checkbox.config(state=tk.NORMAL)

