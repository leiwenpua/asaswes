import os
import json
import threading
import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext
import process_management


# Define the main application class, inheriting from Tkinter's Tk class
class SequenceAssemblyApp(tk.Tk):
    def __init__(self):
        super().__init__()

        # Set the pages, window title and size
        self.pages = {}
        self.title("Automated Sequence Assembly Software")
        self.geometry("1200x800")

        # Initialize the GUI components
        self.create_widgets()

    # Method to create the GUI components (widgets)
    def create_widgets(self):
        # Create the menu bar at the top of the window
        menu_bar = tk.Menu(self)
        self.config(menu=menu_bar)

        # Add menu commands to switch between different pages
        menu_bar.add_command(label="Main", command=lambda: self.show_page("MainPage"))
        menu_bar.add_command(label="Configuration", command=lambda: self.show_page("ConfigurationPage"))
        menu_bar.add_command(label="Program", command=lambda: self.show_page("ProgramPage"))

        # Create a Help menu with a dropdown for different options
        help_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="Help", menu=help_menu)

        # Add dropdown options to the Help menu
        help_menu.add_command(label="About this system", command=lambda: self.show_page("AboutPage"))
        help_menu.add_command(label="Installation", command=lambda: self.show_page("InstallationPage"))
        help_menu.add_command(label="Getting Started", command=lambda: self.show_page("GettingStartedPage"))
        help_menu.add_command(label="Pipeline Workflow", command=lambda: self.show_page("WorkflowPage"))
        help_menu.add_command(label="Common Issues", command=lambda: self.show_page("IssuesPage"))

        # Container to hold all the pages (frames)
        container = tk.Frame(self)
        container.pack(fill="both", expand=True)

        # Loop through the pages and initialize each one
        for page in (MainPage, ConfigurationPage, ProgramPage, AboutPage, InstallationPage, GettingStartedPage, WorkflowPage, IssuesPage):

            page_name = page.__name__  # Get the class name as a string
            frame = page(parent=container, controller=self)  # Create an instance of the page
            self.pages[page_name] = frame  # Store the frame in the dictionary
            frame.grid(row=0, column=0, sticky="nsew")  # Place the frame in the grid

        # Configure the grid to expand with the window size
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # Show the initial page (MainPage)
        self.show_page("MainPage")

    # Method to bring a specific page to the front
    def show_page(self, page_name):
        frame = self.pages[page_name]
        frame.tkraise()


class MainPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        # Label indicating the Main Page
        label = tk.Label(self, text="Main Page", font=("Arial", 24))
        label.pack(pady=10)


        # Create and pack a section for selecting the forward read sequence file
        self.forward_read = self.create_browse_section("Forward Read Sequence", self.browse_file)
        self.forward_read.pack(fill="x", padx=50, pady=5)

        # Create and pack a section for selecting the reverse read sequence file
        self.reverse_read = self.create_browse_section("Reverse Read Sequence", self.browse_file)
        self.reverse_read.pack(fill="x", padx=50, pady=5)

        # Create and pack a section for selecting the reference sequence file
        self.reference_sequence = self.create_browse_section("Reference Genome", self.browse_reference_file)
        self.reference_sequence.pack(fill="x", padx=50, pady=5)

        # Create and pack a section for directory selection
        self.output_directory = self.create_browse_section("Directory", self.browse_directory)
        self.output_directory.pack(fill="x", padx=50, pady=5)

        # Create a frame for Save and Load buttons
        button_frame = ttk.Frame(self)
        button_frame.pack(pady=20)

        # Save button to save the current configuration as default
        save_button = ttk.Button(button_frame, text="Save as Default", command=self.save_configuration)
        save_button.pack(side="left", padx=10)

        # Load button to load the default configuration
        load_button = ttk.Button(button_frame, text="Load Default", command=self.load_configuration)
        load_button.pack(side="left", padx=10)

        # Button to start the analysis process
        run_button = ttk.Button(self, text="Run", command=self.start_analysis_thread)
        run_button.pack(pady=20)

        # ScrolledText widget to display logs or messages
        self.log_text = scrolledtext.ScrolledText(self, wrap=tk.WORD, height=10)
        self.log_text.pack(fill="both", expand=True, padx=50, pady=5)

    def save_configuration(self):
        # Gather the configuration data from the entry widgets
        config = {
            'forward_read': self.forward_read.entry.get(),
            'reverse_read': self.reverse_read.entry.get(),
            'reference_sequence': self.reference_sequence.entry.get(),
            'output_directory': self.output_directory.entry.get(),
        }

        # Define the configuration directory and ensure it exists
        config_directory = os.path.join(os.path.dirname(__file__), 'configuration')
        os.makedirs(config_directory, exist_ok=True)

        # Write the configuration data to a JSON file
        with open(os.path.join(config_directory, 'basic_input.json'), 'w') as f:
            json.dump(config, f, indent=4)

        # Inform the user that the configuration has been saved
        tk.messagebox.showinfo("Save Configuration", "Configuration saved successfully!")

    def load_configuration(self):
        # Define the path to the configuration file
        config_path = os.path.join(os.path.dirname(__file__), 'configuration', 'basic_input.json')

        # Check if the configuration file exists
        if os.path.exists(config_path):
            # Read the configuration file and load the data
            with open(config_path, 'r') as f:
                config = json.load(f)

            # Update the entry widgets with the loaded configuration data
            self.forward_read.entry.delete(0, tk.END)
            self.forward_read.entry.insert(0, config['forward_read'])
            self.reverse_read.entry.delete(0, tk.END)
            self.reverse_read.entry.insert(0, config['reverse_read'])
            self.reference_sequence.entry.delete(0, tk.END)
            self.reference_sequence.entry.insert(0, config['reference_sequence'])
            self.output_directory.entry.delete(0, tk.END)
            self.output_directory.entry.insert(0, config['output_directory'])

            # Inform the user that the configuration has been Loaded
            tk.messagebox.showinfo("Load Configuration", "Configuration loaded successfully!")
        else:
            # Show an error message if the configuration file does not exist
            tk.messagebox.showerror("Load Error", "No configuration file found.")

    # Function to start the analysis in a separate thread (placeholder)
    def start_analysis_thread(self):
        # Clear the log area
        self.log_text.delete('1.0', tk.END)

        # Check if all required fields are filled
        if not self.validate_inputs():
            return  # Exit the function if inputs are not valid

        # Check if program locations are filled and valid
        if not self.validate_program_locations():
            return  # Exit the function if program locations are not valid

        # Start the analysis in a new thread
        analysis_thread = threading.Thread(target=self.start_analysis)
        analysis_thread.start()

    def validate_inputs(self):
        # Check if forward read sequence file entry is filled
        if not self.forward_read.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please select the forward read sequence file.")
            return False

        # Check if reverse read sequence file entry is filled
        if not self.reverse_read.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please select the reverse read sequence file.")
            return False

        # Check if reference sequence file entry is filled
        if not self.reference_sequence.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please select the reference genome file.")
            return False

        # Check if output directory entry is filled
        if not self.output_directory.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please select an output directory.")
            return False

        # If all fields are filled, return True
        return True

    def validate_program_locations(self):
        # Access the ProgramPage instance from the controller
        program_page = self.controller.pages["ProgramPage"]

        # Check if FASTQC location is filled
        if not program_page.fastqc_location.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please specify the location of the FASTQC program.")
            return False

        # Check if GATK location is filled
        if not program_page.gatk_location.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please specify the location of the GATK program.")
            return False

        # Check if BWA location is filled
        if not program_page.bwa_location.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please specify the location of the BWA program.")
            return False

        # Check if ANNOVAR location is filled
        if not program_page.annovar_location.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please specify the location of the ANNOVAR program.")
            return False

        # If all program locations are filled, return True
        return True

    def start_analysis(self):
        try:
            # Log the start of the analysis
            self.log_text.insert(tk.END, "Starting analysis...\n")
            self.update_idletasks()  # Update the GUI to reflect changes immediately

            # Retrieve the output directory from the corresponding entry widget
            output_directory = self.output_directory.entry.get()

            # Retrieve the forward read sequence file path from the corresponding entry widget
            forward_read_path = self.forward_read.entry.get()

            # Retrieve the reverse read sequence file path from the corresponding entry widget
            reverse_read_path = self.reverse_read.entry.get()

            # Retrieve the reference sequence file path from the corresponding entry widget
            reference_sequence_path = self.reference_sequence.entry.get()

            # Access the ProgramPage and ConfigurationPage instances from the controller
            program_page = self.controller.pages["ProgramPage"]
            config_page = self.controller.pages["ConfigurationPage"]

            # Retrieve the paths for the FastQC, GATK, and BWA programs from the ProgramPage
            fastqc_path = program_page.fastqc_location.entry.get()
            gatk_path = program_page.gatk_location.entry.get()
            bwa_path = program_page.bwa_location.entry.get()

            # Define a helper function to filter out optional arguments (those with None or empty values)
            def filter_optional_args(args):
                return {k: v for k, v in args.items() if v}

            # Define a function to log messages to the log_text widget and keep it scrolled to the bottom
            def log_function(log):
                self.log_text.insert(tk.END, log)
                self.log_text.yview(tk.END)  # Scroll to the end of the log
                self.update_idletasks()  # Update the GUI

            # Verify if the forward read file exists in the specified output directory
            if not self.is_file_in_directory(forward_read_path, output_directory):
                self.log_text.insert(tk.END,
                                     f"Error: Forward read file {forward_read_path} is not in the output directory {output_directory}.\n")
                return  # Exit the function if the file is not found

            # Verify if the reverse read file exists in the specified output directory
            if not self.is_file_in_directory(reverse_read_path, output_directory):
                self.log_text.insert(tk.END,
                                     f"Error: Reverse read file {reverse_read_path} is not in the output directory {output_directory}.\n")
                return  # Exit the function if the file is not found

            # Verify if the reference sequence file exists in the specified output directory
            if not self.is_file_in_directory(reference_sequence_path, output_directory):
                self.log_text.insert(tk.END,
                                     f"Error: Reference sequence file {reference_sequence_path} is not in the output directory {output_directory}.\n")
                return  # Exit the function if the file is not found

            # FastQC Analysis
            try:
                # Log the start of the FastQC process
                self.log_text.insert(tk.END, "Running FastQC...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Run FastQC on the forward and reverse read files
                fastqc_result = process_management.run_fastqc(
                    fastqc_path,  # Path to the FastQC executable
                    [forward_read_path, reverse_read_path],  # List of input sequence files
                    output_directory,  # Directory where FastQC output will be stored
                    log_function  # Function to log messages during the process
                )

                # Log the result of the FastQC process
                self.log_text.insert(tk.END, fastqc_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if FastQC fails
                self.log_text.insert(tk.END, f"FastQC error: {str(e)}\n")

            # FastqToSam Conversion
            try:
                # Log the start of the FastqToSam process
                self.log_text.insert(tk.END, "Running FastqToSam...\n")
                self.update_idletasks()  # Update the GUI to reflect changes immediately

                # Retrieve required configuration settings from the ConfigurationPage
                sample_name = config_page.sample_name.get()  # The sample name
                platform = config_page.platform.get()  # The sequencing platform used
                output_bam = os.path.join(output_directory,
                                          config_page.fastq_to_sam_output.get())  # Output BAM file path

                # Collect optional arguments from the ConfigurationPage, filtering out any None or empty values
                optional_args = filter_optional_args({
                    "READ_GROUP_NAME": config_page.read_group_name.get(),
                    "LIBRARY_NAME": config_page.library_name.get(),
                    "DESCRIPTION": config_page.description.get(),
                    "SORT_ORDER": config_page.fastqtosam_sort_order.get(),
                    "MAX_Q": config_page.max_q.get(),
                    "MIN_Q": config_page.min_q.get(),
                    "USE_SEQUENTIAL_FASTQS": str(config_page.use_sequential_fastqs.get()),
                    "ALLOW_AND_IGNORE_EMPTY_LINES": str(config_page.allow_and_ignore_empty_lines.get()),
                    "PLATFORM_MODEL": config_page.platform_model.get(),
                    "PLATFORM_UNIT": config_page.platform_unit.get(),
                    "PREDICTED_INSERT_SIZE": config_page.predicted_insert_size.get(),
                    "PROGRAM_GROUP": config_page.program_group.get(),
                    "QUALITY_FORMAT": config_page.quality_format.get(),
                    "RUN_DATE": config_page.run_date.get(),
                    "SEQUENCING_CENTER": config_page.sequencing_center.get()
                })

                # Run the FastqToSam process with the provided inputs and optional arguments
                fastqtosam_result = process_management.run_fastqtosam(
                    gatk_path,  # Path to the GATK executable
                    forward_read_path,  # Path to the forward read sequence file
                    reverse_read_path,  # Path to the reverse read sequence file
                    output_bam,  # Path to the output BAM file
                    sample_name,  # Sample name
                    platform,  # Sequencing platform
                    optional_args,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the FastqToSam process
                self.log_text.insert(tk.END, fastqtosam_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if FastqToSam fails
                self.log_text.insert(tk.END, f"FastqToSam error: {str(e)}\n")

            # Samtools View to capture header
            try:
                # Log the start of the Samtools View process
                self.log_text.insert(tk.END, "Running Samtools View to capture header...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Run the Samtools View command
                samtools_result = process_management.run_samtools_view(
                    output_bam,  # Path to the BAM file generated by FastqToSam
                    output_directory,  # Directory where any output or logs will be stored
                    log_function  # Function to log messages during the process
                )

                # Log the result of the Samtools View process
                self.log_text.insert(tk.END, samtools_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if Samtools View fails
                self.log_text.insert(tk.END, f"Samtools View error: {str(e)}\n")

            # MarkIlluminaAdapters Process
            try:
                # Log the start of the MarkIlluminaAdapters process
                self.log_text.insert(tk.END, "Running MarkIlluminaAdapters...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Construct the path for the metrics file
                metrics_file = os.path.join(output_directory, config_page.mark_illumina_metrics.get())

                # Construct the path for the output BAM file after marking adapters
                mark_adapters_output_bam = os.path.join(output_directory, config_page.mark_illumina_output.get())

                # Collect optional arguments from the ConfigurationPage
                optional_args = {
                    "ADAPTER_TRUNCATION_LENGTH": config_page.adapter_truncation_length.get(),
                    "ADAPTERS": config_page.adapters.get(),
                    "FIVE_PRIME_ADAPTER": config_page.five_prime_adapter.get(),
                    "MAX_ERROR_RATE_PE": config_page.max_error_rate_pe.get(),
                    "MAX_ERROR_RATE_SE": config_page.max_error_rate_se.get(),
                    "MIN_MATCH_BASES_PE": config_page.min_match_bases_pe.get(),
                    "MIN_MATCH_BASES_SE": config_page.min_match_bases_se.get(),
                    "NUM_ADAPTERS_TO_KEEP": config_page.num_adapters_to_keep.get(),
                    "THREE_PRIME_ADAPTER": config_page.three_prime_adapter.get()
                }

                # Run the MarkIlluminaAdapters process with the provided inputs and optional arguments
                markilluminaadapters_result = process_management.run_markilluminaadapters(
                    gatk_path,  # Path to the GATK executable
                    output_bam,  # Path to the BAM file generated earlier
                    metrics_file,  # Path to the file where metrics will be saved
                    mark_adapters_output_bam,  # Path to the output BAM file after marking adapters
                    optional_args,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the MarkIlluminaAdapters process
                self.log_text.insert(tk.END, markilluminaadapters_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if MarkIlluminaAdapters fails
                self.log_text.insert(tk.END, f"MarkIlluminaAdapters error: {str(e)}\n")

            # SamToFastq Process
            try:
                # Log the start of the SamToFastq process
                self.log_text.insert(tk.END, "Running SamToFastq...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Construct the path for the output FASTQ file
                fastq_output = os.path.join(output_directory, config_page.sam_to_fastq_output.get())

                # Collect optional arguments from the ConfigurationPage
                optional_args = {
                    "CLIPPING_ACTION": config_page.clipping_action.get(),
                    "CLIPPING_ATTRIBUTE": config_page.clipping_attribute.get(),
                    "CLIPPING_MIN_LENGTH": config_page.clipping_min_length.get(),
                    "INCLUDE_NON_PF_READS": config_page.include_non_pf_reads.get(),
                    "INCLUDE_NON_PRIMARY_ALIGNMENTS": config_page.include_non_primary_alignments.get(),
                    "INTERLEAVE": config_page.interleave.get(),
                    "QUALITY": config_page.quality.get(),
                    "READ1_MAX_BASES_TO_WRITE": config_page.read1_max_bases_to_write.get(),
                    "READ1_TRIM": config_page.read1_trim.get(),
                    "READ2_MAX_BASES_TO_WRITE": config_page.read2_max_bases_to_write.get(),
                    "READ2_TRIM": config_page.read2_trim.get(),
                    "RG_TAG": config_page.rg_tag.get()
                }

                # Run the SamToFastq process with the provided inputs and optional arguments
                samtofastq_result = process_management.run_samtofastq(
                    gatk_path,  # Path to the GATK executable
                    mark_adapters_output_bam,  # Path to the BAM file with marked adapters
                    fastq_output,  # Path to the output FASTQ file
                    optional_args,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the SamToFastq process
                self.log_text.insert(tk.END, samtofastq_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if SamToFastq fails
                self.log_text.insert(tk.END, f"SamToFastq error: {str(e)}\n")

            # BWA Index Process
            try:
                # Log the start of the BWA Index process
                self.log_text.insert(tk.END, "Running BWA Index...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Retrieve the index algorithm type from the ConfigurationPage
                index_algorithm = config_page.index_algorithm.get()

                # Define a log function to update the log_text widget with real-time logs
                def log_function(log):
                    self.log_text.insert(tk.END, log)
                    self.log_text.yview(tk.END)  # Scroll to the end of the log
                    self.update_idletasks()  # Update the GUI

                # Run the BWA Index command with the specified parameters
                bwa_index_result = process_management.run_bwa_index(
                    bwa_path,  # Path to the BWA executable
                    reference_sequence_path,  # Path to the reference sequence file
                    index_algorithm,  # Algorithm used for indexing
                    log_function  # Function to log messages during the process
                )

                # Log the result of the BWA Index process
                self.log_text.insert(tk.END, bwa_index_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if BWA Index fails
                self.log_text.insert(tk.END, f"BWA Index error: {str(e)}\n")

            # Update the GUI after the BWA Index process
            self.update_idletasks()


            # BWA Mem Process
            try:
                # Log the start of the BWA Mem process
                self.log_text.insert(tk.END, "Running BWA Mem...\n\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define a log function for real-time updates
                def log_function(log):
                    self.log_text.insert(tk.END, log)
                    self.log_text.yview(tk.END)  # Scroll to the end of the log
                    self.update_idletasks()  # Update the GUI

                # Define the output SAM file path
                sam_output = os.path.join(output_directory, "bwa.sam")

                # Retrieve the number of threads to use from the ConfigurationPage
                threads = config_page.threads.get()

                # Run the BWA Mem command with the provided inputs and optional arguments
                bwa_mem_result = process_management.run_bwa_mem(
                    bwa_path,  # Path to the BWA executable
                    reference_sequence_path,  # Path to the reference sequence file
                    fastq_output,  # Path to the input FASTQ file(s)
                    sam_output,  # Path to the output SAM file
                    threads,  # Number of threads to use for parallel processing
                    log_function  # Function to log messages during the process
                )

                # Log the result of the BWA Mem process
                self.log_text.insert(tk.END, bwa_mem_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if BWA Mem fails
                self.log_text.insert(tk.END, f"BWA Mem error: {str(e)}\n")

            # Update the GUI after the BWA Mem process
            self.update_idletasks()

            # Create Sequence Dictionary Process
            try:
                # Log the start of the Create Sequence Dictionary process
                self.log_text.insert(tk.END, "Creating Sequence Dictionary...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the path for the output sequence dictionary file
                reference_dict = os.path.join(output_directory, os.path.basename(reference_sequence_path).replace('.fa', '.dict'))

                # Collect optional arguments from the ConfigurationPage
                optional_args = filter_optional_args({
                    "GENOME_ASSEMBLY": config_page.genome_assembly.get(),
                    "NUM_SEQUENCES": config_page.number_of_sequences.get(),
                    "SPECIES": config_page.species.get(),
                    "TRUNCATE_NAMES_AT_WHITESPACE": str(config_page.truncate_names_at_whitespace.get()),
                    "URI": config_page.uri.get()
                })

                # Run the Create Sequence Dictionary command with the provided inputs and optional arguments
                sequence_dict_result = process_management.run_create_sequence_dictionary(
                    gatk_path,  # Path to the GATK executable
                    reference_sequence_path,  # Path to the reference sequence file (.fa)
                    reference_dict,  # Path to the output dictionary file (.dict)
                    optional_args,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the Create Sequence Dictionary process
                self.log_text.insert(tk.END, sequence_dict_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the Create Sequence Dictionary process fails
                self.log_text.insert(tk.END, f"Create Sequence Dictionary error: {str(e)}\n")

            # Update the GUI after the Create Sequence Dictionary process
            self.update_idletasks()

            # Create FAI Index Process using samtools faidx
            try:
                # Log the start of the FAI index creation process
                self.log_text.insert(tk.END, "Creating FAI index...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Run the samtools faidx command to create the FAI index
                fai_index_result = process_management.run_samtools_faidx(
                    reference_sequence_path,  # Path to the reference sequence file (FASTA format)
                    log_function  # Function to log messages during the process
                )

                # Log the result of the FAI index creation process
                self.log_text.insert(tk.END, fai_index_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the FAI index creation fails
                self.log_text.insert(tk.END, f"Samtools faidx error: {str(e)}\n")

            # Update the GUI after the FAI index creation process
            self.update_idletasks()

            # MergeBamAlignment Process
            try:
                # Log the start of the MergeBamAlignment process
                self.log_text.insert(tk.END, "Running MergeBamAlignment...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Collect optional arguments from the ConfigurationPage
                optional_args = filter_optional_args({
                    "ADD_MATE_CIGAR": config_page.add_mate_cigar.get(),
                    "ALIGNED_READS_ONLY": config_page.aligned_reads_only.get(),
                    "ALIGNER_PROPER_PAIR_FLAGS": config_page.aligner_proper_pair_flags.get(),
                    "ATTRIBUTES_TO_REMOVE": config_page.attributes_to_remove.get(),
                    "ATTRIBUTES_TO_RETAIN": config_page.attributes_to_retain.get(),
                    "ATTRIBUTES_TO_REVERSE": config_page.attributes_to_reverse.get(),
                    "ATTRIBUTES_TO_REVERSE_COMPLEMENT": config_page.attributes_to_reverse_complement.get(),
                    "CLIP_ADAPTERS": config_page.clip_adapters.get(),
                    "CLIP_OVERLAPPING_READS": config_page.clip_overlapping_reads.get(),
                    "EXPECTED_ORIENTATIONS": config_page.expected_orientations.get(),
                    "INCLUDE_SECONDARY_ALIGNMENTS": config_page.include_secondary_alignments.get(),
                    "IS_BISULFITE_SEQUENCE": config_page.mergebamalignment_is_bisulfite_sequence.get(),
                    "MATCHING_DICTIONARY_TAGS": config_page.matching_dictionary_tags.get(),
                    "MAX_INSERTIONS_OR_DELETIONS": config_page.max_insertions_or_deletions.get(),
                    "MIN_UNCLIPPED_BASES": config_page.min_unclipped_bases.get(),
                    "PRIMARY_ALIGNMENT_STRATEGY": config_page.primary_alignment_strategy.get(),
                    "PROGRAM_GROUP_COMMAND_LINE": config_page.program_group_command_line.get(),
                    "PROGRAM_GROUP_NAME": config_page.program_group_name.get(),
                    "PROGRAM_GROUP_VERSION": config_page.program_group_version.get(),
                    "PROGRAM_RECORD_ID": config_page.program_record_id.get(),
                    "READ1_ALIGNED_BAM": config_page.read1_aligned_bam.get(),
                    "READ1_TRIM": config_page.read1_trim.get(),
                    "READ2_ALIGNED_BAM": config_page.read2_aligned_bam.get(),
                    "READ2_TRIM": config_page.read2_trim.get(),
                    "SORT_ORDER": config_page.mergebamalignment_sort_order.get(),
                    "UNMAP_CONTAMINANT_READS": config_page.unmap_contaminant_reads.get(),
                    "UNMAPPED_READ_STRATEGY": config_page.unmapped_read_strategy.get()
                })

                # Only add PROGRAM_GROUP arguments if all are set
                if not all(optional_args.get(arg) for arg in
                           ['PROGRAM_RECORD_ID', 'PROGRAM_GROUP_NAME', 'PROGRAM_GROUP_VERSION',
                            'PROGRAM_GROUP_COMMAND_LINE']):
                    optional_args = {key: value for key, value in optional_args.items() if
                                     key not in ['PROGRAM_RECORD_ID', 'PROGRAM_GROUP_NAME', 'PROGRAM_GROUP_VERSION',
                                                 'PROGRAM_GROUP_COMMAND_LINE']}

                # Define the path for the merged BAM output file
                merge_bam_output = os.path.join(output_directory, config_page.merge_align_bam_output.get())

                # Run the MergeBamAlignment command with the provided inputs and optional arguments
                mergebamalignment_result = process_management.run_mergebamalignment(
                    gatk_path,  # Path to the GATK executable
                    output_bam,  # Unaligned BAM (original reads)
                    sam_output,  # Aligned BAM (output from BWA Mem)
                    reference_sequence_path,  # Path to the reference sequence file
                    merge_bam_output,  # Path to the merged BAM output file
                    optional_args,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the MergeBamAlignment process
                self.log_text.insert(tk.END, mergebamalignment_result + "\n")
            except ValueError as ve:
                # Handle specific exceptions like ValueError and log the error
                self.log_text.insert(tk.END, f"MergeBamAlignment error: {str(ve)}\n")
            except Exception as e:
                # Handle general exceptions and log the error
                self.log_text.insert(tk.END, f"MergeBamAlignment error: {str(e)}\n")


            # Extract unmapped reads
            try:
                # Log the start of the extraction process for unmapped reads
                self.log_text.insert(tk.END, "Extracting unmapped reads...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the output path for the unmapped reads
                unmapped_output = os.path.join(output_directory, "mergeBA_unmapped.txt")

                # Run the samtools command to extract unmapped reads
                extract_unmapped_result = process_management.run_samtools_view_unmapped(
                    merge_bam_output,  # Path to the merged BAM file
                    unmapped_output,  # Output path for the unmapped reads
                    log_function
                )

                # Log the result of the extraction process
                self.log_text.insert(tk.END, extract_unmapped_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the extraction fails
                self.log_text.insert(tk.END, f"Samtools view error: {str(e)}\n")

            # SortSam For Merged BAM
            try:
                # Log the start of the SortSam process
                self.log_text.insert(tk.END, "Running SortSam...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Retrieve the sort order and output path from the configuration page
                sort_order = config_page.sortsam_for_merged_bam_sort_order.get()
                sorted_output = os.path.join(output_directory, config_page.sortsam_for_merged_bam_output.get())

                # Collect optional arguments from the ConfigurationPage
                optional_args = filter_optional_args({
                    "CREATE_INDEX": config_page.sortsam_for_merged_bam_create_index.get(),
                    "CREATE_MD5_FILE": config_page.sortsam_for_merged_bam_create_md5_file.get()
                })

                # Run the SortSam command with the provided inputs and optional arguments
                sortsam_result = process_management.run_sortsam(
                    gatk_path,  # Path to the GATK executable
                    merge_bam_output,  # Path to the input merged BAM file
                    sorted_output,  # Path to the output sorted BAM file
                    sort_order,  # The desired sort order
                    optional_args,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the SortSam process
                self.log_text.insert(tk.END, sortsam_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the SortSam process fails
                self.log_text.insert(tk.END, f"SortSam error: {str(e)}\n")

            # SetNmMdAndUqTags Process
            try:
                # Log the start of the SetNmMdAndUqTags process
                self.log_text.insert(tk.END, "Running SetNmMdAndUqTags...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the path for the output BAM file with fixed tags
                output_bam_fixed = os.path.join(output_directory, config_page.setnmmduqtags_output.get())

                # Collect optional arguments from the ConfigurationPage
                optional_args = filter_optional_args({
                    "IS_BISULFITE_SEQUENCE": config_page.setnmmduqtags_is_bisulfite_sequence.get(),
                    "SET_ONLY_UQ": config_page.set_only_uq.get(),
                    "CREATE_INDEX": config_page.setnmmduqtags_create_index.get(),
                    "CREATE_MD5_FILE": config_page.setnmmduqtags_create_md5_file.get()
                })

                # Run the SetNmMdAndUqTags command with the provided inputs and optional arguments
                setnmmduqtags_result = process_management.run_setnmmduqtags(
                    gatk_path,  # Path to the GATK executable
                    sorted_output,  # Path to the sorted BAM file (input)
                    output_bam_fixed,  # Path to the output BAM file with fixed tags
                    reference_sequence_path,  # Path to the reference sequence file
                    optional_args,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the SetNmMdAndUqTags process
                self.log_text.insert(tk.END, setnmmduqtags_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the SetNmMdAndUqTags process fails
                self.log_text.insert(tk.END, f"SetNmMdAndUqTags error: {str(e)}\n")

            # MarkDuplicates Process
            try:
                # Log the start of the MarkDuplicates process
                self.log_text.insert(tk.END, "Running MarkDuplicates...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the path for the metrics file and the output BAM file with duplicates marked
                metrics_file = os.path.join(output_directory, config_page.metrics_file.get())
                marked_output_bam = os.path.join(output_directory, config_page.mark_duplicates_output.get())

                # Collect optional arguments from the ConfigurationPage
                optional_args = filter_optional_args({
                    "ASSUME_SORT_ORDER": config_page.assume_sort_order.get(),
                    "BARCODE_TAG": config_page.barcode_tag.get(),
                    "CLEAR_DT": config_page.clear_dt.get(),
                    "COMMENT": config_page.markduplicates_comment.get(),
                    "DUPLICATE_SCORING_STRATEGY": config_page.duplicate_scoring_strategy.get(),
                    "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP": config_page.max_file_handles_for_read_ends_map.get(),
                    "MAX_OPTICAL_DUPLICATE_SET_SIZE": config_page.max_optical_duplicate_set_size.get(),
                    "MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP": config_page.max_sequences_for_disk_read_ends_map.get(),
                    "OPTICAL_DUPLICATE_PIXEL_DISTANCE": config_page.optical_duplicate_pixel_distance.get(),
                    "PROGRAM_GROUP_COMMAND_LINE": config_page.program_group_command_line.get(),
                    "PROGRAM_GROUP_NAME": config_page.program_group_name.get(),
                    "PROGRAM_GROUP_VERSION": config_page.program_group_version.get(),
                    "PROGRAM_RECORD_ID": config_page.program_record_id.get(),
                    "READ_NAME_REGEX": config_page.read_name_regex.get(),
                    "READ_ONE_BARCODE_TAG": config_page.read_one_barcode_tag.get(),
                    "READ_TWO_BARCODE_TAG": config_page.read_two_barcode_tag.get(),
                    "SORTING_COLLECTION_SIZE_RATIO": config_page.sorting_collection_size_ratio.get(),
                    "TAG_DUPLICATE_SET_MEMBERS": config_page.tag_duplicate_set_members.get(),
                    "TAGGING_POLICY": config_page.tagging_policy.get(),
                    "CREATE_INDEX": config_page.markduplicates_create_index.get(),
                    "CREATE_MD5_FILE": config_page.markduplicates_create_md5_file.get()
                })

                # Run the MarkDuplicates command with the provided inputs and optional arguments
                markduplicates_result = process_management.run_markduplicates(
                    gatk_path,  # Path to the GATK executable
                    output_bam_fixed,  # Path to the input BAM file with NM, MD, and UQ tags set
                    marked_output_bam,  # Path to the output BAM file with duplicates marked
                    metrics_file,  # Path to the file where duplication metrics will be saved
                    optional_args,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the MarkDuplicates process
                self.log_text.insert(tk.END, markduplicates_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the MarkDuplicates process fails
                self.log_text.insert(tk.END, f"MarkDuplicates error: {str(e)}\n")

            # SortSam for Marked Duplicate BAM Process
            try:
                # Log the start of the SortSam process for the marked duplicate BAM file
                self.log_text.insert(tk.END, "Running SortSam for Marked Duplicate BAM...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Retrieve the sort order and output path from the configuration page
                sort_order_marked = config_page.sortsam_for_marked_duplicate_sort_order.get()
                sorted_output_marked = os.path.join(output_directory,
                                                    config_page.sortsam_for_marked_duplicate_output.get())

                # Collect optional arguments for sorting
                optional_args_marked = {
                    "CREATE_INDEX": config_page.sortsam_for_marked_duplicate_create_index.get(),
                    "CREATE_MD5_FILE": config_page.sortsam_for_marked_duplicate_create_md5_file.get()
                }

                # Run the SortSam command for the marked duplicate BAM file
                sortsam_marked_result = process_management.run_sortsam_marked(
                    gatk_path,  # Path to the GATK executable
                    marked_output_bam,  # Path to the BAM file with duplicates marked
                    sorted_output_marked,  # Path to the output sorted BAM file
                    sort_order_marked,  # The desired sort order (e.g., coordinate)
                    optional_args_marked,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the SortSam process
                self.log_text.insert(tk.END, sortsam_marked_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the SortSam process fails
                self.log_text.insert(tk.END, f"SortSam for Marked Duplicate BAM error: {str(e)}\n")

            # BaseRecalibrator Process
            try:
                # Log the start of the BaseRecalibrator process
                self.log_text.insert(tk.END, "Running BaseRecalibrator...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Retrieve the output path for the recalibration table and known sites
                base_recal_output = os.path.join(output_directory, config_page.base_recal_output.get())
                known_sites = config_page.base_recal_known_sites.get()

                # Collect optional arguments for the BaseRecalibrator process
                optional_args_base_recal = {
                    "bqsr-baq-gap-open-penalty": config_page.bqsr_baq_gap_open_penalty.get(),
                    "cloud-index-prefetch-buffer": config_page.cloud_index_prefetch_buffer.get(),
                    "cloud-prefetch-buffer": config_page.cloud_prefetch_buffer.get(),
                    "default-base-qualities": config_page.default_base_qualities.get(),
                    "deletions-default-quality": config_page.deletions_default_quality.get(),
                    "disable-bam-index-caching": str(config_page.disable_bam_index_caching.get()).lower(),
                    "gcs-max-retries": config_page.gcs_max_retries.get(),
                    "indels-context-size": config_page.indels_context_size.get(),
                    "insertions-default-quality": config_page.insertions_default_quality.get(),
                    "interval-merging-rule": config_page.interval_merging_rule.get(),
                    "intervals": config_page.intervals.get(),
                    "low-quality-tail": config_page.low_quality_tail.get(),
                    "maximum-cycle-value": config_page.maximum_cycle_value.get(),
                    "mismatches-context-size": config_page.mismatches_context_size.get(),
                    "mismatches-default-quality": config_page.mismatches_default_quality.get(),
                    "preserve-qscores-less-than": config_page.preserve_qscores_less_than.get(),
                    "quantizing-levels": config_page.quantizing_levels.get(),
                    "use-original-qualities": str(config_page.use_original_qualities.get()).lower()
                }

                # Run the BaseRecalibrator command with the provided inputs and optional arguments
                base_recal_result = process_management.run_baserecalibrator(
                    gatk_path,  # Path to the GATK executable
                    sorted_output_marked,  # Input BAM file from the SortSam for Marked Duplicate BAM process
                    reference_sequence_path,  # Path to the reference sequence file
                    base_recal_output,  # Path to the output recalibration table file
                    known_sites,  # Known sites for recalibration
                    optional_args_base_recal,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )
                

                # Log the result of the BaseRecalibrator process
                self.log_text.insert(tk.END, base_recal_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the BaseRecalibrator process fails
                self.log_text.insert(tk.END, f"BaseRecalibrator error: {str(e)}\n")

            # ApplyBQSR Process
            try:
                # Log the start of the ApplyBQSR process
                self.log_text.insert(tk.END, "Running ApplyBQSR...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Retrieve the output path for the recalibrated BAM file and the BQSR recalibration file
                applybqsr_output = os.path.join(output_directory, config_page.applybqsr_output.get())
                bqsr_recal_file = os.path.join(output_directory, config_page.base_recal_output.get())

                # Collect optional arguments for the ApplyBQSR process
                optional_args_applybqsr = {
                    "cloud-index-prefetch-buffer": config_page.applybqsr_cloud_index_prefetch_buffer.get(),
                    "cloud-prefetch-buffer": config_page.applybqsr_cloud_prefetch_buffer.get(),
                    "disable-bam-index-caching": str(config_page.applybqsr_disable_bam_index_caching.get()).lower(),
                    "emit-original-quals": str(config_page.applybqsr_emit_original_quals.get()).lower(),
                    "gcs-max-retries": config_page.applybqsr_gcs_max_retries.get(),
                    "global-qscore-prior": config_page.applybqsr_global_qscore_prior.get(),
                    "interval-merging-rule": config_page.applybqsr_interval_merging_rule.get(),
                    "intervals": config_page.applybqsr_intervals.get(),
                    "preserve-qscores-less-than": config_page.applybqsr_preserve_qscores_less_than.get(),
                    "quantize-quals": config_page.applybqsr_quantize_quals.get(),
                    "use-original-qualities": str(config_page.applybqsr_use_original_qualities.get()).lower()
                }

                # Run the ApplyBQSR command with the provided inputs and optional arguments
                applybqsr_result = process_management.run_applybqsr(
                    gatk_path,  # Path to the GATK executable
                    sorted_output_marked,  # Input BAM file from the SortSam for Marked Duplicate BAM process
                    applybqsr_output,  # Path to the output BAM file with recalibrated base quality scores
                    reference_sequence_path,  # Path to the reference sequence file
                    bqsr_recal_file,  # Path to the BQSR recalibration table file
                    optional_args_applybqsr,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the ApplyBQSR process
                self.log_text.insert(tk.END, applybqsr_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the ApplyBQSR process fails
                self.log_text.insert(tk.END, f"ApplyBQSR error: {str(e)}\n")

            # Collect alignment summary using samtools flagstat
            try:
                # Log the start of the samtools flagstat process
                self.log_text.insert(tk.END, "Collecting alignment summary with samtools flagstat...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the output path for the flagstat report
                flagstat_output = os.path.join(output_directory, "bqsr_flagstat.txt")

                # Run the samtools flagstat command to generate the alignment summary
                flagstat_result = process_management.run_samtools_flagstat(
                    applybqsr_output,  # Path to the BAM file after ApplyBQSR
                    flagstat_output,  # Path to the output flagstat report
                    log_function  # Function to log messages during the process
                )

                # Log the result of the flagstat process
                self.log_text.insert(tk.END, flagstat_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the stat process fails
                self.log_text.insert(tk.END, f"Samtools stat error: {str(e)}\n\n")

            # Collect alignment summary using samtools stat
            try:
                # Log the start of the samtools stat process
                self.log_text.insert(tk.END, "Collecting alignment statistic with samtools stat...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the output path for the stat report
                stat_output = os.path.join(output_directory, "bqsr_stat.txt")

                # Run the samtools stat command to generate the alignment statistic
                stat_result = process_management.run_samtools_stats(
                    applybqsr_output,  # Path to the BAM file after ApplyBQSR
                    stat_output,  # Path to the output stat report
                    log_function  # Function to log messages during the process
                )

                # Log the result of the stat process
                self.log_text.insert(tk.END, stat_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the stat process fails
                self.log_text.insert(tk.END, f"Samtools stat error: {str(e)}\n")


            # HaplotypeCaller Process
            try:
                # Log the start of the HaplotypeCaller process
                self.log_text.insert(tk.END, "Running HaplotypeCaller...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the output path for the variant call file (VCF) generated by HaplotypeCaller
                haplotypecaller_output = os.path.join(output_directory, config_page.haplotypecaller_output.get())

                # Collect optional arguments for the HaplotypeCaller process
                optional_args_haplotypecaller = filter_optional_args({
                    "activity-profile-out": config_page.haplotypecaller_activity_profile_out.get(),
                    "alleles": config_page.haplotypecaller_alleles.get(),
                    "annotate-with-num-discovered-alleles": config_page.haplotypecaller_annotate_with_num_discovered_alleles.get(),
                    "annotation": config_page.haplotypecaller_annotation.get(),
                    "annotation-group": config_page.haplotypecaller_annotation_group.get(),
                    "annotations-to-exclude": config_page.haplotypecaller_annotations_to_exclude.get(),
                    "assembly-region-out": config_page.haplotypecaller_assembly_region_out.get(),
                    "base-quality-score-threshold": config_page.haplotypecaller_base_quality_score_threshold.get(),
                    "cloud-index-prefetch-buffer": config_page.haplotypecaller_cloud_index_prefetch_buffer.get(),
                    "cloud-prefetch-buffer": config_page.haplotypecaller_cloud_prefetch_buffer.get(),
                    "contamination-fraction-to-filter": config_page.haplotypecaller_contamination_fraction_to_filter.get(),
                    "correct-overlapping-quality": config_page.haplotypecaller_correct_overlapping_quality.get(),
                    "dbsnp": config_page.haplotypecaller_dbsnp.get(),
                    "disable-bam-index-caching": config_page.haplotypecaller_disable_bam_index_caching.get(),
                    "disable-sequence-dictionary-validation": config_page.haplotypecaller_disable_sequence_dictionary_validation.get(),
                    "founder-id": config_page.haplotypecaller_founder_id.get(),
                    "gcs-max-retries": config_page.haplotypecaller_gcs_max_retries.get(),
                    "gcs-project-for-requester-pays": config_page.haplotypecaller_gcs_project_for_requester_pays.get(),
                    "graph-output": config_page.haplotypecaller_graph_output.get(),
                    "heterozygosity": config_page.haplotypecaller_heterozygosity.get(),
                    "heterozygosity-stdev": config_page.haplotypecaller_heterozygosity_stdev.get(),
                    "indel-heterozygosity": config_page.haplotypecaller_indel_heterozygosity.get(),
                    "interval-merging-rule": config_page.haplotypecaller_interval_merging_rule.get(),
                    "intervals": config_page.haplotypecaller_intervals.get(),
                    "max-reads-per-alignment-start": config_page.haplotypecaller_max_reads_per_alignment_start.get(),
                    "min-base-quality-score": config_page.haplotypecaller_min_base_quality_score.get(),
                    "native-pair-hmm-threads": config_page.haplotypecaller_native_pair_hmm_threads.get(),
                    "native-pair-hmm-use-double-precision": config_page.haplotypecaller_native_pair_hmm_use_double_precision.get(),
                    "num-reference-samples-if-no-call": config_page.haplotypecaller_num_reference_samples_if_no_call.get(),
                    "output-mode": config_page.haplotypecaller_output_mode.get(),
                    "pedigree": config_page.haplotypecaller_pedigree.get(),
                    "population-callset": config_page.haplotypecaller_population_callset.get(),
                    "sample-name": config_page.haplotypecaller_sample_name.get(),
                    "sample-ploidy": config_page.haplotypecaller_sample_ploidy.get(),
                    "sites-only-vcf-output": config_page.haplotypecaller_sites_only_vcf_output.get(),
                    "standard-min-confidence-threshold-for-calling": config_page.haplotypecaller_stand_call_conf.get()
                })

                # Run the HaplotypeCaller command with the provided inputs and optional arguments
                haplotypecaller_result = process_management.run_haplotypecaller(
                    gatk_path,  # Path to the GATK executable
                    applybqsr_output,  # Input BAM file from the ApplyBQSR process
                    haplotypecaller_output,  # Output VCF file path
                    reference_sequence_path,  # Path to the reference sequence file
                    optional_args_haplotypecaller,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the HaplotypeCaller process
                self.log_text.insert(tk.END, haplotypecaller_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the HaplotypeCaller process fails
                self.log_text.insert(tk.END, f"HaplotypeCaller error: {str(e)}\n")

            # SNV Extraction Process
            try:
                # Log the start of the SNV Extraction process
                self.log_text.insert(tk.END, "Running SNV Extraction by SelectVariants...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the input and output VCF file paths
                input_vcf = os.path.join(output_directory, config_page.haplotypecaller_output.get())
                output_vcf = os.path.join(output_directory, config_page.snv_output.get())

                # Collect optional arguments for the SNV Extraction process
                optional_args_snv = filter_optional_args({
                    "cloud-index-prefetch-buffer": config_page.snv_cloud_index_prefetch_buffer.get(),
                    "cloud-prefetch-buffer": config_page.snv_cloud_prefetch_buffer.get(),
                    "concordance": config_page.snv_concordance.get(),
                    "disable-bam-index-caching": config_page.snv_disable_bam_index_caching.get(),
                    "discordance": config_page.snv_discordance.get(),
                    "exclude-filtered": config_page.snv_exclude_filtered.get(),
                    "exclude-ids": config_page.snv_exclude_ids.get(),
                    "exclude-non-variants": config_page.snv_exclude_non_variants.get(),
                    "exclude-sample-expressions": config_page.snv_exclude_sample_expressions.get(),
                    "exclude-sample-name": config_page.snv_exclude_sample_name.get(),
                    "gcs-max-retries": config_page.snv_gcs_max_retries.get(),
                    "interval-merging-rule": config_page.snv_interval_merging_rule.get(),
                    "intervals": config_page.snv_intervals.get(),
                    "invert-mendelian-violation": config_page.snv_invert_mendelian_violation.get(),
                    "invert-select": config_page.snv_invert_select.get(),
                    "keep-ids": config_page.snv_keep_ids.get(),
                    "keep-original-ac": config_page.snv_keep_original_ac.get(),
                    "keep-original-dp": config_page.snv_keep_original_dp.get(),
                    "max-filtered-genotypes": config_page.snv_max_filtered_genotypes.get(),
                    "max-fraction-filtered-genotypes": config_page.snv_max_fraction_filtered_genotypes.get(),
                    "max-indel-size": config_page.snv_max_indel_size.get(),
                    "max-nocall-fraction": config_page.snv_max_nocall_fraction.get(),
                    "max-nocall-number": config_page.snv_max_nocall_number.get(),
                    "mendelian-violation": config_page.snv_mendelian_violation.get(),
                    "mendelian-violation-qual-threshold": config_page.snv_mendelian_violation_qual_threshold.get(),
                    "min-filtered-genotypes": config_page.snv_min_filtered_genotypes.get(),
                    "min-fraction-filtered-genotypes": config_page.snv_min_fraction_filtered_genotypes.get(),
                    "min-indel-size": config_page.snv_min_indel_size.get(),
                    "pedigree": config_page.snv_pedigree.get(),
                    "preserve-alleles": config_page.snv_preserve_alleles.get(),
                    "remove-fraction-genotypes": config_page.snv_remove_fraction_genotypes.get(),
                    "remove-unused-alternates": config_page.snv_remove_unused_alternates.get(),
                    "restrict-alleles-to": config_page.snv_restrict_alleles_to.get(),
                    "sample-expressions": config_page.snv_sample_expressions.get(),
                    "sample-name": config_page.snv_sample_name.get(),
                    "select-random-fraction": config_page.snv_select_random_fraction.get(),
                    "select-type-to-exclude": config_page.snv_select_type_to_exclude.get(),
                    "select-type-to-include": config_page.snv_select_type_to_include.get(),
                    "select": config_page.snv_select_expressions.get(),
                    "set-filtered-gt-to-nocall": config_page.snv_set_filtered_gt_to_nocall.get()
                })

                # Run the SNV Extraction command with the provided inputs and optional arguments
                snv_extraction_result = process_management.run_snvextraction(
                    gatk_path,  # Path to the GATK executable
                    input_vcf,  # Input VCF file from the HaplotypeCaller output
                    output_vcf,  # Output VCF file path for the extracted SNVs
                    reference_sequence_path,  # Path to the reference sequence file
                    optional_args_snv,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the SNV Extraction process
                self.log_text.insert(tk.END, snv_extraction_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the SNV Extraction process fails
                self.log_text.insert(tk.END, f"SNV Extraction error: {str(e)}\n")

            # Indel Extraction Process
            try:
                # Log the start of the Indel Extraction process
                self.log_text.insert(tk.END, "Running Indel Extraction by SelectVariants...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the input and output VCF file paths
                haplotypecaller_output = os.path.join(output_directory, config_page.haplotypecaller_output.get())
                indel_output = os.path.join(output_directory, config_page.indel_output.get())

                # Collect optional arguments for the Indel Extraction process
                optional_args_indel = filter_optional_args({
                    "cloud-index-prefetch-buffer": config_page.indel_cloud_index_prefetch_buffer.get(),
                    "cloud-prefetch-buffer": config_page.indel_cloud_prefetch_buffer.get(),
                    "concordance": config_page.indel_concordance.get(),
                    "disable-bam-index-caching": str(config_page.indel_disable_bam_index_caching.get()).lower(),
                    "discordance": config_page.indel_discordance.get(),
                    "exclude-filtered": str(config_page.indel_exclude_filtered.get()).lower(),
                    "exclude-ids": config_page.indel_exclude_ids.get(),
                    "exclude-non-variants": str(config_page.indel_exclude_non_variants.get()).lower(),
                    "exclude-sample-expressions": config_page.indel_exclude_sample_expressions.get(),
                    "exclude-sample-name": config_page.indel_exclude_sample_name.get(),
                    "gcs-max-retries": config_page.indel_gcs_max_retries.get(),
                    "interval-merging-rule": config_page.indel_interval_merging_rule.get(),
                    "intervals": config_page.indel_intervals.get(),
                    "invert-mendelian-violation": str(config_page.indel_invert_mendelian_violation.get()).lower(),
                    "invert-select": str(config_page.indel_invert_select.get()).lower(),
                    "keep-ids": config_page.indel_keep_ids.get(),
                    "keep-original-ac": str(config_page.indel_keep_original_ac.get()).lower(),
                    "keep-original-dp": str(config_page.indel_keep_original_dp.get()).lower(),
                    "max-filtered-genotypes": config_page.indel_max_filtered_genotypes.get(),
                    "max-fraction-filtered-genotypes": config_page.indel_max_fraction_filtered_genotypes.get(),
                    "max-indel-size": config_page.indel_max_indel_size.get(),
                    "max-nocall-fraction": config_page.indel_max_nocall_fraction.get(),
                    "max-nocall-number": config_page.indel_max_nocall_number.get(),
                    "mendelian-violation": str(config_page.indel_mendelian_violation.get()).lower(),
                    "mendelian-violation-qual-threshold": config_page.indel_mendelian_violation_qual_threshold.get(),
                    "min-filtered-genotypes": config_page.indel_min_filtered_genotypes.get(),
                    "min-fraction-filtered-genotypes": config_page.indel_min_fraction_filtered_genotypes.get(),
                    "min-indel-size": config_page.indel_min_indel_size.get(),
                    "pedigree": config_page.indel_pedigree.get(),
                    "preserve-alleles": str(config_page.indel_preserve_alleles.get()).lower(),
                    "remove-fraction-genotypes": config_page.indel_remove_fraction_genotypes.get(),
                    "remove-unused-alternates": str(config_page.indel_remove_unused_alternates.get()).lower(),
                    "restrict-alleles-to": config_page.indel_restrict_alleles_to.get(),
                    "sample-expressions": config_page.indel_sample_expressions.get(),
                    "sample-name": config_page.indel_sample_name.get(),
                    "select-random-fraction": config_page.indel_select_random_fraction.get(),
                    "select-type-to-exclude": config_page.indel_select_type_to_exclude.get(),
                    "select-type-to-include": config_page.indel_select_type_to_include.get(),
                    "select": config_page.indel_select_expressions.get(),
                    "set-filtered-gt-to-nocall": str(config_page.indel_set_filtered_gt_to_nocall.get()).lower()
                })

                # Run the Indel Extraction command with the provided inputs and optional arguments
                indelextraction_result = process_management.run_indelextraction(
                    gatk_path,  # Path to the GATK executable
                    haplotypecaller_output,  # Input VCF from HaplotypeCaller
                    indel_output,  # Output VCF for Indel extraction
                    reference_sequence_path,  # Path to the reference sequence file
                    optional_args_indel,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the Indel Extraction process
                self.log_text.insert(tk.END, indelextraction_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the Indel Extraction process fails
                self.log_text.insert(tk.END, f"Indel Extraction error: {str(e)}\n")

            # CollectVariantCallingMetrics Process
            try:
                # Log the start of the CollectVariantCallingMetrics process
                self.log_text.insert(tk.END, "Running CollectVariantCallingMetrics...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Define the output path for the metrics file
                collectvariantcallingmetrics_output = os.path.join(output_directory,
                                                                   config_page.collectvariantcallingmetrics_output.get())

                # Define the path to the dbSNP file for known variants
                dbsnp_file = config_page.collectvariantcallingmetrics_dbsnp.get()

                # Collect optional arguments for the CollectVariantCallingMetrics process
                optional_args_collectvariantcallingmetrics = filter_optional_args({
                    "GVCF_INPUT": config_page.collectvariantcallingmetrics_gvcf_input.get(),
                    "SEQUENCE_DICTIONARY": config_page.collectvariantcallingmetrics_sequence_dictionary.get(),
                    "TARGET_INTERVALS": config_page.collectvariantcallingmetrics_target_intervals.get(),
                    "THREAD_COUNT": config_page.collectvariantcallingmetrics_thread_count.get()
                })

                # Run the CollectVariantCallingMetrics command with the provided inputs and optional arguments
                collectvariantcallingmetrics_result = process_management.run_collectvariantcallingmetrics(
                    gatk_path,  # Path to the GATK executable
                    haplotypecaller_output,  # Input VCF from HaplotypeCaller
                    collectvariantcallingmetrics_output,  # Output path for metrics
                    dbsnp_file,  # dbSNP file for known variants
                    optional_args_collectvariantcallingmetrics,  # Dictionary of optional arguments
                    log_function  # Function to log messages during the process
                )

                # Log the result of the CollectVariantCallingMetrics process
                self.log_text.insert(tk.END, collectvariantcallingmetrics_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the CollectVariantCallingMetrics process fails
                self.log_text.insert(tk.END, f"CollectVariantCallingMetrics error: {str(e)}\n")

            # ANNOVAR Process
            try:
                # Log the start of the ANNOVAR process
                self.log_text.insert(tk.END, "Running ANNOVAR...\n")
                self.update_idletasks()  # Update the GUI to show the log message immediately

                # Retrieve the path to the ANNOVAR executable
                annovar_path = program_page.annovar_location.entry.get()

                # Define the input VCF files for annotation
                input_vcfs = [
                    os.path.join(output_directory, config_page.haplotypecaller_output.get()),
                    # VCF from HaplotypeCaller
                    os.path.join(output_directory, config_page.snv_output.get()),  # VCF with extracted SNVs
                    os.path.join(output_directory, config_page.indel_output.get())  # VCF with extracted Indels
                ]

                # Get the database directory and genome build version
                db_dir = config_page.database_directory.get()  # Directory containing ANNOVAR databases
                buildver = config_page.genome_build_version.get()  # Genome build version (e.g., hg19, hg38)
                output_prefix = config_page.output_prefix.get()  # Prefix for the output files

                # Collect optional arguments for the ANNOVAR process
                optional_args_annovar = filter_optional_args({
                    "protocol": config_page.protocol.get(),  # Protocols for annotation (e.g., refGene, cytoBand)
                    "operation": config_page.operation.get(),  # Operations corresponding to the protocols
                    "nastring": config_page.nastring.get(),  # String to use for missing values
                    "nopolish": config_page.no_polish.get(),  # Whether to disable post-processing
                    "remove": config_page.remove_intermediate_files.get()  # Whether to remove intermediate files
                })

                # Run ANNOVAR for each input VCF file
                for input_vcf in input_vcfs:
                    # Define the output prefix for the current VCF file
                    output_prefix_vcf = os.path.join(output_directory,
                                                     f"{output_prefix}_{os.path.basename(input_vcf).split('.')[0]}")

                    # Run the ANNOVAR command with the provided inputs and optional arguments
                    annovar_result = process_management.run_annovar(
                        annovar_path,  # Path to the ANNOVAR executable
                        input_vcf,  # Input VCF file
                        output_prefix_vcf,  # Prefix for the output files
                        buildver,  # Genome build version
                        db_dir,  # Database directory for ANNOVAR
                        optional_args_annovar,  # Dictionary of optional arguments
                        log_function  # Function to log messages during the process
                    )

                    # Log the result of the ANNOVAR process for the current VCF file
                    self.log_text.insert(tk.END, annovar_result + "\n")
            except Exception as e:
                # Handle exceptions and log the error if the ANNOVAR process fails
                self.log_text.insert(tk.END, f"ANNOVAR error: {str(e)}\n")

            # Indicate the completion of the analysis
            self.log_text.insert(tk.END, "Analysis completed.\n")
            self.log_text.yview(tk.END)  # Scroll the log to the end to show the latest message

        except Exception as e:
            # Handle any unexpected exceptions that occur during the process
            self.log_text.insert(tk.END, f"Unexpected error: {str(e)}\n")
            self.update_idletasks() # Update the GUI to show the error message immediately

    def create_browse_section(self, label_text, browse_command):
        # Create a frame to contain the label, entry box, and browse button
        frame = ttk.Frame(self)

        # Create a label with the specified text
        label = ttk.Label(frame, text=label_text, width=20, anchor="w")
        label.pack(side="left", padx=5)  # Place the label on the left side of the frame with padding

        # Create an entry box where the user can type or see the selected path
        entry = tk.Entry(frame, width=50)
        entry.pack(side="left", fill="x", expand=True, padx=5)  # Allow the entry box to expand horizontally

        # Create a "Browse" button that calls the provided browse_command with the entry as an argument
        button = ttk.Button(frame, text="Browse", command=lambda: browse_command(entry))
        button.pack(side="left", padx=5)  # Place the button on the left side of the frame with padding

        # Attach the entry widget to the frame for later reference
        frame.entry = entry
        return frame  # Return the complete frame with its children

    def browse_directory(self, entry):
        # Open a dialog for selecting a directory
        directory = filedialog.askdirectory()
        if directory:
            # If a directory was selected, update the entry box with the directory path
            entry.delete(0, tk.END)  # Clear the current contents of the entry box
            entry.insert(0, directory)  # Insert the selected directory path into the entry box

    def browse_file(self, entry):
        # Open a dialog for selecting a file and get the selected file's path
        file_path = filedialog.askopenfilename()
        output_directory = self.output_directory.entry.get()  # Retrieve the output directory from another entry

        # Ensure file_path is a string, not a tuple (older tkinter versions might return a tuple)
        if isinstance(file_path, tuple):
            file_path = file_path[
                0] if file_path else ""  # Use the first file if multiple selected, or empty string if none.

        if file_path:  # Proceed only if file_path is not empty.
            # Check if the selected file is within the output directory
            if output_directory and not file_path.startswith(output_directory):
                # Show an error message if the file is not in the output directory
                tk.messagebox.showerror("Invalid Selection",
                                        f"Please select a file from the directory: {output_directory}")
            else:
                # Clear the current entry and insert the selected file path
                entry.delete(0, tk.END)
                entry.insert(0, file_path)
        else:
            # Inform the user that no file was selected
            tk.messagebox.showinfo("No File Selected", "No file was selected.")

    def browse_reference_file(self, entry):
        # Open a dialog for selecting a FASTA file
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fa")])
        if file_path:
            # Clear the current entry and insert the selected file path
            entry.delete(0, tk.END)
            entry.insert(0, file_path)

    def is_file_in_directory(self, file_path, directory):
        # Get absolute paths for both the directory and the file
        abs_directory = os.path.abspath(directory)
        abs_file_path = os.path.abspath(file_path)
        # Check if the file's absolute path starts with the directory's absolute path
        return abs_file_path.startswith(abs_directory)


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
                ("FASTQ to SAM Output", "unaligned_read_pairs.bam", False, tk.StringVar),
                ("Sample Name", "sample001", False, tk.StringVar),
                ("Platform", "", False, tk.StringVar),
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
                ("SAM to FASTQ Output", "output.fastq", False, tk.StringVar)
            ]),
            ("Optional Arguments for SamToFastq", [
                ("Clipping Action", "", False, tk.StringVar),
                ("Clipping Attribute", "", False, tk.StringVar),
                ("Clipping Min Length", 0, False, tk.IntVar),
                ("Include Non PF Reads", False, True, tk.BooleanVar),
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
                ("Threads", 4, False, tk.IntVar),
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
                ("Add Mate Cigar", False, True, tk.BooleanVar),
                ("Aligned Reads Only", False, True, tk.BooleanVar),
                ("Aligner Proper Pair Flags", False, True, tk.BooleanVar),
                ("Attributes To Remove", "", False, tk.StringVar),
                ("Attributes To Retain", "", False, tk.StringVar),
                ("Attributes To Reverse", "", False, tk.StringVar),
                ("Attributes To Reverse Complement", "", False, tk.StringVar),
                ("Clip Adapters", True, True, tk.BooleanVar),
                ("Clip Overlapping Reads", True, True, tk.BooleanVar),
                ("Expected Orientations", "", False, tk.StringVar),
                ("Include Secondary Alignments", True, True, tk.BooleanVar),
                ("MergeBamAlignment Is Bisulfite Sequence", False, True, tk.BooleanVar),
                ("Matching Dictionary Tags", "", False, tk.StringVar),
                ("Max Insertions Or Deletions", 1, False, tk.IntVar),
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
                ("Unmap Contaminant Reads", False, True, tk.BooleanVar),
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
                ("SetnmmdUqtags Output", "fixed.bam")
            ]),
            ("Optional Arguments for SetNmMdAndUqTags", [
                ("SetnmmdUqtags Is Bisulfite Sequence", False, True),
                ("Set Only Uq", False, True),
                ("SetnmmdUqtags Create Index", False, True),
                ("SetnmmdUqtags Create MD5 File", False, True)
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
                ("MarkDuplicates Create MD5 File", False, True)
            ])
        ])

        # SortSam (Picard) For Marked Duplicate BAM section
        self.create_picard_section(scrollable_frame, "SortSam (Picard) For Marked Duplicate BAM", [
            ("Required Arguments for SortSam", [
                ("SortSam For Marked Duplicate Output", "sorted_marked.bam", False, tk.StringVar),
                ("SortSam For Marked Duplicate Sort Order", "coordinate", False, tk.StringVar)
            ]),
            ("Optional Arguments for SortSam", [
                ("SortSam For Marked Duplicate Create Index", False, True),
                ("SortSam For Marked Duplicate Create MD5 File", False, True)
            ])
        ])


        # BaseRecalibrator section
        self.create_picard_section(scrollable_frame, "BaseRecalibrator (Picard)", [
            ("Required Arguments for BaseRecalibrator", [
                ("Base Recal Output", "recal_data.table", False, tk.StringVar),
                ("Base Recal Known Sites", "", False, tk.StringVar)
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
                ("ApplyBQSR Use Original Qualities", False, True)
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
                ("Collectvariantcallingmetrics Dbsnp", "", False, tk.StringVar)
            ]),
            ("Optional Arguments for CollectVariantCallingMetrics", [
                ("Collectvariantcallingmetrics Gvcf Input", False, True),
                ("Collectvariantcallingmetrics Sequence Dictionary", "", False, tk.StringVar),
                ("Collectvariantcallingmetrics Target Intervals", "", False, tk.StringVar),
                ("Collectvariantcallingmetrics Thread Count", 10, False, tk.IntVar)
            ])
        ])

        # ANNOVAR Section
        self.create_picard_section(scrollable_frame, "ANNOVAR", [
            ("Required Arguments for ANNOVAR", [
                ("Database Directory", "", False, tk.StringVar),
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


class ProgramPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        # Create a label for the page title
        label = tk.Label(self, text="Program Page", font=("Arial", 24))
        label.pack(pady=10)  # Add some padding around the label

        # Create sections for browsing and selecting program locations
        self.fastqc_location = self.create_browse_section("FASTQC Location", self.browse_file)
        self.fastqc_location.pack(fill="x", padx=50, pady=5)

        self.gatk_location = self.create_browse_section("GATK Location", self.browse_file)
        self.gatk_location.pack(fill="x", padx=50, pady=5)

        self.bwa_location = self.create_browse_section("BWA Location", self.browse_file)
        self.bwa_location.pack(fill="x", padx=50, pady=5)

        self.annovar_location = self.create_browse_section("ANNOVAR Location", self.browse_directory)
        self.annovar_location.pack(fill="x", padx=50, pady=5)

        # Create a frame for Save and Load buttons
        button_frame = ttk.Frame(self)
        button_frame.pack(pady=10)

        # Save button to save the current configuration as default
        save_button = ttk.Button(button_frame, text="Save as Default", command=self.save_configuration)
        save_button.pack(side="left", padx=10)

        # Load button to load the default configuration
        load_button = ttk.Button(button_frame, text="Load Default", command=self.load_configuration)
        load_button.pack(side="left", padx=10)

    def create_browse_section(self, label_text, browse_function):
        # Create a frame to hold the label, entry, and browse button
        frame = ttk.Frame(self)

        # Create a label for the section
        label = ttk.Label(frame, text=label_text, width=20, anchor="w")
        label.pack(side="left", padx=5)

        # Create an entry widget for displaying the selected path
        entry = tk.Entry(frame, width=50)
        entry.pack(side="left", fill="x", expand=True, padx=5)

        # Create a browse button that triggers the provided browse_function
        button = ttk.Button(frame, text="Browse", command=lambda: browse_function(entry))
        button.pack(side="left", padx=5)

        # Store the entry widget in the frame for later access
        frame.entry = entry
        return frame

    def browse_file(self, entry):
        # Open a file selection dialog
        file_path = filedialog.askopenfilename()
        if file_path:
            # If a file was selected, update the entry widget with the file path
            entry.delete(0, tk.END)
            entry.insert(0, file_path)

    def browse_directory(self, entry):
        # Open a directory selection dialog
        directory_path = filedialog.askdirectory()
        if directory_path:
            # If a directory was selected, update the entry widget with the directory path
            entry.delete(0, tk.END)
            entry.insert(0, directory_path)

    def save_configuration(self):
        # Gather the configuration data from the entry widgets
        config = {
            'fastqc_location': self.fastqc_location.entry.get(),
            'gatk_location': self.gatk_location.entry.get(),
            'bwa_location': self.bwa_location.entry.get(),
            'annovar_location': self.annovar_location.entry.get(),
        }

        # Define the configuration directory and ensure it exists
        config_directory = os.path.join(os.path.dirname(__file__), 'configuration')
        os.makedirs(config_directory, exist_ok=True)

        # Write the configuration data to a JSON file
        with open(os.path.join(config_directory, 'external_program_location.json'), 'w') as f:
            json.dump(config, f, indent=4)

        # Inform the user that the configuration has been saved
        tk.messagebox.showinfo("Save Configuration", "Configuration saved successfully!")

    def load_configuration(self):
        # Define the path to the configuration file
        config_path = os.path.join(os.path.dirname(__file__), 'configuration', 'external_program_location.json')

        # Check if the configuration file exists
        if os.path.exists(config_path):
            # Read the configuration file and load the data
            with open(config_path, 'r') as f:
                config = json.load(f)

            # Update the entry widgets with the loaded configuration data
            self.fastqc_location.entry.delete(0, tk.END)
            self.fastqc_location.entry.insert(0, config['fastqc_location'])
            self.gatk_location.entry.delete(0, tk.END)
            self.gatk_location.entry.insert(0, config['gatk_location'])
            self.bwa_location.entry.delete(0, tk.END)
            self.bwa_location.entry.insert(0, config['bwa_location'])
            self.annovar_location.entry.delete(0, tk.END)
            self.annovar_location.entry.insert(0, config['annovar_location'])

            # Inform the user that the configuration has been Loaded
            tk.messagebox.showinfo("Load Configuration", "Configuration loaded successfully!")
        else:
            # Show an error message if the configuration file does not exist
            tk.messagebox.showerror("Load Error", "No configuration file found.")


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
        self.create_section("Overview", overview_content, 0)
        self.create_section("Purpose", purpose_content, 1)
        self.create_section("Audience", audience_content, 2)

    def create_section(self, title, content, row):
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
        self.create_section("Software Requirements", software_requirements, 0)
        self.create_section("Hardware Requirements", hardware_requirements, 1)

    def create_section(self, title, content, row):
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
        self.create_section("Launching the Application", launching_content, 0)
        self.create_section("User Interface Overview", ui_overview_content, 1)

    def create_section(self, title, content, row):
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

        self.create_section("Workflow Overview", workflow_content, 0)

    def create_section(self, title, content, row):
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
        self.create_section("Common Issues", common_issues_content, 0)
        self.create_section("Logs and Debugging", logs_debugging_content, 1)

    def create_section(self, title, content, row):
        frame = tk.Frame(self, borderwidth=2, relief="groove")
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        label = tk.Label(frame, text=title, font=("Arial", 16, "bold"))
        label.pack(pady=(5, 10))

        text_widget = tk.Text(frame, wrap="word", height=10, width=80, state="disabled", bg="#f0f0f0", bd=0)
        text_widget.pack(padx=10, pady=10, fill="both", expand=True)
        text_widget.configure(state="normal")
        text_widget.insert("1.0", content)
        text_widget.configure(state="disabled")


    # Main Application Execution
if __name__ == "__main__":
    app = SequenceAssemblyApp()
    app.mainloop()
