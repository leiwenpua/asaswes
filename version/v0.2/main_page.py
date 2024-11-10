import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
import os
import json
import threading
import process_management


class MainPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        # Label indicating the Main Page
        label = tk.Label(self, text="Main Page", font=("Arial", 24))
        label.pack(pady=10)

        # Create and pack a section for selecting the folder containing the sequences
        self.sequence_folder = self.create_browse_section("Sequence Folder", self.browse_seq)
        self.sequence_folder.pack(fill="x", padx=50, pady=5)

        # Create and pack a section for selecting the reference sequence file
        self.reference_sequence = self.create_browse_section("Reference Genome", self.browse_reference_file)
        self.reference_sequence.pack(fill="x", padx=50, pady=5)

        # Create and pack a section for directory selection
        self.base_directory = self.create_browse_section("Directory", self.browse_directory)
        self.base_directory.pack(fill="x", padx=50, pady=5)

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
            'sequence_folder': self.sequence_folder.entry.get(),
            'reference_sequence': self.reference_sequence.entry.get(),
            'base_directory': self.base_directory.entry.get(),
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
            self.sequence_folder.entry.delete(0, tk.END)
            self.sequence_folder.entry.insert(0, config['sequence_folder'])
            self.reference_sequence.entry.delete(0, tk.END)
            self.reference_sequence.entry.insert(0, config['reference_sequence'])
            self.base_directory.entry.delete(0, tk.END)
            self.base_directory.entry.insert(0, config['base_directory'])

            # Inform the user that the configuration has been loaded
            tk.messagebox.showinfo("Load Configuration", "Configuration loaded successfully!")
        else:
            # Show an error message if the configuration file does not exist
            tk.messagebox.showerror("Load Error", "No configuration file found.")

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

    def browse_seq(self, entry):
        selected_folders = []

        while True:
            # Open a dialog for selecting a directory
            folder_path = filedialog.askdirectory()

            if folder_path:
                # List all files in the selected directory
                files = os.listdir(folder_path)

                # Filter out only files (ignore subdirectories)
                files = [f for f in files if os.path.isfile(os.path.join(folder_path, f))]

                # Check if the folder contains exactly 2 files
                if len(files) != 2:
                    # Show an error message if the folder does not contain exactly 2 files
                    tk.messagebox.showerror("Invalid Folder", "The selected folder must contain exactly 2 files.")
                else:
                    # Filter to include only .gz, .fastq, or .fq files
                    valid_files = [f for f in files if f.endswith(('.gz', '.fastq', '.fq'))]

                    # Check if both files are in valid formats
                    if len(valid_files) == 2:
                        # Add the selected folder path to the list
                        selected_folders.append(folder_path)
                    else:
                        # Show an error message if the files are not in valid formats
                        tk.messagebox.showerror("Invalid File Format",
                                                "The files in the selected folder must be in .gz, .fastq, or .fq format. Please ensure your files meet this requirement.")

            # Ask if the user wants to select another folder
            add_another = tk.messagebox.askyesno("Add Another Folder", "Do you want to add another folder?")
            if not add_another:
                break

        # Display the selected folders in the entry (comma-separated)
        if selected_folders:
            entry.delete(0, tk.END)
            entry.insert(0, ", ".join(selected_folders))
        else:
            # Show an info message if no folders were selected
            tk.messagebox.showinfo("No Folders Selected", "No folders were selected.")

    def validate_inputs(self):
        # Check if the sequence folder is selected
        if not self.sequence_folder.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please select the sequence folder.")
            return False

        # Check if reference sequence file entry is filled
        if not self.reference_sequence.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please select the reference genome file.")
            return False

        # Check if output directory entry is filled
        if not self.base_directory.entry.get().strip():
            tk.messagebox.showerror("Input Error", "Please select a directory.")
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

    def start_analysis(self):
        try:
            # Log the start of the analysis
            self.log_text.insert(tk.END, "Starting analysis...\n")
            self.update_idletasks()  # Update the GUI to reflect changes immediately

            # Retrieve the directories containing the forward and reverse reads
            sequence_folders = self.sequence_folder.entry.get().split(", ")

            for sequence_folder in sequence_folders:
                # List all files in the selected directory
                files = os.listdir(sequence_folder)

                # Filter out only files (ignore subdirectories)
                files = [f for f in files if os.path.isfile(os.path.join(sequence_folder, f))]

                # Sort files alphabetically to ensure the first is forward and the second is reverse
                files.sort()

                # Assign the first file as the forward read and the second as the reverse read
                forward_read_path = os.path.join(sequence_folder, files[0])
                reverse_read_path = os.path.join(sequence_folder, files[1])

                # Log the selected folder and file names
                self.log_text.insert(tk.END, f"\nRunning for {sequence_folder}\n")
                self.log_text.insert(tk.END, f"Forward read: {forward_read_path}\n")
                self.log_text.insert(tk.END, f"Reverse read: {reverse_read_path}\n")
                self.update_idletasks()

                # Retrieve the directory from the corresponding entry widget
                base_directory = self.base_directory.entry.get()

                # Create the output directory name based on the selected folder name
                output_directory_name = f"output_{os.path.basename(sequence_folder)}"

                # Create the full path for the output directory
                output_directory = os.path.join(base_directory, output_directory_name)

                # Create the output directory if it doesn't exist
                os.makedirs(output_directory, exist_ok=True)

                # Log the creation of the output directory
                self.log_text.insert(tk.END, f"\nCreated output directory: {output_directory}\n")
                self.update_idletasks()

                # Retrieve the reference sequence file path from the corresponding entry widget
                reference_sequence_path = self.reference_sequence.entry.get()

                # Access the ProgramPage and ConfigurationPage instances from the controller
                program_page = self.controller.pages["ProgramPage"]
                config_page = self.controller.pages["ConfigurationPage"]

                # Retrieve the paths for the FastQC, GATK, BWA, and ANNOVAR programs from the ProgramPage
                fastqc_path = program_page.fastqc_location.entry.get()
                gatk_path = program_page.gatk_location.entry.get()
                bwa_path = program_page.bwa_location.entry.get()
                annovar_path = program_page.annovar_location.entry.get()

                # Define a helper function to filter out optional arguments (those with None or empty values)
                def filter_optional_args(args):
                    return {k: v for k, v in args.items() if v}

                # Define a function to log messages to the log_text widget and keep it scrolled to the bottom
                def log_function(log):
                    self.log_text.insert(tk.END, log)
                    self.log_text.yview(tk.END)  # Scroll to the end of the log
                    self.update_idletasks()  # Update the GUI

                # Verify if the reference sequence file exists in the specified output directory
                if not self.is_file_in_directory(reference_sequence_path, base_directory):
                    self.log_text.insert(tk.END,
                                         f"Error: Reference sequence file {reference_sequence_path} is not in the selected directory {base_directory}.\n")
                    return  # Exit the function if the file is not found
                
                # FastQC Analysis
                try:
                    # Log the start of the FastQC process
                    self.log_text.insert(tk.END, "\nRunning FastQC...\n")
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
                    self.log_text.insert(tk.END, f"FastQC error: {str(e)}\n\n")

                # FastqToSam Conversion
                try:
                    # Log the start of the FastqToSam process
                    self.log_text.insert(tk.END, "\nRunning FastqToSam...\n")
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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the FastqToSam process
                    self.log_text.insert(tk.END, fastqtosam_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if FastqToSam fails
                    self.log_text.insert(tk.END, f"FastqToSam error: {str(e)}\n\n")

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
                    self.log_text.insert(tk.END, f"Samtools View error: {str(e)}\n\n")

                # MarkIlluminaAdapters Process
                try:
                    # Log the start of the MarkIlluminaAdapters process
                    self.log_text.insert(tk.END, "\nRunning MarkIlluminaAdapters...\n")
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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the MarkIlluminaAdapters process
                    self.log_text.insert(tk.END, markilluminaadapters_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if MarkIlluminaAdapters fails
                    self.log_text.insert(tk.END, f"MarkIlluminaAdapters error: {str(e)}\n\n")

                # SamToFastq Process
                try:
                    # Log the start of the SamToFastq process
                    self.log_text.insert(tk.END, "\nRunning SamToFastq...\n")
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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the SamToFastq process
                    self.log_text.insert(tk.END, samtofastq_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if SamToFastq fails
                    self.log_text.insert(tk.END, f"SamToFastq error: {str(e)}\n\n")

                # BWA Index Process
                try:
                    # Log the start of the BWA Index process
                    self.log_text.insert(tk.END, "\nChecking for existing BWA index files...\n")
                    self.update_idletasks()  # Update the GUI to show the log message immediately

                    # List of expected BWA index files based on the reference sequence filename
                    index_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa"]
                    index_files = [f"{reference_sequence_path}{ext}" for ext in index_extensions]

                    # Check if all expected BWA index files exist in the directory
                    if all(os.path.exists(index_file) for index_file in index_files):
                        self.log_text.insert(tk.END, "BWA index files already exist. Skipping BWA index process.\n\n")
                    else:
                        # Log the start of the BWA Index process
                        self.log_text.insert(tk.END, "Running BWA Index...\n")
                        self.update_idletasks()  # Update the GUI to show the log message immediately

                        # Retrieve the index algorithm type from the ConfigurationPage
                        index_algorithm = config_page.index_algorithm.get()

                        # Run the BWA Index command with the specified parameters
                        bwa_index_result = process_management.run_bwa_index(
                            bwa_path,  # Path to the BWA executable
                            reference_sequence_path,  # Path to the reference sequence file
                            index_algorithm,  # Algorithm used for indexing
                            log_function,  # Function to log messages during the process
                            output_directory
                        )

                        # Log the result of the BWA Index process
                        self.log_text.insert(tk.END, bwa_index_result + "\n")

                except Exception as e:
                    # Handle exceptions and log the error if BWA Index fails
                    self.log_text.insert(tk.END, f"BWA Index error: {str(e)}\n\n")

                # Update the GUI after the BWA Index process
                self.update_idletasks()

                # BWA Mem Process
                try:
                    # Log the start of the BWA Mem process
                    self.log_text.insert(tk.END, "Running BWA Mem...\n")
                    self.update_idletasks()  # Update the GUI to show the log message immediately

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the BWA Mem process
                    self.log_text.insert(tk.END, bwa_mem_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if BWA Mem fails
                    self.log_text.insert(tk.END, f"BWA Mem error: {str(e)}\n\n")

                # Update the GUI after the BWA Mem process
                self.update_idletasks()

                # Create Sequence Dictionary Process
                try:
                    # Log the start of the Create Sequence Dictionary process
                    self.log_text.insert(tk.END, "Creating Sequence Dictionary...\n")
                    self.update_idletasks()  # Update the GUI to show the log message immediately

                    # Define the path for the output sequence dictionary file
                    reference_dict = os.path.join(base_directory, os.path.basename(reference_sequence_path).replace('.fa', '.dict'))

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the Create Sequence Dictionary process
                    self.log_text.insert(tk.END, sequence_dict_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the Create Sequence Dictionary process fails
                    self.log_text.insert(tk.END, f"Create Sequence Dictionary error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the FAI index creation process
                    self.log_text.insert(tk.END, fai_index_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the FAI index creation fails
                    self.log_text.insert(tk.END, f"Samtools faidx error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the MergeBamAlignment process
                    self.log_text.insert(tk.END, mergebamalignment_result + "\n")
                except ValueError as ve:
                    # Handle specific exceptions like ValueError and log the error
                    self.log_text.insert(tk.END, f"MergeBamAlignment error: {str(ve)}\n")
                except Exception as e:
                    # Handle general exceptions and log the error
                    self.log_text.insert(tk.END, f"MergeBamAlignment error: {str(e)}\n\n")


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
                        log_function,
                        output_directory
                    )

                    # Log the result of the extraction process
                    self.log_text.insert(tk.END, extract_unmapped_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the extraction fails
                    self.log_text.insert(tk.END, f"Samtools view error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the SortSam process
                    self.log_text.insert(tk.END, sortsam_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the SortSam process fails
                    self.log_text.insert(tk.END, f"SortSam error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the SetNmMdAndUqTags process
                    self.log_text.insert(tk.END, setnmmduqtags_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the SetNmMdAndUqTags process fails
                    self.log_text.insert(tk.END, f"SetNmMdAndUqTags error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the MarkDuplicates process
                    self.log_text.insert(tk.END, markduplicates_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the MarkDuplicates process fails
                    self.log_text.insert(tk.END, f"MarkDuplicates error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the SortSam process
                    self.log_text.insert(tk.END, sortsam_marked_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the SortSam process fails
                    self.log_text.insert(tk.END, f"SortSam for Marked Duplicate BAM error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the BaseRecalibrator process
                    self.log_text.insert(tk.END, base_recal_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the BaseRecalibrator process fails
                    self.log_text.insert(tk.END, f"BaseRecalibrator error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the ApplyBQSR process
                    self.log_text.insert(tk.END, applybqsr_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the ApplyBQSR process fails
                    self.log_text.insert(tk.END, f"ApplyBQSR error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the stat process
                    self.log_text.insert(tk.END, stat_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the stat process fails
                    self.log_text.insert(tk.END, f"Samtools stat error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the HaplotypeCaller process
                    self.log_text.insert(tk.END, haplotypecaller_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the HaplotypeCaller process fails
                    self.log_text.insert(tk.END, f"HaplotypeCaller error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the SNV Extraction process
                    self.log_text.insert(tk.END, snv_extraction_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the SNV Extraction process fails
                    self.log_text.insert(tk.END, f"SNV Extraction error: {str(e)}\n\n")

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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the Indel Extraction process
                    self.log_text.insert(tk.END, indelextraction_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the Indel Extraction process fails
                    self.log_text.insert(tk.END, f"Indel Extraction error: {str(e)}\n\n")
                
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
                        log_function,  # Function to log messages during the process
                        output_directory
                    )

                    # Log the result of the CollectVariantCallingMetrics process
                    self.log_text.insert(tk.END, collectvariantcallingmetrics_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the CollectVariantCallingMetrics process fails
                    self.log_text.insert(tk.END, f"CollectVariantCallingMetrics error: {str(e)}\n\n")


                # Convert2ANNOVAR Process
                try:
                    # Log the start of the Convert2ANNOVAR process
                    self.log_text.insert(tk.END, "Running Convert2ANNOVAR ...\n")
                    self.update_idletasks()  # Update the GUI to show the log message immediately

                    # Define the input VCF files for conversion
                    input_vcfs = [
                        os.path.join(output_directory, config_page.haplotypecaller_output.get()),
                        # VCF from HaplotypeCaller
                        os.path.join(output_directory, config_page.snv_output.get()),  # VCF with extracted SNVs
                        os.path.join(output_directory, config_page.indel_output.get())  # VCF with extracted Indels
                    ]

                    # Define the output files for conversion
                    output_avinputs = [
                        input_vcfs[0].replace('.vcf', '.avinput'),  # Output AVInput for HaplotypeCaller VCF
                        input_vcfs[1].replace('.vcf', '.avinput'),  # Output AVInput for SNV VCF
                        input_vcfs[2].replace('.vcf', '.avinput')  # Output AVInput for Indel VCF
                    ]

                    # Optional arguments from Convert2ANNOVAR section in config_page
                    optional_args = filter_optional_args({
                        "includeinfo": config_page.convert2annovar_include_info.get(),
                        "withzyg": config_page.convert2annovar_with_zyg.get(),
                        "withfreq": config_page.convert2annovar_with_freq.get(),
                        "comment": config_page.convert2annovar_comment.get()
                    })

                    # Run Convert2ANNOVAR for each VCF file
                    for input_vcf, output_avinput in zip(input_vcfs, output_avinputs):
                        conversion_result = process_management.run_convert2annovar(
                            annovar_path,
                            input_vcf,
                            output_avinput,
                            optional_args,
                            log_function,
                            output_directory
                        )

                        # Log the result of the Convert2ANNOVAR process for the current VCF file
                        self.log_text.insert(tk.END, conversion_result + "\n")
                except Exception as e:
                    # Handle exceptions and log the error if the Convert2ANNOVAR process fails
                    self.log_text.insert(tk.END, f"Convert2ANNOVAR error: {str(e)}\n\n")

                # ANNOVAR Process
                try:
                    # Log the start of the ANNOVAR process
                    self.log_text.insert(tk.END, "Running ANNOVAR...\n")
                    self.update_idletasks()  # Update the GUI to show the log message immediately

                    # Define the AVInput files generated from Convert2ANNOVAR process
                    input_avinputs = [
                        os.path.join(output_directory,
                                     config_page.haplotypecaller_output.get().replace('.vcf', '.avinput')),
                        os.path.join(output_directory, config_page.snv_output.get().replace('.vcf', '.avinput')),
                        os.path.join(output_directory, config_page.indel_output.get().replace('.vcf', '.avinput'))
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

                    # Run ANNOVAR for each input AVInput file
                    for input_avinput in input_avinputs:
                        # Define the output prefix for the current AVInput file
                        output_prefix_avinput = os.path.join(output_directory,
                                                             f"{output_prefix}_{os.path.basename(input_avinput).split('.')[0]}")

                        # Run the ANNOVAR command with the provided inputs and optional arguments
                        annovar_result = process_management.run_annovar(
                            annovar_path,  # Path to the ANNOVAR executable
                            input_avinput,  # Input AVInput file
                            output_prefix_avinput,  # Prefix for the output files
                            buildver,  # Genome build version
                            db_dir,  # Database directory for ANNOVAR
                            optional_args_annovar,  # Dictionary of optional arguments
                            log_function,  # Function to log messages during the process
                            output_directory
                        )

                        # Log the result of the ANNOVAR process for the current AVInput file
                        self.log_text.insert(tk.END, annovar_result + "\n")

                except Exception as e:
                    # Handle exceptions and log the error if the ANNOVAR process fails
                    self.log_text.insert(tk.END, f"ANNOVAR error: {str(e)}\n\n")

                # Indicate the completion of the analysis
                self.log_text.insert(tk.END, "Analysis completed.\n\n")
                self.log_text.yview(tk.END)  # Scroll the log to the end to show the latest message

        except Exception as e:
            # Handle any unexpected exceptions that occur during the process
            self.log_text.insert(tk.END, f"Unexpected error: {str(e)}\n\n")
            self.update_idletasks()  # Update the GUI to show the error message immediately
