import subprocess
import os

def run_fastqc(fastqc_path, input_files, output_directory, log_function):
    # Set up logging directory and files
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'fastqc_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'fastqc_error_file.txt')

    command = [fastqc_path, '-o', output_directory] + input_files

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        # Log the command being run
        log_function(f"Running command: {' '.join(command)}\n")
        # Running the command and capturing the output
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)

        # Wait for the process to complete and get the exit code
        return_code = process.wait()

        # Check the process return code to determine if there were errors
        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "FastQC completed successfully.\n"

def run_fastqtosam(gatk_path, forward_read, reverse_read, output_bam, sample_name, platform, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'fastqtosam_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'fastqtosam_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "FastqToSam",
        "--FASTQ", forward_read,
        "--OUTPUT", output_bam,
        "--SAMPLE_NAME", sample_name,
        "--PLATFORM", platform
    ]
    if reverse_read:
        command += ["--FASTQ2", reverse_read]
    for key, value in optional_args.items():
        if value is not None:
            command += [f"--{key}", str(value)]

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        log_function(f"Running command: {' '.join(command)}\n")
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "FastqToSam completed successfully.\n"

def run_samtools_view(bam_file, output_directory, log_function):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'samtools_view_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'samtools_view_error_file.txt')

    output_header = os.path.join(output_directory, "header.txt")
    command = ["samtools", "view", "-H", bam_file]

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        log_function(f"Running command: {' '.join(command)}\n")
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    with open(output_header, "w") as f:
        with open(log_file_path) as log_f:
            f.write(log_f.read())

    return "Samtools View completed successfully.\n"


def run_markilluminaadapters(gatk_path, input_bam, metrics_file, output_bam, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'markilluminaadapters_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'markilluminaadapters_error_file.txt')

    command = ["java", "-jar", gatk_path, "MarkIlluminaAdapters", "--INPUT", input_bam, "--OUTPUT", output_bam,
               "--METRICS", metrics_file]

    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key}", str(value)])

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        log_function(f"Running command: {' '.join(command)}\n")
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "MarkIlluminaAdapters completed successfully.\n"

def run_samtofastq(gatk_path, input_bam, fastq_output, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'samtofastq_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'samtofastq_error_file.txt')

    command = ["java", "-jar", gatk_path, "SamToFastq", "--INPUT", input_bam, "--FASTQ", fastq_output]

    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key}", str(value)])

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        log_function(f"Running command: {' '.join(command)}\n")
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "SamToFastq completed successfully.\n"

def run_bwa_index(bwa_path, reference_sequence, algorithm, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'bwa_index_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'bwa_index_error_file.txt')

    command = [bwa_path, "index", "-a", algorithm, reference_sequence]

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "BWA Indexing completed successfully.\n"


def run_bwa_mem(bwa_path, reference_sequence, fastq_output, sam_output, threads, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    error_file_path = os.path.join(logging_dir, 'bwa_mem_error_file.txt')

    command = [bwa_path, "mem", "-t", str(threads), "-p", "-M", reference_sequence, fastq_output]

    log_function(f"Running command: {' '.join(command)}\n")


    with open(sam_output, "w") as sam_file, open(error_file_path, "w") as err_f:
        # Create the process
        process = subprocess.Popen(command, stdout=sam_file, stderr=err_f, text=True)
        return_code = process.wait()

    if return_code:
        raise subprocess.CalledProcessError(return_code, command, output=f"Check error log at: {error_file_path}")


    return "BWA Mem alignment completed successfully.\n"


def run_create_sequence_dictionary(gatk_path, reference_sequence, reference_dict, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'create_sequence_dictionary_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'create_sequence_dictionary_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "CreateSequenceDictionary",
        "-R", reference_sequence,
        "-O", reference_dict
    ]

    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "CreateSequenceDictionary completed successfully.\n"

def run_samtools_faidx(reference_sequence, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'samtools_faidx_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'samtools_faidx_error_file.txt')

    command = ["samtools", "faidx", reference_sequence]

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "Samtools faidx completed successfully.\n"

def run_mergebamalignment(gatk_path, unaligned_bam, aligned_bam, reference_sequence, output_bam, optional_args,
                          log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'mergebamalignment_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'mergebamalignment_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "MergeBamAlignment",
        "-UNMAPPED_BAM", unaligned_bam,
        "-ALIGNED_BAM", aligned_bam,
        "-R", reference_sequence,
        "-O", output_bam
    ]

    # Handling list-type arguments separately
    list_args = ["ATTRIBUTES_TO_RETAIN", "ATTRIBUTES_TO_REVERSE", "ATTRIBUTES_TO_REVERSE_COMPLEMENT",
                 "EXPECTED_ORIENTATIONS", "READ1_ALIGNED_BAM", "READ2_ALIGNED_BAM"]

    for key, value in optional_args.items():
        if value is True:
            command.extend([f"--{key}"])
        elif value:
            if key in list_args and isinstance(value, str):
                values = value.split(',')
                for val in values:
                    command.extend([f"--{key}", val])
            else:
                command.extend([f"--{key}", str(value)])

    # Check for mandatory group arguments
    program_group_required = ['PROGRAM_RECORD_ID', 'PROGRAM_GROUP_NAME', 'PROGRAM_GROUP_VERSION',
                              'PROGRAM_GROUP_COMMAND_LINE']
    if all(optional_args.get(arg) for arg in program_group_required):
        for arg in program_group_required:
            command.extend([f"--{arg}", str(optional_args[arg])])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "MergeBamAlignment completed successfully.\n"

def run_samtools_view_unmapped(merge_bam_output, unmapped_output, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'samtools_view_unmapped_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'samtools_view_unmapped_error_file.txt')

    command = ["samtools", "view", "-f", "4", merge_bam_output]

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return f"Unmapped reads extracted to {unmapped_output}\n"

def run_sortsam(gatk_path, merge_bam_output, sorted_output, sort_order, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'sortsam_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'sortsam_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "SortSam",
        "--INPUT", merge_bam_output,
        "--OUTPUT", sorted_output,
        "--SORT_ORDER", sort_order
    ]

    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "SortSam completed successfully.\n"

def run_setnmmduqtags(gatk_path, sorted_output, output_bam_fixed, reference_sequence, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'setnmmduqtags_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'setnmmduqtags_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "SetNmMdAndUqTags",
        "--INPUT", sorted_output,
        "--OUTPUT", output_bam_fixed,
        "--REFERENCE_SEQUENCE", reference_sequence
    ]

    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "SetNmMdAndUqTags completed successfully.\n"

def run_markduplicates(gatk_path, output_bam_fixed, marked_output_bam, metrics_file, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'markduplicates_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'markduplicates_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "MarkDuplicates",
        "--INPUT", output_bam_fixed,
        "--OUTPUT", marked_output_bam,
        "--METRICS_FILE", metrics_file
    ]

    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "MarkDuplicates completed successfully.\n"

def run_sortsam_marked(gatk_path, input_bam, sorted_output, sort_order, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'sortsam_marked_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'sortsam_marked_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "SortSam",
        "--INPUT", input_bam,
        "--OUTPUT", sorted_output,
        "--SORT_ORDER", sort_order
    ]

    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "SortSam for marked duplicates completed successfully.\n"


def run_baserecalibrator(gatk_path, input_bam, reference_sequence, output_table, known_sites, optional_args, log_function, output_directory):
    # Ensure logging directory exists
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)

    # Define log and error file paths
    log_file_path = os.path.join(logging_dir, 'baserecalibrator_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'baserecalibrator_error_file.txt')

    # Construct the base command
    command = [
        "java", "-jar", gatk_path, "BaseRecalibrator",
        "--input", input_bam.strip(),
        "--reference", reference_sequence.strip(),
        "--output", output_table.strip()
    ]

    # Add known sites, ensuring no leading or trailing whitespace
    if known_sites:
        known_sites_list = [site.strip() for site in known_sites.split(',')]
        for site in known_sites_list:
            if site:  # Ensure the site is not empty
                command.extend(["--known-sites", site])

    # Add optional arguments, filtering out empty or None values
    for key, value in optional_args.items():
        if value not in [None, ""]:
            command.extend([f"--{key}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    # Execute the command, capturing output and errors
    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        # Check for errors
        if return_code:
            with open(error_file_path, "r") as error_f:
                error_messages = error_f.read()
            raise Exception(f"Error running BaseRecalibrator: {error_messages}")

    return "BaseRecalibrator completed successfully.\n"

def run_applybqsr(gatk_path, input_bam, output_bam, reference_sequence, bqsr_recal_file, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'applybqsr_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'applybqsr_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "ApplyBQSR",
        "--input", input_bam,
        "--output", output_bam,
        "--reference", reference_sequence,
        "--bqsr-recal-file", bqsr_recal_file
    ]

    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "ApplyBQSR completed successfully.\n"

def run_samtools_flagstat(input_bam, output_flagstat, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'samtools_flagstat_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'samtools_flagstat_error_file.txt')

    command = ["samtools", "flagstat", input_bam]

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        log_function(f"Running command: {' '.join(command)}\n")
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    with open(output_flagstat, 'w') as f:
        with open(log_file_path) as log_f:
            f.write(log_f.read())

    return "Samtools flagstat completed successfully.\n"

def run_samtools_stats(input_bam, output_stats, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'samtools_stats_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'samtools_stats_error_file.txt')

    command = ["samtools", "stats", input_bam]

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        log_function(f"Running command: {' '.join(command)}\n")
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err_f, text=True)

        # Handle the output with grep in another subprocess
        grep = subprocess.Popen(["grep", "^SN"], stdin=process.stdout, stdout=open(output_stats, 'w'), text=True)
        process.stdout.close()  # Correctly close stdout after piping it
        grep_stdout, grep_stderr = grep.communicate()  # Ensure we collect any output and errors

        return_code = process.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, command, output=grep_stderr)

    return "Samtools stats completed successfully.\n"

def run_haplotypecaller(gatk_path, input_bam, output_vcf, reference_sequence, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'haplotypecaller_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'haplotypecaller_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "HaplotypeCaller",
        "--input", input_bam,
        "--output", output_vcf,
        "--reference", reference_sequence
    ]

    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key.replace('_', '-')}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "HaplotypeCaller completed successfully.\n"

def run_snvextraction(gatk_path, input_vcf, output_vcf, reference, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'snvextraction_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'snvextraction_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "SelectVariants",
        "--V", input_vcf,
        "--O", output_vcf,
        "--R", reference
    ]

    for key, value in optional_args.items():
        if value is not None and value != "":  # Only include non-empty arguments
            command.extend([f"--{key.replace('_', '-')}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "SNV Extraction completed successfully.\n"

def run_indelextraction(gatk_path, input_vcf, output_vcf, reference, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'indelextraction_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'indelextraction_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "SelectVariants",
        "--variant", input_vcf,
        "--output", output_vcf,
        "--R", reference
    ]

    for key, value in optional_args.items():
        if value is not None and value != "":  # Only include non-empty arguments
            command.extend([f"--{key.replace('_', '-')}", str(value)])

    # Log the command being run
    log_function(f"Running command: {' '.join(command)}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "Indel extraction completed successfully.\n"

def run_collectvariantcallingmetrics(gatk_path, input_vcf, output_path, dbsnp_file, optional_args, log_function, output_directory):
    # Set up logging directory and files
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'collectvariantcallingmetrics_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'collectvariantcallingmetrics_error_file.txt')

    command = [
        "java", "-jar", gatk_path, "CollectVariantCallingMetrics",
        "--INPUT", input_vcf,
        "--OUTPUT", output_path,
        "--DBSNP", dbsnp_file
    ]

    # Adding optional arguments to the command
    for key, value in optional_args.items():
        if value:  # Only include non-empty arguments
            command.extend([f"--{key}", str(value)])

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        # Running the command and capturing the output
        process = subprocess.Popen(command, stdout=log_f, stderr=err_f, text=True)

        # Log the command being run
        log_function(f"Running command: {' '.join(command)}\n")

        # Wait for the process to complete and get the exit code
        return_code = process.wait()

        # Check the process return code to determine if there were errors
        if return_code:
            raise subprocess.CalledProcessError(return_code, command)

    return "CollectVariantCallingMetrics completed successfully.\n"


def run_convert2annovar(annovar_path, input_vcf, output_avinput, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'convert2annovar_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'convert2annovar_error_file.txt')

    # Build the basic command to execute Convert2ANNOVAR
    command = [
        "perl", os.path.join(annovar_path, "convert2annovar.pl"),
        "-format", "vcf4",
        f'"{input_vcf}"',
        "-outfile", f'"{output_avinput}"'
    ]

    # Adding optional arguments
    if optional_args.get("includeinfo"):
        command.append("-includeinfo")
    if optional_args.get("withzyg"):
        command.append("-withzyg")
    if optional_args.get("withfreq"):
        command.append("-withfreq")
    if optional_args.get("comment"):
        command.append("-comment")

    # Convert command list to string with proper quoting
    command_str = ' '.join(command)

    # Log the command being run
    log_function(f"Running command: {command_str}\n")

    # Run the Convert2ANNOVAR process
    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command_str, shell=True, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command_str)

    return "Convert2ANNOVAR completed successfully.\n"


def run_annovar(annovar_path, input_vcf, output_prefix, buildver, db_dir, optional_args, log_function, output_directory):
    logging_dir = os.path.join(output_directory, 'log')
    os.makedirs(logging_dir, exist_ok=True)
    log_file_path = os.path.join(logging_dir, 'annovar_log_file.txt')
    error_file_path = os.path.join(logging_dir, 'annovar_error_file.txt')

    command = [
        "perl", os.path.join(annovar_path, "table_annovar.pl"),
        f'"{input_vcf}"',  # Input VCF path
        f'"{db_dir}"',  # Database directory path
        "-buildver", f'"{buildver}"',
        "-out", f'"{output_prefix}"',
        "-csvout"
    ]

    # Adding optional arguments
    if optional_args.get("protocol"):
        command.extend(["-protocol", optional_args["protocol"]])
    if optional_args.get("operation"):
        command.extend(["-operation", optional_args["operation"]])
    if optional_args.get("nastring"):
        command.extend(["-nastring", optional_args["nastring"]])
    if optional_args.get("remove"):
        command.append("-remove")
    if optional_args.get("nopolish"):
        command.append("-nopolish")

    # Convert command list to string with proper quoting
    command_str = ' '.join(command)

    # Log the command being run
    log_function(f"Running command: {command_str}\n")

    with open(log_file_path, "w") as log_f, open(error_file_path, "w") as err_f:
        process = subprocess.Popen(command_str, shell=True, stdout=log_f, stderr=err_f, text=True)
        return_code = process.wait()

        if return_code:
            raise subprocess.CalledProcessError(return_code, command_str)

    return "ANNOVAR completed successfully.\n"
