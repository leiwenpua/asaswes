import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import json

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
        messagebox.showinfo("Save Configuration", "Configuration saved successfully!")

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
            messagebox.showinfo("Load Configuration", "Configuration loaded successfully!")
        else:
            # Show an error message if the configuration file does not exist
            messagebox.showerror("Load Error", "No configuration file found.")
