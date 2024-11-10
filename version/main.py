import os
import json
import threading
import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext
from main_page import MainPage
from configuration_page import ConfigurationPage
from program_page import ProgramPage
from manual_pages import AboutPage, InstallationPage, GettingStartedPage, WorkflowPage, IssuesPage

# Define the main application class, inheriting from Tkinter's Tk class
class SequenceAssemblyApp(tk.Tk):
    def __init__(self):
        super().__init__()

        # Set the pages, window title and size
        self.pages = {}
        self.title("Automated Sequence Assembly Software")
        self.geometry("1200x800")
        self.iconphoto(False, tk.PhotoImage(file='ASAWES.png'))

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
        help_menu.add_command(label="Workflow", command=lambda: self.show_page("WorkflowPage"))
        help_menu.add_command(label="Common Issues", command=lambda: self.show_page("IssuesPage"))

        # Container to hold all the pages (frames)
        container = tk.Frame(self)
        container.pack(fill="both", expand=True)

        # Loop through the pages and initialize each one
        for page in (
        MainPage, ConfigurationPage, ProgramPage, AboutPage, InstallationPage, GettingStartedPage, WorkflowPage, IssuesPage):

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

    # Main Application Execution
if __name__ == "__main__":
    app = SequenceAssemblyApp()
    app.mainloop()

