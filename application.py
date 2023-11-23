import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd

from encoder import Encoder


class Application(tk.Frame):
    def __init__(self, master):
        self.master = master
        self.master.title("Storing Digital data in DNA")
        self.master.geometry("500x300")

        self.encoder = Encoder()

        self.setUpPages()

    def setUpPages(self):
        self.tabControl = ttk.Notebook(self.master)

        self.tab1 = ttk.Frame(self.tabControl)
        self.tab2 = ttk.Frame(self.tabControl)

        self.tabControl.add(self.tab1, text="Tab 1")
        self.tabControl.add(self.tab2, text="Tab 2")
        self.tabControl.pack(expand=1, fill="both")

        self.encoderPage()

    def encoderPage(self):
        label = tk.Label(self.tab1, text="Encoder")
        label.pack()

        label = tk.Label(self.tab1, text="input file")
        label.pack()

        container1 = tk.Label(self.tab1)
        container1.pack()

        self.inputFilePathLabel = tk.Label(
            container1, text="...", width=50, height=1, borderwidth=1, relief="solid"
        )
        self.inputFilePathLabel.pack(side=tk.LEFT, padx=20)

        greet_button = tk.Button(
            container1, text="Select File", command=self.getInputFilePath, width=10
        )
        greet_button.pack(side=tk.LEFT, padx=0)

        label = tk.Label(self.tab1, text="output file")
        label.pack()

        container2 = tk.Label(self.tab1)
        container2.pack()

        self.outputFilePathLabel = tk.Label(
            container2, text="...", width=50, height=1, borderwidth=1, relief="solid"
        )
        self.outputFilePathLabel.pack(side=tk.LEFT, padx=20)

        greet_button = tk.Button(
            container2, text="Select File", command=self.getOutputFilePath, width=10
        )
        greet_button.pack(side=tk.LEFT, padx=0)

        encode_button = tk.Button(
            self.tab1, text="Select File", command=self.encode, width=10
        )
        encode_button.pack()

    def greet(self):
        print("Greet")

    def getInputFilePath(self):
        self.inputFilePath = fd.askopenfilename()
        self.inputFilePathLabel.configure(text=self.inputFilePath)

    def getOutputFilePath(self):
        self.outputFilePath = fd.askopenfilename()
        self.outputFilePathLabel.configure(text=self.outputFilePath)

    def encode(self):
        if self.inputFilePath & self.outputFilePath:
            self.encoder.encode(self.inputFilePath, self.outputFilePath)
