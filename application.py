import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd

from encoder import Encoder
from decoder import Decoder
from errorSimulator import ErrorSimulator


class Application(tk.Frame):
    def __init__(self, master):
        self.master = master
        self.master.title("Storing Digital data in DNA")
        self.master.geometry("500x300")

        self.encoder = Encoder()
        self.decoder = Decoder()
        self.errorSim = ErrorSimulator()

        self.setUpPages()

    def setUpPages(self):
        self.tabControl = ttk.Notebook(self.master)

        self.tab1 = ttk.Frame(self.tabControl)
        self.tab2 = ttk.Frame(self.tabControl)
        self.tab3 = ttk.Frame(self.tabControl)

        self.tabControl.add(self.tab1, text="Encoder")
        self.tabControl.add(self.tab2, text="Decoder")
        self.tabControl.add(self.tab3, text="Error simulation")
        self.tabControl.pack(expand=1, fill="both")

        self.encoderPage()
        self.decoderPage()
        self.errorPage()

    def encoderPage(self):
        self.BasePage("Encoder", self.encode, self.tab1)

    def decoderPage(self):
        self.BasePage("Decoder", self.decode, self.tab2)

    def errorPage(self):
        self.BasePage("Error simulation", self.induceErrors, self.tab3)

    def BasePage(self, pageName, opperation, tab):
        self.inputFilePath = None
        self.outputFilePath = None

        label = tk.Label(tab, text=pageName)
        label.pack()

        label = tk.Label(tab, text="input file")
        label.pack()

        container1 = tk.Label(tab)
        container1.pack()

        self.inputFilePathLabel = tk.Label(
            container1, text="...", width=50, height=1, borderwidth=1, relief="solid"
        )
        self.inputFilePathLabel.pack(side=tk.LEFT, padx=20)

        greet_button = tk.Button(
            container1, text="Select File", command=self.getInputFilePath, width=10
        )
        greet_button.pack(side=tk.LEFT, padx=0)

        label = tk.Label(tab, text="output file")
        label.pack()

        container2 = tk.Label(tab)
        container2.pack()

        self.outputFilePathLabel = tk.Label(
            container2, text="...", width=50, height=1, borderwidth=1, relief="solid"
        )
        self.outputFilePathLabel.pack(side=tk.LEFT, padx=20)

        greet_button = tk.Button(
            container2, text="Select File", command=self.getOutputFilePath, width=10
        )
        greet_button.pack(side=tk.LEFT, padx=0)

        encode_button = tk.Button(tab, text=pageName, command=opperation, width=20)
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
        if self.inputFilePath and self.outputFilePath:
            self.encoder.encode(self.inputFilePath, self.outputFilePath)

    def decode(self):
        if self.inputFilePath and self.outputFilePath:
            self.decoder.decode(self.inputFilePath, self.outputFilePath)

    def induceErrors(self):
        if self.inputFilePath and self.outputFilePath:
            self.errorSim.induceErrors(self.inputFilePath, self.outputFilePath)
