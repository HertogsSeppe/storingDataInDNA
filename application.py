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

        self.inputFilePath = None
        self.outputFilePath = None

        self.setUpPages()

    def setUpPages(self):
        self.tabControl = ttk.Notebook(self.master)

        self.tab1 = ttk.Frame(self.tabControl)
        self.tab2 = ttk.Frame(self.tabControl)
        self.tab3 = ttk.Frame(self.tabControl)

        self.tabControl.add(self.tab1, text="Home")
        self.tabControl.add(self.tab2, text="Settings")
        self.tabControl.pack(expand=1, fill="both")
        # self.tabControl.bind("<<NotebookTabChanged>>", self.OnTabChange)

        self.homePage()

    def homePage(self):
        #  --  Selection Box  --  #
        self.operation = tk.StringVar(self.tab1)
        self.operation.set("Decoder")  # default value

        w = tk.OptionMenu(
            self.tab1, self.operation, "Encoder", "Decoder", "Error Simulator"
        )
        w.pack()

        #  --  Input file field  --  #
        tk.Label(self.tab1, text="input file").pack()

        container1 = tk.Label(self.tab1)
        container1.pack()

        self.inputFilePathLabel = tk.Label(
            container1, text="...", width=50, height=1, borderwidth=1, relief="solid"
        )
        self.inputFilePathLabel.pack(side=tk.LEFT, padx=20)

        tk.Button(
            container1, text="Select File", command=self.getInputFilePath, width=10
        ).pack(side=tk.LEFT, padx=0)

        #  --  Output file field  --  #
        tk.Label(self.tab1, text="output file").pack()

        container2 = tk.Label(self.tab1)
        container2.pack()

        self.outputFilePathLabel = tk.Label(
            container2, text="...", width=50, height=1, borderwidth=1, relief="solid"
        )
        self.outputFilePathLabel.pack(side=tk.LEFT, padx=20)

        tk.Button(
            container2, text="Select File", command=self.getOutputFilePath, width=10
        ).pack(side=tk.LEFT, padx=0)

        #  --  Submit button  --  #
        encode_button = tk.Button(
            self.tab1,
            text="Start",
            command=self.submit,
            width=20,
        )
        encode_button.pack()

    # def OnTabChange(self, *args):
    #     self.tabIndex = int(self.tabControl.index(self.tabControl.select()))
    #     print(self.tabIndex)

    def getInputFilePath(self):
        self.inputFilePath = fd.askopenfilename()
        self.inputFilePathLabel.configure(text=self.inputFilePath)

    def getOutputFilePath(self):
        self.outputFilePath = fd.askopenfilename()
        self.outputFilePathLabel.configure(text=self.outputFilePath)

    def submit(self):
        if not self.inputFilePath or not self.outputFilePath:
            print("No file paths selected")
            return
        operation = self.operation.get()
        if operation == "Encoder":
            self.encoder.encode(self.inputFilePath, self.outputFilePath)
        if operation == "Decoder":
            self.decoder.decode(self.inputFilePath, self.outputFilePath)
        if operation == "Error Simulator":
            self.errorSim.induceErrors(self.inputFilePath, self.outputFilePath)
