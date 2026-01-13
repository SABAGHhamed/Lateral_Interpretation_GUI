from .tsp import *
import cigsegy

import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog   # new
import shutil  # new
import tkinter.font as tkfont
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_svg
from PIL import Image
from tkinter import filedialog, simpledialog, messagebox
#from tkinterdnd2 import TkinterDnD, DND_FILES
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RectangleSelector
from matplotlib.figure import Figure
import numpy as np
import hashlib
from skimage.filters.rank import entropy
from skimage.morphology import disk
from skimage.util import img_as_ubyte
import sys
 

   

class SeismicApp:
    def __init__(self, root):
        self.root = root

        # =====================
        # Internal state
        # =====================
        self.seismic_profile_full_relief = None
        self.data = None
        self.inline = None
        self.x_start = None
        self.y_start = None
        self.crossline_val = None
        self.twt_spacing = None
        self.initial_inline = None
        self.selection_finalized = False
        self.is_image = False  # for 2D image imports not 3D SEG=Y files
        # ============================================================
        # Canvas placeholder
        # ============================================================
        self.canvas_widget = None
        self.seismic_profile = None        

        # =====================
        # Window config
        # =====================
        self.root.title("Lateral Interpretation GUI")
        center_window(self.root, 700, 700)

        # ============================================================
        # TOP STATUS LABEL
        # ============================================================
        self.label = ttk.Label(
            root,
            text="Enter the file name (in the main directory) or click “Browse”"
        )
        self.label.pack(padx=10, pady=10)

        # ============================================================
        # ENTRY + BROWSE
        # ============================================================
        entry_frame = ttk.Frame(root)
        entry_frame.pack(pady=5)

        self.entry = ttk.Entry(entry_frame, width=60)
        self.entry.insert(0, "Seismic_data.sgy")
        self.entry.pack(side=tk.LEFT)

        ttk.Button(
            entry_frame,
            text="Browse",
            command=self.browse_sgy_file
        ).pack(side=tk.LEFT, padx=5)

        # ============================================================
        # RELIEF CONTROLS + LOAD BUTTON
        # ============================================================
        relief_frame = ttk.Frame(root)
        relief_frame.pack(pady=5)

        self.load_button = ttk.Button(
            relief_frame,
            text="Load Selected File",
            command=self.load_from_entry
        )
        self.load_button.grid(row=0, column=0, padx=5)

        self.relief_var = tk.IntVar(value=99)
        self.relief_spinbox = ttk.Spinbox(
            relief_frame,
            from_=51,
            to=100,
            textvariable=self.relief_var,
            width=5
        )
        self.relief_spinbox.grid(row=0, column=1, padx=5)

        self.update_relief_button = ttk.Button(
            relief_frame,
            text="Update Relief",
            command=self.update_seismic_relief
        )
        self.update_relief_button.grid(row=0, column=2, padx=5)

        # ============================================================
        # PARAMETERS
        # ============================================================
        params_frame = ttk.Frame(root)
        params_frame.pack(pady=5)

        ttk.Label(params_frame, text="Initial inline/crossline:")\
            .grid(row=0, column=0, sticky=tk.W)

        self.l_threshold_var = tk.DoubleVar(value=0)
        self.l_threshold_entry = ttk.Entry(
            params_frame, textvariable=self.l_threshold_var, width=8
        )
        self.l_threshold_entry.grid(row=0, column=1, padx=5)

        ttk.Label(params_frame, text="Time Sample Interval:")\
            .grid(row=0, column=2, sticky=tk.W)

        self.rescaling_constant_var = tk.DoubleVar(value=4.0)
        self.rescaling_constant_entry = ttk.Entry(
            params_frame, textvariable=self.rescaling_constant_var, width=6
        )
        self.rescaling_constant_entry.grid(row=0, column=3, padx=5)

        ttk.Label(params_frame, text="x_start:")\
            .grid(row=0, column=4, sticky=tk.W)

        self.x_start_var = tk.DoubleVar(value=0.0)
        self.x_start_entry = ttk.Entry(
            params_frame, textvariable=self.x_start_var, width=6
        )
        self.x_start_entry.grid(row=0, column=5, padx=5)

        ttk.Label(params_frame, text="y_start:")\
            .grid(row=0, column=6, sticky=tk.W)

        self.y_start_var = tk.DoubleVar(value=0.0)
        self.y_start_entry = ttk.Entry(
            params_frame, textvariable=self.y_start_var, width=6
        )
        self.y_start_entry.grid(row=0, column=7, padx=5)

        # ============================================================
        # MAIN BUTTONS
        # ============================================================
        button_frame = ttk.Frame(root)
        button_frame.pack(pady=8)

        self.manual_interp_var = tk.IntVar(value=1)

        ttk.Label(button_frame, text="Manual Interpretation:")\
            .grid(row=0, column=1, padx=(10, 2))

        ttk.Radiobutton(
            button_frame, text="Yes",
            variable=self.manual_interp_var, value=0
        ).grid(row=0, column=2, padx=(0, 1))

        ttk.Radiobutton(
            button_frame, text="No",
            variable=self.manual_interp_var, value=1
        ).grid(row=0, column=3, padx=(1, 0))

        self.crossline_var = tk.IntVar(value=0)

        ttk.Radiobutton(
            button_frame, text="inline",
            variable=self.crossline_var, value=0
        ).grid(row=1, column=0, padx=2, pady=(0, 2))

        ttk.Radiobutton(
            button_frame, text="crossline",
            variable=self.crossline_var, value=1
        ).grid(row=1, column=1, padx=2, pady=(0, 2))

        self.ok_button = ttk.Button(
            button_frame,
            text="Proceed to Lateral Interpretation",
            command=self.finalize_selection,
            state="disabled"
        )
        self.ok_button.grid(row=0, column=4, padx=5)

        self.more_button = ttk.Button(
            button_frame,
            text="More",
            command=self.show_more_menu
        )
        self.more_button.grid(row=0, column=5, padx=5)

        self.help_button = ttk.Button(
            button_frame,
            text="Help",
            command=self.open_help_window
        )
        self.help_button.grid(row=0, column=6, padx=5)

        # ============================================================
        # DEVELOPER OPTIONS
        # ============================================================

        self.dev_container = tk.Frame(self.root)
        
        self.dev_container.pack(fill="x")

        self.dev_frame = tk.Frame(self.dev_container, borderwidth=1, relief="groove")
        
        
        self.dev_options_visible = False

        #self.dev_frame = ttk.Frame(root)

        self.dev_toggle = ttk.Button(
            root,
            text="Lateral Interpretation Steps ▼",
            command=self.toggle_dev_options
        )
        self.dev_toggle.pack(anchor="e", padx=10)

        ttk.Button(
            self.dev_frame, text="Contour Detection",
            command=self.skip_to_contours_processing
        ).grid(row=0, column=0, padx=4, pady=3)

        ttk.Button(
            self.dev_frame, text="TSP",
            command=self.skip_to_tsp_processing
        ).grid(row=0, column=1, padx=4, pady=3)

        ttk.Button(
            self.dev_frame, text="Onlap/Downlap Connection",
            command=self.skip_to_onlap_processing
        ).grid(row=0, column=2, padx=4, pady=3)

        ttk.Button(
            self.dev_frame, text="Vertical Smoothing",
            command=self.skip_to_smooth_processing
        ).grid(row=0, column=3, padx=4, pady=3)

        ttk.Button(
            self.dev_frame, text="Feature Extraction",
            command=self.skip_to_feature_processing
        ).grid(row=0, column=4, padx=4, pady=3)

        
    # === Help menu ===    
    def open_help_window(self):
        # Prevent multiple help windows
        if hasattr(self, "help_window") and self.help_window.winfo_exists():
            self.help_window.lift()
            return

        # 1️⃣ Create window FIRST
        self.help_window = tk.Toplevel(self.root)
        style = ttk.Style(self.help_window)
#        style.theme_use("default")
        self.help_window.configure(background="white")        
        #center_window(self.help_window)                
        self.help_window.resizable(True, True)
        self.help_window.title("Help")

        # 2️⃣ Create notebook
        notebook = ttk.Notebook(self.help_window)
        notebook.pack(fill=tk.BOTH, expand=True)
        
        # Style
        style = ttk.Style(self.help_window)
#        style.theme_use("default")   # important for background control

        style.configure(
            "TNotebook",
            background="white",
            borderwidth=0
        )

        style.configure(
            "TNotebook.Tab",
            background="white",
            padding=(10, 5)
        )

        style.map(
            "TNotebook.Tab",
            background=[
                ("selected", "white"),
                ("active", "white")
            ]
        )

        # 3️⃣ Create tabs
        params_tab = tk.Frame(notebook, bg="white")
        tutorial_tab = tk.Frame(notebook, bg="white")
        manual_tab = tk.Frame(notebook, bg="white")
        plotting_tab = tk.Frame(notebook, bg="white")
        credit_tab = tk.Frame(notebook, bg="white")        


        notebook.add(tutorial_tab, text="Quick Tutorial")        
        notebook.add(params_tab, text="Parameters")
        notebook.add(manual_tab, text="Manual Interpretation")                
        notebook.add(plotting_tab, text="Plotting")                        
        notebook.add(credit_tab, text="Credit")                                
        
        
        
        # below till tutorial tab is needed for tutorial tab and parameters tab
        def build_paged_section(
            parent_tab,
            container,
            items,
            items_per_page,
            index_attr_name,
            bg="white"
        ):
            # ---------- build pages ----------
            def build_page(parent, items, bg=bg):
                page = tk.Frame(parent, bg=bg)

                for name, desc in items:
                    section = tk.Frame(page, bg=bg)
                    section.pack(fill=tk.X, pady=6)

                    has_title = bool(name.strip())

                    if has_title:
                        tk.Label(
                            section,
                            text=name,
                            font=("Arial", 12, "bold"),
                            anchor="w",
                            bg=bg
                        ).pack(anchor="w")

                    tk.Label(
                        section,
                        text=desc,
                        font=("Arial", 10),
                        wraplength=700,
                        justify=tk.LEFT,
                        anchor="w",
                        bg=bg,
                        padx=12,
                        pady=(10 if has_title else 0)   # ← key line
                    ).pack(anchor="w", padx=(12, 0))


                return page

            num_pages = max(1, math.ceil(len(items) / items_per_page))

            pages_data = [
                items[i * items_per_page:(i + 1) * items_per_page]
                for i in range(num_pages)
            ]

            pages = [build_page(container, data) for data in pages_data]

            setattr(self, index_attr_name, 0)
            pages[0].pack(fill=tk.BOTH, expand=True)

            # ---------- navigation ----------
            nav = tk.Frame(parent_tab, bg=bg)
            nav.pack(pady=8)

            page_label = tk.Label(
                nav,
                text=f"Page 1 / {len(pages)}",
                bg=bg
            )
            page_label.pack(side=tk.LEFT, padx=10)

            prev_btn = tk.Button(nav, text="◀ Previous")
            prev_btn.pack(side=tk.LEFT, padx=5)

            next_btn = tk.Button(nav, text="Next ▶")
            next_btn.pack(side=tk.LEFT, padx=5)

            # ---------- callbacks (BOUND SAFELY) ----------
            def update_nav_buttons(
                pages=pages,
                prev_btn=prev_btn,
                next_btn=next_btn
            ):
                idx = getattr(self, index_attr_name)
                prev_btn.config(state=tk.NORMAL if idx > 0 else tk.DISABLED)
                next_btn.config(state=tk.NORMAL if idx < len(pages) - 1 else tk.DISABLED)

            def show_page(
                i,
                pages=pages,
                page_label=page_label
            ):
                idx = getattr(self, index_attr_name)
                if 0 <= i < len(pages):
                    pages[idx].pack_forget()
                    setattr(self, index_attr_name, i)
                    pages[i].pack(fill=tk.BOTH, expand=True)
                    page_label.config(text=f"Page {i + 1} / {len(pages)}")
                    update_nav_buttons()

            prev_btn.config(command=lambda: show_page(getattr(self, index_attr_name) - 1))
            next_btn.config(command=lambda: show_page(getattr(self, index_attr_name) + 1))

            update_nav_buttons()
                


        # =========================
        # TUTORIAL TAB
        # =========================

        bg = "white"

        title = tk.Label(
            tutorial_tab,
            text="Tutorial",
            font=("Arial", 12, "bold"),
            pady=10,
            bg=bg
        )
        title.pack()

        container = tk.Frame(tutorial_tab, bg=bg)
        container.pack(fill=tk.BOTH, expand=True, padx=12, pady=5)

        tutorial_items = [
            ("Step 1 – Load Data",
             "• Enter the file name (SEG-Y or image profile) if it exists in the current directory, or select it via \"Browse\".\n"
             "• Adjust initial parameters according to your file.\n"
             "• Press \"Load File\" to choose an inline or crossline profile from the 3D data.\n"
             "• Select the desired study area and press \"Proceed to Lateral Interpretation\"."
            ),

            ("Step 2 – Manual Interpretation (optional)",
             "• If enabled, you can guide the interpretation manually.\n"
             "• You may revisit manual interpretation of the cropped section via the \"More\" option.\n"
             "• See the Manual Interpretation tab for details."
            ),

            ("Step 3 – Contour Detection",
             "• If Manual Interpretation is enabled, manually interpreted layers or masks will be used.\n"
             "• Otherwise, contours are detected automatically."
            ),

            ("Step 4 – TSP Onlap/Downlap Connection & Vertical Smoothing",
             "• Using default parameters is highly recommended."
            ),

            ("Step 5 – Feature Extraction (optional)",
             "• You may extract features by selecting top and bottom layers."
            ),
            
            ("Notes", 
             "• This method was originally developed for seismic (SEG-Y) data, but it may also be used with seismic or non-seismic image files. Users should apply discretion when interpreting the results."        
            ),            
            
            ("", 
             "• When an image (not SEG-Y) is loaded, the software automatically emulates seismic polarity. In RGB mode, red-dominant pixels become positive amplitudes and blue-dominant pixels become negative, weighted by their intensity. In Grayscale mode, pixels are split into positive or negative values based on a mid-point brightness threshold."
            ),  

            ("", 
             "• This program is designed to detect laterally continuous layers (seismic-like horizons) and does not reliably recognize closed curves."
            ),  

            ("", 
             "• In Manual Interpretation → Drawing Layers, the \"Merge\" button does not outline shapes with similar y-values, which is typical of closed geometries."             
            ),  

            ("", 
             "• Developer Options allow loading intermediate stages (.npz files)."
            ),  
            
            ("Notes (Continued)", 
             "• cropped.npz can be read in Python as:\n"
             "  cropped = (np.load(\"cropped.npz\")[\"data\"]).T"
            ),  

            ("", 
             "• Other .npz files contain interpreted layers and can be read in Python as:\n"
             "  layers = [list(map(tuple, file_npz[key])) for key in file_npz if key != \"flag\"]"
            ),  

            ("", 
             "• Inline/Crossline numbers are stored in cropped.npz and can be read in Python as:\n"
             "  inline = int(cropped_npz[\"inline\"][0]) if \"inline\" in cropped_npz else -999"
            ),  

            ("", 
             "• You can mark individual layers via \"More\" → \"Layer Selection\".\n"
             " (In \"Layer Selection\" window, press on \"Select Layer\" to mark a layer and then press \"Save Selection\" to add that line to the list of saved lines. Finally export saved line(s) by pressing \"Export layer.npz\".)"             
            ),              
            

            ("",
             "• You can generate mask of your profile using minimum and maximum image data thresholds via \"More\" → \"Polarity\".\n"
             "  Polarity mask can be used as Manual Interpretation Type → Mask in the Contour Detection stage."             
            )

        ]


        build_paged_section(
            parent_tab=tutorial_tab,
            container=container,
            items=tutorial_items,
            items_per_page=5,
            index_attr_name="page0_index",
            bg="white"
        )

                
                
        # =========================
        # PARAMETERS TAB
        # =========================

        bg = "white"

        title = tk.Label(
            params_tab,
            text="Parameters",
            font=("Arial", 12, "bold"),
            pady=10,
            bg=bg
        )
        title.pack()

        container = tk.Frame(params_tab, bg=bg)
        container.pack(fill=tk.BOTH, expand=True, padx=12, pady=5)

        param_items = [

            ("Time Sample Interval",
             "Vertical sampling interval of the seismic data (e.g., two-way travel time per sample). "
             "If you use an image profile, the vertical spacing between pixels is used as the reference; therefore, setting this value to 1 is recommended."
             "If you use a non-seismic image profile, make sure to set this parameter to 1 (consecutive vertical spacing)"),

            ("x_start",
             "Starting x-index of the seismic profile. This value defines the horizontal reference origin. "
             "If you use an image profile, please notice the x_start is the zero pixel coordinate of the bottom-left-corner of your image."),

            ("y_start",
             "Starting y-index of the seismic profile. This parameter represents the vertical reference origin (e.g., two-way travel time). "
             "If you use an image profile, please notice the y_start is the zero pixel coordinate of the bottom-left-corner of your image."),

            ("Inline vs. Crossline",
             "Specifies whether the selected profile is an inline or a crossline section. "
             "This setting determines how the 3D seismic volume is sliced into a 2D profile."
             "If you use a non-seismic image, ignore this setting."),

            ("Initial Inline / Crossline",
             "Reference inline or crossline number used as the logical origin of the displayed profile. "
             "This value is preserved in saved results and allows correct spatial referencing in later processing stages."
             "If you use a non-seismic image, ignore this setting."),

            ("Manual Interpretation",
             "Enables or disables Manual Interpretation. If set to Yes, manually interpreted layers or masks "
             "will be incorporated into the Lateral Interpretation process. If set to No, the workflow proceeds fully automatically."),

            ("Relief",
             "Parameter relief manipulates (seismic) image constrast by percentile: low values strongly clip the data, while 100 preserves the full range."),   

            ("Line Length Threshold (l_threshold)",
             "Layers with lengths shorter than l_threshold points are filtered. Very low l_threshold values "
             "produce more detailed interpretations, but increase processing time and the likelihood of errors after the TSP stage."),
             
            ("flood-fill option",
             "Enables flood-fill contour detection. Refer to the author's publication for the underlying connectivity logic."),             

            ("Tiling Window Size",
             "Window size, in points, for each tile of the TSP canvas."),

            ("Tiling Overlap",
             "(important: try to keep this value at 0) Overlap, in points, between two consecutive tiles in the TSP canvas."),

            ("TSP Angle Threshold",
             "Maximum absolute angle, in degrees, between two merging layers."),

            ("TSP Distance Threshold",
             "Maximum distance, in points, between two merging layers."),

            ("N",
             "Number of consecutive points from the start or end of a layer used to calculate the onlap/downlap connection slope."),

            ("r",
             "(important: try to keep this value at 1) Square-offset radius, in points, that increases the chance of merging in onlap/downlap connections."),

            ("x_tol Threshold",
             "(important: try to keep this value at 0) Tolerance along the x-axis for vertical smoothing."),

            ("Sigma Threshold",
             "Gaussian vertical smoothing sigma value."),

            ("Similarity Threshold",
             "Similarity threshold, in percentage, between any two layers. "
             "Layers that are too similar will be merged into an average layer.")
        ]


        build_paged_section(
            parent_tab=params_tab,
            container=container,
            items=param_items,
            items_per_page=5,
            index_attr_name="page_index",
            bg="white"
        )



        # =========================
        # MANUAL TAB
        # =========================
        step_font = tkfont.Font(family="Arial", size=12, weight="bold")
        body_font = tkfont.Font(family="Arial", size=10)
        
        bg = manual_tab.cget("background")

        text = tk.Text(
            manual_tab,
            wrap="word",
            borderwidth=0,
            highlightthickness=0,
            background=bg,
            relief="flat"
        )

        
        
        text.pack(fill=tk.BOTH, expand=True, padx=15, pady=15)
        
        text.tag_configure("step", font=step_font)
        text.tag_configure("body", font=body_font)
        

        
        text.insert("end", "\n", "step")
        
        text.insert("end", "\n", "step")

        text.insert("end", "Location\n", "step")
        text.insert("end",
                    "• You can manually interpret the cropped profile via \"More\" → \"Manual Interpretation\".\n"
                    "• If Manual Interpretation is enabled at startup, you will be automatically directed to the Manual Interpretation window after the cropping stage.\n\n",
                    "body")

        text.insert("end", "Drawing Mode → Layer Drawing\n", "step")
        text.insert("end",
                    "• Click on profile to outline a layer, then press \"Merge\".\n"
                    "• Continue interpretation by clicking on the canvas to add more layers.\n"
                    "• Press \"Finish\" to complete interpretation.\n\n",
                    "body")

        text.insert("end", "Drawing Mode → Mask Drawing\n", "step")
        text.insert("end",
                    "• Click any number of points to define a mask using the minimum and maximum amplitudes of the selected points.\n"
                    "• If interpreted points are spatially separated, later interpretation and plotting may not display them correctly.\n"
                    "• Press \"Finish\" to complete interpretation.\n\n",
                    "body")

        text.insert("end", "Notes:\n", "step")
        text.insert("end",
                    "• Layer Drawing and Mask Drawing results cannot be combined.\n"
                    "• Mask Drawing is mainly intended for development purposes; the recommended manual interpretation method is Layer Drawing.\n"
                    "• You must proceed to the \"Contours\" stage to export or apply Manual Interpretation results.\n"
                    "• To export Manual Interpretation results in the \"Contours\" stage, set Interpretation Setting to \"None\".\n"
                    "• After exporting manual.npz, you can plot manual results via \"More\" → \"Plotting\".\n\n",
                    "body")

        
        
        text.config(state="disabled")



        # =========================
        # PLOTTING TAB
        # =========================
        step_font = tkfont.Font(family="Arial", size=12, weight="bold")
        body_font = tkfont.Font(family="Arial", size=10)
        
        bg = plotting_tab.cget("background")

        text = tk.Text(
            plotting_tab,
            wrap="word",
            borderwidth=0,
            highlightthickness=0,
            background=bg,
            relief="flat"
        )

        
        
        text.pack(fill=tk.BOTH, expand=True, padx=15, pady=15)
        
        text.tag_configure("step", font=step_font)
        text.tag_configure("body", font=body_font)
        

        
        text.insert("end", "\n", "step")
        
        text.insert("end", "\n", "step")        
        
        text.insert("end", "Location\n", "step")
        text.insert("end",
                    "• You can plot all generated profiles via \"More\" → \"Plotting\".\n"
                    "• By following the interpretation workflow, you can export .png plots from smooth.npz or feature.npz; however, these plots do not include axis information.\n\n",
                    "body")

        text.insert("end", "Drawing Plots\n", "step")
        text.insert("end",
                    "• Selecting cropped.npz plots the cropped section.\n"
                    "• Other plots require cropped.npz, so ensure interpretation is completed for the corresponding profile.\n"
                    "• Overlay options allow original profile and interpreted layers to be plotted together.\n\n",
                    "body")

        text.insert("end", "Drawing Interpreted Plots\n", "step")
        text.insert("end",
                    "• All .npz files except cropped.npz and masks.npz store interpreted layers in a consistent format: a list of lists of (x, y) coordinates.\n"
                    "• Exporting plots as .svg allows each layer to be edited as a vector.\n\n",
                    "body")

        text.insert("end", "Overlay Option\n", "step")
        text.insert("end",
                    "• None: No overlays.\n"
                    "• Original: The selected interpreted plot is overlaid on the original profile (cropped.npz).\n"
                    "• Marked Layers: Results from \"Layer Selection\" (layer.npz) are overlaid on the selected (original or interpreted) profile.\n\n",
                    "body")

        text.insert("end", "Notes:\n", "step")
        text.insert("end",
                    "• You must proceed to the \"Contours\" stage to export or apply Manual Interpretation results.\n"
                    "• To use the Marked Layers overlay option, first export marked layers via \"More\" → \"Layer Selection\".\n\n",
                    "body")

        
        
        text.config(state="disabled")


        # =========================
        # CREDIT TAB
        # =========================
        
        bg = "white"

        title = tk.Label(
            credit_tab,
            text="Credits",
            font=("Arial", 12, "bold"),
            pady=10,
            bg=bg
        )
        title.pack()

        container = tk.Frame(credit_tab, bg=bg)
        container.pack(fill=tk.BOTH, expand=True, padx=12, pady=5)

        
        credit_items = [

            (
                "Author and Developer",
                "Hamed Sabagh\n"
                "Department of Geology and Mineralogy\n"
                "Kyoto University, Japan\n"
                "Email: sabagh.hamed.57f@st.kyoto-u.ac.jp\n"
                "GitHub: https://github.com/SABAGHhamed"
            ),

            (
                "Academic Context",
                "This graphical user interface (GUI) was developed during the author’s PhD research "
                "at Kyoto University, titled:\n\n"
                "“Reconstruction of Sea-Level Change from 2D Stratigraphic Sections: "
                "An Inverse Model Using Convolutional Neural Network”.\n\n"
                "The research is conducted under the supervision of Professor Hajime Naruse "
                "(naruse@kueps.kyoto-u.ac.jp), Department of Geology and Mineralogy, "
                "Kyoto University, Kitashirakawaoiwake-cho, Sakyo-ku, Kyoto, "
                "606-8502, Japan."
            ),

            (
                "Related Manuscript",
                "The methodology implemented in this software is described in the manuscript:\n\n"
                "“Lateral Interpretation: A New Method for High-Resolution Seismic Stratigraphy”.\n\n"
                "This software was developed to enable high-resolution stratigraphic interpretation, "
                "which is required for practical application of sea-level reconstruction "
                "on real seismic data."
            ),

            (
                "Development Status and Scope",
                "This project is under active development. The current implementation is primarily "
                "designed for lateral interpretation of seismic-style stratigraphic layers.\n\n"
                "Support for more complex image-based and non-seismic interpretations, "
                "including closed or irregular geometries, may be expanded in future versions, "
                "subject to available research time."
            ),

            (
                "Licensing and Availability",
                "Upon acceptance of the associated manuscript, this software is intended to be "
                "released publicly under the MIT License via the author’s GitHub repository.\n\n"
                "The software may also be provided as supplementary material to the publication."
            ),

            (
                "Third-Party Libraries",
                "This software uses the CIGSEGY library for reading SEG-Y seismic data:\n\n"
                "Li, J. “CIGSEGY: A tool for exchanging data between SEG-Y format and NumPy array "
                "inside Python environment”.\n\n"
                "Project URL: https://github.com/JintaoLee-Roger/cigsegy\n"
                "Documentation: https://cigsegy.readthedocs.io/en/latest/"
            ),

            (
                "Contributions",
                "Contributions, suggestions, and extensions are welcome.\n\n"
                "Users are free to modify or extend the software under the terms of the MIT License, "
                "with or without notifying the author, particularly for applications beyond "
                "seismic interpretation."
            ),
        ]



        build_paged_section(
            parent_tab=credit_tab,
            container=container,
            items=credit_items,
            items_per_page=3,
            index_attr_name="page_index",
            bg="white"
        )








    def browse_sgy_file(self):
        """
        Open native file dialog to select a .sgy/.segy file from any drive.
        """
        file_path = filedialog.askopenfilename(
            title="Select Seismic SEG-Y File or Profile Images",
            filetypes=[
                ("SEG-Y files", "*.sgy *.segy"),
                ("Image files", "*.png *.jpg"),
                ("All files", "*.*")
            ]
        )

        if not file_path:
            return  # user canceled

        # Store full path internally
        self.segy_path = file_path

        # Update entry box
        self.entry.delete(0, tk.END)
        self.entry.insert(0, file_path)


        
        
        
    # === Toggle Developer Options ===
    def toggle_dev_options(self):
        if self.dev_options_visible:
            self.dev_frame.pack_forget()
            self.dev_toggle.config(text="Lateral Interpretation Steps ▼")
        else:
            self.dev_frame.pack(pady=5)
            self.dev_toggle.config(text="Lateral Interpretation Steps ▲")
        self.dev_options_visible = not self.dev_options_visible



        

    # Add this method to your GUI class
    def show_more_menu(self):
        """Display a popup menu with additional options"""
        menu = tk.Menu(self.root, tearoff=0)
        menu.add_command(label="Manual Interpretation", command=self.skip_launch_seismic_mask_gui)        
        menu.add_command(label="Layer Selection", command=self.skip_launch_affine_selector)
        menu.add_command(label="Polarity", command=self.skip_launch_seismic_polarity_gui)
        menu.add_command(label="Plotting", command=self.skip_launch_plotting_gui)

        # Display the menu at the button's location
        x = self.more_button.winfo_rootx()
        y = self.more_button.winfo_rooty() + self.more_button.winfo_height()
        menu.post(x, y)        


    def drop(self, event):
        file_path = event.data.strip("{}")
        self.entry.delete(0, tk.END)
        self.entry.insert(0, file_path)
        self.load_file(file_path)

    def load_from_entry(self):
        filename = self.entry.get().strip()
        if os.path.isfile(filename):
            self.load_file(filename)
        else:
            self.label.config(text="File not found.")

    def load_file(self, file_path):
        if not os.path.isabs(file_path):
            file_path = os.path.join(os.getcwd(), file_path)

        # ---------- Inline input (always) ----------
        temp_root = tk.Toplevel(self.root)
        temp_root.withdraw()

        inline_num = simpledialog.askinteger( 
            "Inline Input",
            "Enter inline / crossline number to view",
            initialvalue=100,
            parent=temp_root
        )

        temp_root.destroy()

        if inline_num is None:   # inline_num is segy friendly not arr friendly inline but self.inline is array friendly
            return

        self.inline = inline_num - int(self.l_threshold_entry.get())

        # ---------- Try seismic first ----------
        try:
            d = cigsegy.fromfile(file_path)
            self.data = d
            self.is_image = False
            self.label.config(text=f"Loaded profile: {os.path.basename(file_path)}")

        except Exception:
            # ---------- Fallback to image ----------
            try:
                img_raw = Image.open(file_path)
                mode = img_raw.mode  # Check if it's RGB, L (grayscale), etc.

                # RGB Logic
                if mode == 'RGB':
                    rgb_arr = np.asarray(img_raw, dtype=np.float32)
                    R = rgb_arr[:, :, 0]
                    B = rgb_arr[:, :, 2]
                    
                    arr = np.zeros(R.shape, dtype=np.float32)
                    not_both_zero = (R != 0) | (B != 0)
                    
                    # This covers both your cases: 
                    # If R > B, result is positive. If B > R, result is negative.
                    arr[not_both_zero] = (R[not_both_zero] - B[not_both_zero]) * (R[not_both_zero] + B[not_both_zero]) / 2


                # Grayscale Logic
                else:
                    img_l = img_raw.convert("L")
                    gray_arr = np.asarray(img_l, dtype=np.float32)
                    
                    # If >= 127.5, keep positive value. If < 127.5, make it negative.
                    arr = np.where(gray_arr >= 127.5, gray_arr, -gray_arr)


                # Assign to your GUI attributes
                self.image_full_relief = arr
                self.seismic_profile_full_relief = arr
                self.seismic_profile = arr.T  # Transposed for display
                self.is_image = True

                self.show_plot()
                self.label.config(text=f"Loaded image: {os.path.basename(file_path)}")
                return

            except Exception as e:
                self.label.config(text=f"Error loading file: {e}")
                return

        # ---------- Seismic slicing ----------
        try:
            crossline_val = self.crossline_var.get()
            seismic_relief = np.clip(self.relief_var.get(), 51, 100)

            if crossline_val == 1:
                slice_2d = d[:, self.inline, :]
            else:
                slice_2d = d[self.inline, :, :]

            self.seismic_profile_full_relief = slice_2d
            vm = np.percentile(slice_2d, seismic_relief)
            self.seismic_profile = np.clip(slice_2d, -vm, vm)

            self.show_plot()

        except Exception as e:
            self.label.config(text=f"Invalid inline number or read error: {e}")

            
    def update_seismic_relief(self): # Thank God, .T is fixed
        if self.seismic_profile_full_relief is None:
            self.label.config(text="Load a file first.")
            return

        seismic_relief = np.clip(self.relief_var.get(), 51, 100)

        try:
            vm = np.percentile(self.seismic_profile_full_relief, seismic_relief)
            self.seismic_profile = np.clip(
                self.seismic_profile_full_relief.T, -vm, vm
            )
            self.show_plot()
        except Exception as e:
            self.label.config(text=f"Error updating relief: {e}")


    def show_plot(self):
        if self.canvas_widget:
            self.canvas_widget.get_tk_widget().destroy()
            self.canvas_widget = None

        if hasattr(self, 'plot_frame') and self.plot_frame:
            self.plot_frame.destroy()

        seismic_T = self.seismic_profile.T

        fig, ax = plot_seismic_GUI(seismic_T, pixels=512)

        self.plot_frame = tk.Frame(self.root)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

        self.canvas_widget = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas_widget.draw()
        self.canvas_widget.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # RectangleSelector setup
        self.selector_coords_y = None
        self.selector_coords_x = None
        self.selector = RectangleSelector(
            ax, self.on_select, useblit=True,
            button=[1], minspanx=5, minspany=5,
            spancoords='pixels', interactive=True
        )

    def on_select(self, eclick, erelease):
        # Note: RectangleSelector returns coordinates in the plotted image coordinate system
        y1 = int(eclick.ydata)
        y2 = int(erelease.ydata)
        x1 = int(eclick.xdata)
        x2 = int(erelease.xdata)
        self.selector_coords_y = (min(y1, y2), max(y1, y2))
        self.selector_coords_x = (min(x1, x2), max(x1, x2))
        self.ok_button.config(state=tk.NORMAL)
        x_start = self.x_start_var.get()
        y_start = self.y_start_var.get()
        
        
        twt_spacing = float(self.rescaling_constant_var.get())
        presented_range_y = (
            y_start + (self.selector_coords_y[0])*twt_spacing,
            y_start + (self.selector_coords_y[1])*twt_spacing
        )
        presented_range_x = (
            self.selector_coords_x[0] + x_start,
            self.selector_coords_x[1] + x_start
        )        

        #messagebox.showinfo("Saved", f"{self.selector_coords_x},{x_start}")
        self.label.config(text=f"Y-range: {presented_range_y}, X-range: {presented_range_x}")

    def finalize_selection(self):
        if not self.selector_coords_y or not self.selector_coords_x:
            self.label.config(text="Please select a region first.")
            return

        plot_y_min, plot_y_max = self.selector_coords_y
        plot_x_min, plot_x_max = self.selector_coords_x

        # Crop on the already-clipped seismic_profile (this mirrors the pipeline)
        cropped = self.seismic_profile[plot_x_min:plot_x_max, plot_y_min:plot_y_max]
        
        
        # writing variables for plotting to self
        self.twt_spacing = float(self.rescaling_constant_var.get())
        self.x_start = self.x_start_var.get() + plot_x_min # this is very important, because not only data might be sliced by general x_start, we sliced a second time by plot_min
        self.y_start = self.y_start_var.get() + (plot_y_min*self.twt_spacing)
        self.crossline_val = self.crossline_var.get()
        self.initial_inline = int(self.l_threshold_entry.get())

        cropped = cropped.T

        # Save cropped data and PNG
        np.savez(f"cropped.npz", inline=np.array([self.inline]), data=cropped)
#         np.savez(f"cropped_{self.inline}.npz", inline=np.array([self.inline]), data=cropped) # for record   # makes overwriting mistakes for crossline
        try:
            plot_seismic_CROPPED(
                                cropped.T,
                                inline=self.inline, # raw inline
                                twt_spacing=self.twt_spacing,
                                x_start=self.x_start,
                                y_start=self.y_start,
                                initial_inline=self.initial_inline,
                                crossline=self.crossline_val,
                                save_option=True,
                                close_option=True,
                                return_fig=False
                                )
        except Exception as e:
            # Non-fatal if plotting helper isn't available
            pass

        # Close this GUI and proceed to the next
        radio_val = self.manual_interp_var.get()
        
        if radio_val == 0:
            try:
                self.root.withdraw()
                launch_seismic_mask_gui(self)
            except Exception as e:
                # If launcher not available, just close
                try:
                    self.root.withdraw()
                except Exception as e:
                    self.root.deiconify()
        else:
            try:
                self.root.withdraw()
                launch_second_gui(self)
            except Exception as e:
                # If launcher not available, just close
                try:
                    self.root.withdraw()
                except Exception as e:
                    self.root.deiconify()


                
    # The skip functions simply close this GUI and launch the downstream GUIs
    def skip_to_contours_processing(self):
        try:
            self.root.withdraw()
            launch_second_gui(self)
        except Exception as e:
            print("Error:", e)
            self.root.deiconify()

    def skip_to_tsp_processing(self):
        try:
            self.root.withdraw()
            launch_second_and_half_gui(self)
        except Exception as e:
            print("Error:", e)
            self.root.deiconify()

    def skip_to_onlap_processing(self):
        try:
            self.root.withdraw()
            launch_third_gui(self)
        except Exception as e:
            print("Error:", e)
            self.root.deiconify()

    def skip_to_smooth_processing(self):
        try:
            self.root.withdraw()
            launch_fourth_gui(self)
        except Exception as e:
            print("Error:", e)
            self.root.deiconify()

    def skip_to_feature_processing(self):
        try:
            self.root.withdraw()
            launch_fifth_gui(self)
        except Exception as e:
            print("Error:", e)
            self.root.deiconify()
         
    def skip_launch_affine_selector(self):
        try:
            # Load and map data            
            self.root.withdraw()
            launch_Layer_selector(self) 
        except Exception as e:
            print("Error:", e)
            self.root.deiconify()

    def skip_launch_seismic_mask_gui(self):
        try:
            # Load and map data            
            self.root.withdraw()
            launch_seismic_mask_gui(self) 
        except Exception as e:
            print("Error:", e)
            self.root.deiconify()
            
    def skip_launch_seismic_polarity_gui(self):
        try:
            # Load and map data            
            self.root.withdraw()
            launch_seismic_polarity_gui(self) 
        except Exception as e:
            print("Error:", e)
            self.root.deiconify()
            
    def skip_launch_plotting_gui(self):
        try:
            # below will run if click and drag event is not performed
            if self.x_start is None:
                self.x_start = self.x_start_var.get()
            if self.y_start is None:    
                self.y_start = self.y_start_var.get()
            if self.crossline_val is None:
                self.crossline_val = self.crossline_var.get()
            if self.twt_spacing is None:
                self.twt_spacing = self.rescaling_constant_var.get()
            if self.initial_inline is None:
                self.initial_inline = int(self.l_threshold_entry.get())

            self.root.withdraw()
            launch_plotting_gui(self)

        except Exception as e:
            print("Error:", e)
            self.root.deiconify()



        
#contours        
def launch_second_gui(app):
    """Second GUI: full contour extraction pipeline (updated).

    Loads cropped.npz (expects the array saved as `data` — this should be the transposed cropped block
    consistent with the rest of the pipeline). The GUI accepts:
      - Relief (used for p12 percentile: relief - 10 and relief)
      - Line length threshold (l_threshold)
      - Rescaling constant (rescaling_constant)

    The processing follows the pipeline provided: extract p/n masks, extract p12 interval masks,
    filter p12 by length, combine coordinates as p_p12 + p + n_p12 + n, rescale by `dim * rescaling_constant`,
    integer-fit, filter by length (`l_threshold * rescaling_constant`) and save contours.npz.
    """

    app.contour_root = tk.Toplevel(app.root)
    app.contour_root.title("Contour Detection")    
    center_window(app.contour_root, 300, 400)                
    frm = tk.Frame(app.contour_root)
    frm.pack(padx=10, pady=8)

    tk.Label(frm, text="").grid(row=0, column=0, sticky=tk.W)   # Relief (percentile, e.g. 99): <= hidden
    relief_entry = tk.Entry(frm, width=8)
    relief_entry.insert(0, "99")
    relief_entry.grid(row=0, column=1, padx=6, pady=3)
    # hide it
    relief_entry.grid_remove()    

    tk.Label(frm, text="Line Length Threshold (l_threshold):").grid(row=1, column=0, sticky=tk.W)
    l_threshold_entry = tk.Entry(frm, width=8)
    l_threshold_entry.insert(0, "50")
    l_threshold_entry.grid(row=1, column=1, padx=6, pady=3)

    tk.Label(frm, text="").grid(row=2, column=0, sticky=tk.W)  # Rescaling constant: was removed from text # placeholder
    rescaling_entry = tk.Entry(frm, width=8)
    rescaling_entry.insert(0, "1.0")
    rescaling_entry.grid(row=2, column=1, padx=6, pady=3)
    # hide it
    rescaling_entry.grid_remove()
    
    # First set of radio buttons
    tk.Label(app.contour_root, text="Use Manual Interpretation?").pack(pady=(10, 2))
    use_mask_var = tk.IntVar(value=0)  # default No
    frame_mask = tk.Frame(app.contour_root)
    frame_mask.pack(pady=2)
    tk.Radiobutton(frame_mask, text="Yes", variable=use_mask_var, value=1).pack(side="left", padx=5)
    tk.Radiobutton(frame_mask, text="No", variable=use_mask_var, value=0).pack(side="left", padx=5)

    # Second set - Manual Interpretation Type
    tk.Label(app.contour_root, text="Manual Interpretation Type").pack(pady=(10, 2))
    use_mask_var2 = tk.IntVar(value=0)  # default Layer
    frame_mask2 = tk.Frame(app.contour_root)
    frame_mask2.pack(pady=2)
    rb2_mask = tk.Radiobutton(frame_mask2, text="Mask", variable=use_mask_var2, value=1, state="disabled")
    rb2_mask.pack(side="left", padx=5)
    rb2_layer = tk.Radiobutton(frame_mask2, text="Layer", variable=use_mask_var2, value=0, state="disabled")
    rb2_layer.pack(side="left", padx=5)

    # Third set - Interpretation Setting
    tk.Label(app.contour_root, text="Lateral Interpretation Setting").pack(pady=(10, 2))
    use_mask_var3 = tk.IntVar(value=0)  # default Full
    frame_mask3 = tk.Frame(app.contour_root)
    frame_mask3.pack(pady=2)
    rb3_none = tk.Radiobutton(frame_mask3, text="None", variable=use_mask_var3, value=2, state="disabled")
    rb3_none.pack(side="left", padx=5)
    rb3_simple = tk.Radiobutton(frame_mask3, text="Simple", variable=use_mask_var3, value=1, state="disabled")
    rb3_simple.pack(side="left", padx=5)
    rb3_full = tk.Radiobutton(frame_mask3, text="Full", variable=use_mask_var3, value=0, state="disabled")
    rb3_full.pack(side="left", padx=5)
    
    
    # --- New Checkbox: Flood-fill Option ---
    app.floodfill_option = tk.BooleanVar(value=True)

    # Create the checkbutton
    chk_flood = tk.Checkbutton(
        app.contour_root, 
        text="flood-fill option", 
        variable=app.floodfill_option
    )
    chk_flood.pack(pady=(10, 2))    
    
    
    def run_contour_pipeline():
        # Parse inputs with validation
        try:
            relief = int(relief_entry.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Relief must be an integer (e.g. 99).")
            return

        try:
            l_threshold = float(l_threshold_entry.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Threshold must be a number.")
            return

        try:
            rescaling_constant = float(rescaling_entry.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Rescaling constant must be a number (e.g. 1.0).")
            return

        # Load cropped data
        try:
            data_npz = np.load("cropped.npz")
            if 'data' not in data_npz:
                messagebox.showerror("File error", "cropped.npz does not contain key 'data'.")
                return
            arr = data_npz['data']
            if arr.ndim != 2:
                messagebox.showerror("Data error", "cropped.npz['data'] must be a 2D array.")
                return
        except FileNotFoundError:
            messagebox.showerror("File missing", "cropped.npz not found. Run the first GUI crop step first.")
            return
        except Exception as e:
            messagebox.showerror("Load error", str(e))
            return

        try:
            # --- Contour extraction for main masks ---

            
            length_limit = l_threshold * rescaling_constant                     
            mask_p, mask_n = get_positive_negative_masks(arr)
            floodfill_option = app.floodfill_option.get()


            def floodfill_verticalpolarity(arr, length_limit, floodfill_option=True, option='arr', layer_arr=None):
                if option=='arr':
                    ### vertical polarity method
                    compacted_image   = arrRB(arr)
                    digofile_positive = digofile(compacted_image, negative_sign = False)
                    digofile_negative = digofile(compacted_image, negative_sign = True)                
                    stripe_coordinates_p = cv2_layers(digofile_positive)
                    stripe_coordinates_n = cv2_layers(digofile_negative)
                    stripe_coordinates = stripe_coordinates_p + stripe_coordinates_n
                    # mask_n
                    mask_p, mask_n = get_positive_negative_masks(arr)
                    stripe_coordinates = filter_lines_by_length(stripe_coordinates, length_limit) 
                    stripe_coordinates = split_lines_by_y(stripe_coordinates)                    

                    if floodfill_option:
                        ### flood-fill method
                        outlines_and_blobs = find_puddles_recursive(mask_n, connectivity=4)  
                        stripe_coordinates_n = outlines_and_blobs

                        stripe_coordinates_n = split_lines_by_y(stripe_coordinates_n)
                        stripe_coordinates_n = integer_fitting(stripe_coordinates_n)
                        stripe_coordinates_n = flip_layers(stripe_coordinates_n)
                        min_value = min(arr.shape[0], arr.shape[1])
                        stripe_coordinates_n = filter_lines_by_length(
                            stripe_coordinates_n, min(min_value, length_limit)
                        )
                        stripe_coordinates_n = remove_subset_lines(stripe_coordinates_n)
                        

                        ## for paper
                        # extracting pure outlines and blobs
            #             outlines = find_puddles_cv2(mask_n, connectivity=4)
            #             outlines = custome_floodfill_process(outlines, arr, length_limit)
            #             outlines_and_blobs = custome_floodfill_process(outlines_and_blobs, arr, length_limit=2)  # length_limit=2 to keep blobs
            #             blobs = exclude_sub_lines(outlines_and_blobs, outlines) 

            #           ##
            #             np_data = [np.array(line) for line in stripe_coordinates]
            #             np.savez("classic.npz", *np_data)   
            #             np_data = [np.array(line) for line in stripe_coordinates_n]
            #             np.savez("floodfill.npz", *np_data)   
                        ##                     
                        stripe_coordinates_n += stripe_coordinates           
                    
                elif option=='mask':  # manual mask
                    ### vertical polarity method                
                    compacted_image   = arrRB(arr)
                    digofile_negative = digofile(compacted_image, negative_sign = True)                
                    stripe_coordinates = cv2_layers(digofile_negative)
                    # mask_n
                    mask_n = arr
                    stripe_coordinates = split_lines_by_y(stripe_coordinates)                    

                    if floodfill_option:
                        ### flood-fill method
                        outlines_and_blobs = find_puddles_recursive(mask_n, connectivity=4)  
                        stripe_coordinates_n = outlines_and_blobs

                        stripe_coordinates_n = split_lines_by_y(stripe_coordinates_n)
                        stripe_coordinates_n = integer_fitting(stripe_coordinates_n)
                        stripe_coordinates_n = flip_layers(stripe_coordinates_n)
                        stripe_coordinates_n = remove_subset_lines(stripe_coordinates_n)
                        
                        stripe_coordinates_n += stripe_coordinates           

                elif option=='layer':  # manual layer
                    # No contour detection, only polishing things up
                    stripe_coordinates = arr # arr is actually stripe_coordinates type not an array                   
                    stripe_coordinates = flip_layers(stripe_coordinates)  # flipping vertically for better plotting
                    stripe_coordinates = sort_lines_by_x(stripe_coordinates)  # new; tricky little think I lost!
                    stripe_coordinates = split_lines_by_x_trend(stripe_coordinates)
                    stripe_coordinates = correct_lines_trend(stripe_coordinates)  # reverse any reversely incrementing layer
                    stripe_coordinates = remove_duplicate_x_points(stripe_coordinates)
                    stripe_coordinates = bresenham_line(stripe_coordinates)
                    stripe_coordinates = filter_lines_by_length(stripe_coordinates, length_limit) 
                    stripe_coordinates_n = split_lines_by_y(stripe_coordinates)                    

                return stripe_coordinates_n


            ## debug~
#             root = tk.Tk()
#             root.title("Plot Display")
#             # Create your figure
#             fig = plt.Figure(figsize=(6, 6))
#             ax = fig.add_subplot(111)
#             ax.clear()
#             ax.imshow(mask_n, cmap='RdBu', origin='upper', aspect='auto')

#             # Embed the figure in tkinter
#             canvas = FigureCanvasTkAgg(fig, master=root)
#             canvas.draw()
#             canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)        
            # ~

            ## old cv2 contour method            
#             alternation_mask = get_polarity_alternation_mask(arr) # new            
#             stripe_coordinates = cv2_layers(alternation_mask)

            use_mask = use_mask_var.get()   # 1 = Yes, 0 = No
            use_mask2 = use_mask_var2.get()  # 1 = Mask, 0 = Layer
            use_mask3 = use_mask_var3.get()  # 2 = None,1 = Simple, 0 = Full
            
            if  use_mask3 == 0:
                
                stripe_coordinates_n = floodfill_verticalpolarity(arr, length_limit, floodfill_option=floodfill_option, option='arr')
            
            try:
                if use_mask == 0:
                    np_data = [np.array(line) for line in stripe_coordinates_n]
                    np.savez("contours.npz", flag=np.array([0]), *np_data) 
                    messagebox.showinfo("Saved", f"Saved contours.npz with {len(stripe_coordinates_n)} lines (no mask used).")
#                     ## for paper
#                     np_data = [np.array(line) for line in outlines]
#                     np.savez("outlines.npz", *np_data) 
#                     np_data = [np.array(line) for line in blobs]
#                     np.savez("blobs.npz", *np_data)                     
#                    ##

                elif use_mask == 1 and use_mask3 == 0:  # Full method
                    np_data = [np.array(line) for line in stripe_coordinates_n]
                    np.savez("contours.npz", flag=np.array([1]), *np_data)   # notice 1 flag of contours_n
                
                    if use_mask2 == 1:  # Masks
                        data_npz = np.load("masks.npz")
                        mask = data_npz["masks"] 
                        stripe_coordinates = floodfill_verticalpolarity(mask, length_limit, floodfill_option=floodfill_option, option='mask')
                        
                    elif use_mask2 == 0:  # Layers
                        np_data = np.load("layers.npz") 
                        layers = [list(map(tuple, np_data[key])) for key in np_data]
                        stripe_coordinates = floodfill_verticalpolarity(
                                                                        layers,
                                                                        length_limit,
                                                                        option='layer',
                                                                        layer_arr=arr
                                                                        )
                        
                    np_data = [np.array(line) for line in stripe_coordinates]
                    np.savez("manual.npz", *np_data)
                    messagebox.showinfo("Saved", "Saved manual.npz with manual data for full interpretaion.")


                elif use_mask == 1 and use_mask3 == 1:  # Simple method

                    if use_mask2 == 1:  # Masks
                        data_npz = np.load("masks.npz")
                        mask = data_npz["masks"] 
                        stripe_coordinates = floodfill_verticalpolarity(mask, length_limit, floodfill_option=floodfill_option, option='mask')
                        
                    elif use_mask2 == 0:  # Layers
                        np_data = np.load("layers.npz") 
                        layers = [list(map(tuple, np_data[key])) for key in np_data]
                        stripe_coordinates = floodfill_verticalpolarity(
                                                                        layers,
                                                                        length_limit,
                                                                        option='layer',
                                                                        layer_arr=arr
                                                                        )


                    np_data = [np.array(line) for line in stripe_coordinates]  # <= notice n_mask is not used, contours.npz is just high amplitudes
                    np.savez("contours.npz", flag=np.array([0]), *np_data)   # notice 0 flag of contours (only highs will be processed)

                    np_data = [np.array(line) for line in stripe_coordinates] # duplicate of above for safe keepings
                    np.savez("manual.npz", *np_data)
                    messagebox.showinfo("Saved", "Saved only manual data for simple interpretaion.")


                elif use_mask == 1 and use_mask3 == 2:  # None; fast export of manual results

                    if use_mask2 == 1:  # Masks
                        data_npz = np.load("masks.npz")
                        mask = data_npz["masks"] 
                        stripe_coordinates = floodfill_verticalpolarity(mask, length_limit, floodfill_option=floodfill_option, option='mask')
                        
                    elif use_mask2 == 0:  # Layers
                        np_data = np.load("layers.npz") 
                        layers = [list(map(tuple, np_data[key])) for key in np_data]
                        stripe_coordinates = floodfill_verticalpolarity(
                                                                        layers,
                                                                        length_limit,
                                                                        option='layer',
                                                                        layer_arr=arr
                                                                        )

                    flbm = stripe_coordinates
                    np_data = [np.array(line) for line in flbm]
                    np.savez(f"manual.npz", *np_data)

                    print(f"Contour detection saved {len(flbm)} lines")

                    dpi = 100
                    width = 512
                    height = 512
                    fig_size = (width / dpi, height / dpi)
                    fig = Figure(figsize=fig_size, dpi=dpi)
                    ax = fig.add_subplot(111)

                    all_x, all_y = [], []

                    def process(layer):
                        unique = list(set(layer))
                        return sorted(unique, key=lambda pt: pt[0])

                    for line in flbm:
                        cleaned = process(line)
                        if cleaned:
                            x = [pt[0] for pt in cleaned]
                            y = [pt[1] for pt in cleaned]
                            all_x.extend(x)
                            all_y.extend(y)
                            ax.plot(x, y, c="black", linewidth=1)

                    if all_x and all_y:
                        ax.set_xlim(min(all_x), max(all_x))
                        ax.set_ylim(min(all_y), max(all_y))

                    ax.axis("off")
                    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
                    fig.savefig(f"manual.png", bbox_inches="tight", pad_inches=0, dpi=dpi)
                    print("Saved manual.png")

                    messagebox.showinfo("Saved", "manual.npz and image of profile has been exported.")



            except Exception as e:
                messagebox.showerror("Error", f"Failed to save contours.npz\n\n{e}")
                                


            app.contour_root.destroy()
            app.root.deiconify()

            # Try to launch next GUI (try several known names in order) (if not None interpretation)
            if use_mask3 != 2:
                try:
                    app.root.after(0, lambda: launch_second_and_quarter_gui(app))
                except Exception as e:
                    pass


        # except Exception as e:
            # messagebox.showerror("Processing Error", str(e))
            
        except Exception as e:
            # Check if the error is specifically a list/array index error
            if isinstance(e, IndexError):
                custom_error = (
                    "The Length Threshold is probably too large to detect any contour. "
                    "Please reduce Length Threshold. If even l_threshold = 2 did not work, "
                    "the profile cannot be interpreted."
                )
                messagebox.showerror("Threshold Error", f"Failed to save contours.npz\n\n{custom_error}")
            else:
                # Fallback for all other types of errors
                messagebox.showerror("Error", f"Failed to save contours.npz\n\n{e}")            

    

    # Function to enable/disable dependent radio buttons
    def toggle_dependent_buttons(*args):
        if use_mask_var.get() == 1:
            # Enable buttons when "Yes" is selected
            state = "normal"
        else:
            # Disable buttons and reset to defaults when "No" is selected
            state = "disabled"
            use_mask_var2.set(0)  # Reset to Layer
            use_mask_var3.set(0)  # Reset to Full

        rb2_mask.config(state=state)
        rb2_layer.config(state=state)
        rb3_none.config(state=state)
        rb3_simple.config(state=state)
        rb3_full.config(state=state)
    
    # Attach the callback to the first variable
    use_mask_var.trace_add("write", toggle_dependent_buttons)

    # Initialize the state
    toggle_dependent_buttons()

    run_button = tk.Button(app.contour_root, text="Run contour extraction & save contours.npz", command=run_contour_pipeline)
    run_button.pack(pady=10)
    
    # --- Back to Home handling ---
    def back_to_home():
        app.contour_root.destroy()
        app.root.deiconify()

    back_btn = tk.Button(app.contour_root, text="⬅ Back to Home", command=back_to_home)
    back_btn.pack(pady=8)    




#quarter
def launch_second_and_quarter_gui(app):
    merging_indices = []
    contours_npz = np.load("contours.npz")
    contours_lines = [list(map(tuple, contours_npz[key])) for key in contours_npz if key != "flag"]
    flag = contours_npz["flag"][0] if "flag" in contours_npz else None    
    
    if flag == 1:    
        manual_npz = np.load("manual.npz")
        manual = [list(map(tuple, manual_npz[key])) for key in manual_npz if key != "flag"]            
        if contours_lines != manual:  # this ensures Simple mode is not activated because in Simple contours.npz and manual.npz are the same
            contours_lines+= manual
            
    # making contours_lines bresenham for easier click event
    contours_lines = bresenham_line(contours_lines)



    def update_canvas():
        fig, ax = plot_layers_GUI(contours_lines, marked_layers=[contours_lines[i] for i in merging_indices])
        canvas.figure = fig
        canvas._tkcanvas.destroy()
        new_canvas = FigureCanvasTkAgg(fig, master=app.quarter_root)
        new_canvas.get_tk_widget().grid(row=0, column=1, columnspan=4)
        new_canvas.draw()
        canvas.__dict__.update(new_canvas.__dict__)
        canvas.mpl_connect("button_press_event", on_canvas_click)

    def on_canvas_click(event):
        inv = ax.transData.inverted()
        x, y = inv.transform((event.x, event.y))
        clicked = (round(x), round(y))
        for idx, line in enumerate(contours_lines):
            if any(abs(pt[0] - clicked[0]) <= 1 and abs(pt[1] - clicked[1]) <= 1 for pt in line):
                merging_indices.append(idx)
                merging_indices[:] = list(np.unique(merging_indices))
                update_canvas()
                break

    def merge_layers():
        if not merging_indices:
            messagebox.showinfo("No Selection", "No layers selected for merging.")
            return
        first_idx = merging_indices[0]
        for idx in merging_indices[1:]:
            contours_lines[first_idx].extend(contours_lines[idx])
        contours_lines[first_idx].sort(key=lambda pt: pt[0]) # sorting x order
        contours_lines[first_idx] = remove_duplicate_x_points(contours_lines[first_idx])
        contours_lines[first_idx] = bresenham_line(contours_lines[first_idx])
        for idx in sorted(merging_indices[1:], reverse=True):
            del contours_lines[idx]
        merging_indices.clear()
        update_canvas()

    def undo_last_merge():
        if merging_indices:
            merging_indices.pop()
            update_canvas()

    def finalize_and_proceed():
        np_data = [np.array(line) for line in contours_lines]
        if flag == 1: # it works fine with "Simple" mode in Contour detection of manual interpretation
            np.savez("contours.npz", flag=np.array([1]), *np_data) 
        else:
            np.savez("contours.npz", flag=np.array([0]), *np_data) 
            
        app.quarter_root.destroy()
        launch_second_and_half_gui(app)
        
    def remove_layer():
        if not merging_indices:
            messagebox.showinfo("No Selection", "No layers selected for merging.")
            return
        for idx in sorted(merging_indices, reverse=True):
            del contours_lines[idx]
        merging_indices.clear()
        update_canvas()        
        

    app.quarter_root = tk.Toplevel(app.root)
    app.quarter_root.title("Contour Edit GUI")
    #center_window(app.quarter_root)                    
    
    fig, ax = plot_layers_GUI(contours_lines)
    canvas = FigureCanvasTkAgg(fig, master=app.quarter_root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.grid(row=0, column=1, columnspan=4)
    canvas.mpl_connect("button_press_event", on_canvas_click)

    merge_btn = tk.Button(app.quarter_root, text="Remove Layers", command=remove_layer)  # 
    merge_btn.grid(row=1, column=1)

    undo_btn = tk.Button(app.quarter_root, text="Undo Last Merge", command=undo_last_merge)
    undo_btn.grid(row=1, column=2)

    perform_merge_btn = tk.Button(app.quarter_root, text="Perform Merge", command=merge_layers)
    perform_merge_btn.grid(row=1, column=3)

    finalize_btn = tk.Button(app.quarter_root, text="Finish Merging", command=finalize_and_proceed)
    finalize_btn.grid(row=1, column=4)
    
    # --- Back to Home handling ---
    def back_to_home():
        app.quarter_root.destroy()
        app.root.deiconify()

    back_btn = tk.Button(app.quarter_root, text="⬅ Back to Home", command=back_to_home)
    back_btn.grid(row=2, column=1, columnspan=4, pady=8)

    

        #messagebox.showinfo("Imporant", "This stage does not show manual interpretation results!")



    
    


#tsp    
def launch_second_and_half_gui(app):
    def run_tsp_pipeline():
        try:
            magnet_threshold = int(magnet_entry.get())
            tiling_window_threshold = int(tiling_window_entry.get())
            tiling_overlap_threshold = int(tiling_overlap_entry.get())
            tsp_angle_thresh = int(tsp_angle_thresh_entry.get())
            tsp_threshold_distance = int(tsp_dist_thresh_entry.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Parameters must be integers.")
            return

        try:
            contours = np.load("contours.npz")
            contours_ = [list(map(tuple, contours[k])) for k in contours if k != "flag"]
                                
            flag = contours["flag"][0] if "flag" in contours else None
    
            if flag == 1 or flag == 0:  # or flag == 0 makes it contain all cases; originally it was if flag == 1 but we decided that quarter gui does that
                # manual = np.load("manual.npz")
                # manual_ = [list(map(tuple, manual[key])) for key in manual]

                # M = magnet(manual_, magnet_threshold, distance_threshold = 1)
                # M = remove_duplicate_x_points(M)
                # M_bre = bresenham_line(M)

                # recon_lines = tiled_and_tsp_merge(
                    # M_bre,
                    # window_size=tiling_window_threshold,
                    # overlap=tiling_overlap_threshold,
                    # angle_thresh=tsp_angle_thresh,
                    # threshold_distance=tsp_threshold_distance
                # )

                # recon_lines = remove_duplicate_x_points(recon_lines)
                # recon_lines = [line for line in recon_lines if line]
                # recon_lines = conjoin_lines(recon_lines)
                # recon_lines = remove_subset_lines(recon_lines)
# #                 plot_layers(recon_lines, saving_directory="debugging_gif")     # for debug
                # print(f"recon {len(recon_lines)}")
                
                # contours_ += recon_lines   # (important) adding manual tsp to contours
                
                # M = magnet(contours_, magnet_threshold, distance_threshold = 1)
                # M = remove_duplicate_x_points(M)
                # M_bre = bresenham_line(M)

                # recon_lines = tiled_and_tsp_merge(
                    # M_bre,
                    # window_size=tiling_window_threshold,
                    # overlap=tiling_overlap_threshold,
                    # angle_thresh=tsp_angle_thresh,
                    # threshold_distance=tsp_threshold_distance
                # )

                # recon_lines = remove_duplicate_x_points(recon_lines)
                # recon_lines = [line for line in recon_lines if line]
                # recon_lines = conjoin_lines(recon_lines)
                # recon_lines = remove_subset_lines(recon_lines)
                # print(f"recon {len(recon_lines)}")
                                
            # else:            
                # print("no flag")
                M = magnet(contours_, magnet_threshold, distance_threshold = 1)
                M = remove_duplicate_x_points(M)
                M_bre = bresenham_line(M)

                recon_lines = tiled_and_tsp_merge(
                    M_bre,
                    window_size=tiling_window_threshold,
                    overlap=tiling_overlap_threshold,
                    angle_thresh=tsp_angle_thresh,
                    threshold_distance=tsp_threshold_distance
                )

                recon_lines = remove_duplicate_x_points(recon_lines)
                recon_lines = [line for line in recon_lines if line]
                recon_lines = conjoin_lines(recon_lines)
                recon_lines = remove_subset_lines(recon_lines)
                print(f"recon {len(recon_lines)}")
                 
                    
###                     
            np_data = [np.array(line) for line in recon_lines]  # modified
            np.savez("tsp.npz", *np_data)


            messagebox.showinfo("Done", f"Saved tsp.npz with {len(recon_lines)} lines.")
            app.tsp_root.destroy()
            launch_third_gui(app)

        except Exception as e:
            messagebox.showerror("Processing Error", str(e))
            

    app.tsp_root = tk.Toplevel(app.root)
    app.tsp_root.title("TSP")
    center_window(app.tsp_root, 200, 370)

    tk.Label(app.tsp_root, text="").pack(pady=5)  # Magnet Threshold placeholder
    magnet_entry = tk.Entry(app.tsp_root)
    magnet_entry.insert(0, "0")
    magnet_entry.pack(pady=5)
    # hide from user
    magnet_entry.pack_forget()

    tk.Label(app.tsp_root, text="Tiling Window Size").pack(pady=2)
    tiling_window_entry = tk.Entry(app.tsp_root)
    tiling_window_entry.insert(0, "200")
    tiling_window_entry.pack(pady=2)

    tk.Label(app.tsp_root, text="Tiling Overlap").pack(pady=2)
    tiling_overlap_entry = tk.Entry(app.tsp_root)
    tiling_overlap_entry.insert(0, "0")
    tiling_overlap_entry.pack(pady=2)

    tk.Label(app.tsp_root, text="TSP Angle Threshold").pack(pady=2)
    tsp_angle_thresh_entry = tk.Entry(app.tsp_root)
    tsp_angle_thresh_entry.insert(0, "80")
    tsp_angle_thresh_entry.pack(pady=2)

    tk.Label(app.tsp_root, text="TSP Distance Threshold").pack(pady=2)
    tsp_dist_thresh_entry = tk.Entry(app.tsp_root)
    tsp_dist_thresh_entry.insert(0, "210")
    tsp_dist_thresh_entry.pack(pady=2)

    run_button = tk.Button(app.tsp_root, text="Run TSP and Save tsp.npz", command=run_tsp_pipeline)
    run_button.pack(pady=10)
    
    # --- Back to Home handling ---
    def back_to_home():
        app.tsp_root.destroy()
        app.root.deiconify()

    back_btn = tk.Button(app.tsp_root, text="⬅ Back to Home", command=back_to_home)
    back_btn.pack(pady=8)    






###################################

#onlap
def launch_third_gui(app):  # spstR
    def run_pipeline_2():
        try:
            # Get threshold values from entry fields
            N = int(N_entry.get())
            r = int(r_entry.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Both N and r must be numeric values.")
            return

        try:
            # === tsp detection ===
            np_data = np.load("tsp.npz")
            recon_lines = [list(map(tuple, np_data[key])) for key in np_data]

            # Run spstR with user-defined parameters
            spstR_result = spstR(
                recon_lines,
                N=N,
                r=r,
                debug=False,
                DEBUG_it=50,
                debug_tpath=False
            )

            # Save the result
            np_data = [np.array(line) for line in spstR_result]
            np.savez("onlap.npz", *np_data)

            messagebox.showinfo("All Done", f"Saved onlap.npz with {len(spstR_result)} lines.")
            app.onlap_root.destroy()  # Close GUI

            # Launch next GUI
            launch_fourth_gui(app)

        except Exception as e:
            messagebox.showerror("Processing Error", str(e))

    # === GUI window ===
    app.onlap_root = tk.Toplevel(app.root)
    app.onlap_root.title("Onlap/Downlap Connection")
    center_window(app.onlap_root, 400, 230)    

    # N input
    tk.Label(app.onlap_root, text="N (default: 100)").pack(pady=2)
    N_entry = tk.Entry(app.onlap_root)
    N_entry.insert(0, "100")
    N_entry.pack(pady=2)

    # r input
    tk.Label(app.onlap_root, text="r (default: 1)").pack(pady=2)
    r_entry = tk.Entry(app.onlap_root)
    r_entry.insert(0, "1")
    r_entry.pack(pady=2)

    # Run button
    run_button = tk.Button(app.onlap_root, text="Run onlap/downlap connection and Save", command=run_pipeline_2)
    run_button.pack(pady=10)
    
    # --- Back to Home handling ---
    def back_to_home():
        app.onlap_root.destroy()
        app.root.deiconify()

    back_btn = tk.Button(app.onlap_root, text="⬅ Back to Home", command=back_to_home)
    back_btn.pack(pady=8)    


    
    
########################################

#smooth
def launch_fourth_gui(app):
    def run_pipeline_3():
        try:
            # Get threshold from entry field
            x_tol_threshold = float(x_tol_entry.get())
            sigma_threshold = float(sigma_entry.get())
            similarity_threshold = float(similarity_entry.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Threshold must be a number.")
            return

        # Load and process data
        try:
            # === onlap detection ===
            np_data = np.load("onlap.npz")
            R = [list(map(tuple, np_data[key])) for key in np_data]

            msl = remove_duplicate_x_points(R)
            msl = bresenham_line(msl)
            msl = merge_similar_lines_keep_longest(msl, threshold=similarity_threshold)
            msl = extend_critical_lines(msl)
            msl = remove_duplicate_x_points(msl)
            smsl = smooth_lines(msl, x_tol=x_tol_threshold, sigma=sigma_threshold)

            print(f"smoothed {len(smsl)}")

            # Save smooth result
            np_data = [np.array(line) for line in smsl]
            np.savez("smooth.npz", *np_data)

            messagebox.showinfo("All Done", f"Saved smooth.npz with {len(smsl)} lines.")

            # Check radio selection
            if action_var.get() == 0:
                dpi = 100
                width = 512
                height = 512
                flbm = smsl
                try:
                    data_npz = np.load("cropped.npz")
                    inline = data_npz["inline"][0]
                except:
                    inline = "no_recorded_inline"

                fig_size = (width / dpi, height / dpi)
                fig = Figure(figsize=fig_size, dpi=dpi)
                ax = fig.add_subplot(111)

                all_x, all_y = [], []

                def process(layer):
                    unique = list(set(layer))
                    return sorted(unique, key=lambda pt: pt[0])

                for line in flbm:
                    cleaned = process(line)
                    if cleaned:
                        x = [pt[0] for pt in cleaned]
                        y = [pt[1] for pt in cleaned]
                        all_x.extend(x)
                        all_y.extend(y)
                        ax.plot(x, y, c="black", linewidth=1)

                if all_x and all_y:
                    ax.set_xlim(min(all_x), max(all_x))
                    ax.set_ylim(min(all_y), max(all_y))

                ax.axis("off")
                fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
                fig.savefig(f"feature_{inline}.png", bbox_inches="tight", pad_inches=0, dpi=dpi)
                print("Saved feature.png")                  

                # Save smooth result
                np_data = [np.array(line) for line in flbm]
                np.savez(f"feature.npz", *np_data)                
#                 np.savez(f"feature_{inline}.npz", *np_data)  # makes overwriting mistakes for crossline
                
                app.smooth_root.destroy()
                app.root.deiconify()
                return
            else:
                app.smooth_root.destroy()
                launch_fifth_gui(app)

        except Exception as e:
            messagebox.showerror("Processing Error", str(e))

    # === GUI Setup ===
    app.smooth_root = tk.Toplevel(app.root)
    app.smooth_root.title("Vertical Smoothing")
    center_window(app.smooth_root, 200, 370)                        
    tk.Label(app.smooth_root, text="x_tol Threshold (e.g., 0):").pack(pady=5)
    x_tol_entry = tk.Entry(app.smooth_root)
    x_tol_entry.insert(0, "0")
    x_tol_entry.pack(pady=5)

    tk.Label(app.smooth_root, text="Sigma Threshold (e.g., 2.0):").pack(pady=5)
    sigma_entry = tk.Entry(app.smooth_root)
    sigma_entry.insert(0, "2.0")
    sigma_entry.pack(pady=5)

    tk.Label(app.smooth_root, text="Similarity Threshold (default: 0.75)").pack(pady=2)
    similarity_entry = tk.Entry(app.smooth_root)
    similarity_entry.insert(0, "0.75")
    similarity_entry.pack(pady=2)

    # === Radio Button for Export/Proceed ===
    tk.Label(app.smooth_root, text="Choose Action:").pack(pady=(10, 2))
    action_var = tk.IntVar(value=0)  # Default: export smooth profile
    tk.Radiobutton(app.smooth_root, text="Export smooth profile", variable=action_var, value=0).pack(anchor="w", padx=20)
    tk.Radiobutton(app.smooth_root, text="Proceed to feature selection", variable=action_var, value=1).pack(anchor="w", padx=20)

    # === Run Button ===
    run_button = tk.Button(app.smooth_root, text="Run smoothing and Save", command=run_pipeline_3)
    run_button.pack(pady=10)
    
    # --- Back to Home handling ---
    def back_to_home():
        app.smooth_root.destroy()
        app.root.deiconify()

    back_btn = tk.Button(app.smooth_root, text="⬅ Back to Home", command=back_to_home)
    back_btn.pack(pady=8)    



#####################################################################################



#feature
def launch_fifth_gui(app):

    app.feature_root = tk.Toplevel(app.root)
    app.feature_root.title("Feature Extraction")
    #center_window(app.feature_root)                            


    flbm_global = None
    # Load and map data
    np_data = np.load("smooth.npz") 
    smsl = [list(map(tuple, np_data[key])) for key in np_data]

    # Point-to-line mapping
    point_to_lines = {}
    for idx, line in enumerate(smsl):
        for pt in line:
            if pt not in point_to_lines:
                point_to_lines[pt] = []
            point_to_lines[pt].append(idx)

    # GUI state
    selection_mode = tk.StringVar(value="")  # "top" or "bottom" or ""
    top_and_bottom_idx = [[], []]  # [top_idx, bottom_idx]

    # Initial plot
    fig, ax = plot_layers_GUI(smsl)
    canvas = FigureCanvasTkAgg(fig, master=app.feature_root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.grid(row=0, column=0, columnspan=3)

    # Redraw helper
    def update_plot_with_highlight(selected_idx):
        new_fig, new_ax = plot_layers_GUI(smsl, marked_layers=smsl[selected_idx])
        # Replace the entire canvas with a new figure
        canvas.figure = new_fig
        canvas._tkcanvas.destroy()
        new_canvas_widget = FigureCanvasTkAgg(new_fig, master=app.feature_root)
        new_canvas_widget.get_tk_widget().grid(row=0, column=0, columnspan=3)
        new_canvas_widget.draw()
        canvas.__dict__.update(new_canvas_widget.__dict__)  # Sync internal state
        canvas.mpl_connect("button_press_event", on_canvas_click)



    # Event: canvas click
    def on_canvas_click(event):
        if selection_mode.get() not in ["top", "bottom", "remove"]:
            return  # do nothing if no mode is active

        # Convert canvas to data coordinates
        inv = ax.transData.inverted()
        x, y = inv.transform((event.x, event.y))
        clicked = (round(x), round(y))

        # Find lines containing this point
        candidates = [idx for pt, indices in point_to_lines.items() if abs(pt[0] - clicked[0]) <= 1 and abs(pt[1] - clicked[1]) <= 1 for idx in indices]
        if not candidates:
            messagebox.showinfo("Invalid", "Please try again.")
            return

        # Choose the longest line
        selected_idx = max(set(candidates), key=lambda idx: len(smsl[idx]))

        # Update state and GUI
        if selection_mode.get() == "top":
            top_and_bottom_idx[0] = selected_idx
            top_btn.config(state="disabled")
            remove_btn.config(state="disabled")
        elif selection_mode.get() == "bottom":
            top_and_bottom_idx[1] = selected_idx
            bottom_btn.config(state="disabled")
            remove_btn.config(state="disabled")
        elif selection_mode.get() == "remove":
            update_plot_with_highlight(selected_idx)
            confirm = messagebox.askyesno("Confirm Deletion", f"Do you want to delete layer {selected_idx}?")
            if confirm:
                smsl[selected_idx] = []
            

        update_plot_with_highlight(selected_idx)
        selection_mode.set("")

        if all(isinstance(i, int) for i in top_and_bottom_idx):
            export_btn.config(state="normal")

    canvas.mpl_connect("button_press_event", on_canvas_click)

    # Event handlers for buttons
    def activate_top_select():
        selection_mode.set("top")

    def activate_bottom_select():
        selection_mode.set("bottom")
        
    def activate_remove_select():
        selection_mode.set("remove")        
        
    def save_final_figure():
        try:
            width = int(width_entry.get())
            height = int(height_entry.get())
        except ValueError:
            messagebox.showerror("Invalid Input", "Width and height must be integers.")
            return

        dpi = 100
        fig_size = (width / dpi, height / dpi)
        fig = Figure(figsize=fig_size, dpi=dpi)
        ax = fig.add_subplot(111)

        if not flbm_global:
            messagebox.showinfo("No Feature", "No feature to save yet.")
            return

        # Plot flbm without marked layers
        all_x = []
        all_y = []

        def process(layer):
            unique = list(set(layer))
            return sorted(unique, key=lambda pt: pt[0])

        for line in flbm_global:
            cleaned = process(line)
            if cleaned:
                x = [pt[0] for pt in cleaned]
                y = [pt[1] for pt in cleaned]
                all_x.extend(x)
                all_y.extend(y)
                ax.plot(x, y, c='black', linewidth=1)

        if all_x and all_y:
            ax.set_xlim(min(all_x), max(all_x))
            ax.set_ylim(min(all_y), max(all_y))

        ax.axis('off')
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
        fig.savefig("feature.png", bbox_inches='tight', pad_inches=0, dpi=dpi)
        messagebox.showinfo("Saved", "Saved as feature.png")
        

        
    def export_result():
        nonlocal flbm_global
        print("Top and Bottom Layer Indices:", top_and_bottom_idx)
        messagebox.showinfo("Exported", f"Top: {top_and_bottom_idx[0]}, Bottom: {top_and_bottom_idx[1]}")

        # === Step 1: Extract and filter ===
        maxdownlaps = [smsl[idx] for idx in top_and_bottom_idx]
        flbm = filter_lines_between_maxdownlaps(smsl, maxdownlaps, threshold=0.90)
        flbm_global = flbm
        print(f"Extracted feature {len(flbm)}")

        # === Step 2: Save result ===
        np_data = [np.array(line) for line in flbm]
        np.savez("feature.npz", *np_data)

        # === Step 3: Update canvas with new data ===
        new_fig, new_ax = plot_layers_GUI(flbm)
        canvas.figure = new_fig
        canvas._tkcanvas.destroy()
        new_canvas_widget = FigureCanvasTkAgg(new_fig, master=app.feature_root)
        new_canvas_widget.get_tk_widget().grid(row=0, column=0, columnspan=3)
        new_canvas_widget.draw()
        canvas.__dict__.update(new_canvas_widget.__dict__)
        canvas.mpl_connect("button_press_event", on_canvas_click)
        save_btn.config(state="normal")  # Enable save button
        save_btn.config(command=save_final_figure)



    # Buttons
    top_btn = tk.Button(app.feature_root, text="Select Top Layer", command=activate_top_select)
    top_btn.grid(row=1, column=0)

    bottom_btn = tk.Button(app.feature_root, text="Select Bottom Layer", command=activate_bottom_select)
    bottom_btn.grid(row=1, column=1)
    
    remove_btn = tk.Button(app.feature_root, text="Select Removing Layer", command=activate_remove_select)
    remove_btn.grid(row=1, column=3)    

    export_btn = tk.Button(app.feature_root, text="Export Top and Bottom", command=export_result, state="disabled")
    export_btn.grid(row=1, column=2)
    
    tk.Label(app.feature_root, text="Width:").grid(row=2, column=0)
    width_entry = tk.Entry(app.feature_root)
    width_entry.insert(0, "512")
    width_entry.grid(row=2, column=1)

    tk.Label(app.feature_root, text="Height:").grid(row=3, column=0)
    height_entry = tk.Entry(app.feature_root)
    height_entry.insert(0, "512")
    height_entry.grid(row=3, column=1)

    save_btn = tk.Button(app.feature_root, text="Save Feature Image", state="disabled")
    save_btn.grid(row=4, column=0, columnspan=3)
    
    # --- Back to Home handling ---
    def back_to_home():
        app.feature_root.destroy()
        app.root.deiconify()

    back_btn = tk.Button(app.feature_root, text="⬅ Back to Home", command=back_to_home)
    back_btn.grid(row=4, column=1, columnspan=3)
    
    
#Layer
def launch_Layer_selector(app):
    """
    Launch a GUI to visually select lines from `smsl` and build an Layer list.

    Parameters
    ----------
    smsl : list of lists of (x,y)
        Input lines to display.
    Layer : list, optional
        If provided, selections are appended to this list. If None, an internal list is used.
    select_tolerance : int, optional
        Pixel-distance tolerance when matching a clicked coordinate to points in smsl (in data coords).
    """
    app.layer_root = tk.Toplevel(app.root)
    app.layer_root.title("Layer Selection")    
    #center_window(app.layer_root)                            
    
    np_data = np.load("feature.npz") 
    smsl = [list(map(tuple, np_data[key])) for key in np_data]
    smsl = sorted(smsl, key=lambda l: l[0][1] if l else 0)
    
    Layer=None
    select_tolerance=2
    
    if Layer is None:
        Layer = []

        
    # Build quick point->line index mapping for fast lookup
    point_to_lines = {}
    for idx, line in enumerate(smsl):
        for pt in line:
            point_to_lines.setdefault(pt, []).append(idx)

    # GUI state
    selection_mode = tk.BooleanVar(value=False)  # True while waiting for a canvas click to select
    last_selected_idx = [None]  # store last selected index (list so closures can modify it)
    selected_stack = []  # keep a stack of selected indices (for undo/visualization)

    fig, ax = plot_layers_GUI(smsl)  # user-provided plotting function
    canvas = FigureCanvasTkAgg(fig, master=app.layer_root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.grid(row=0, column=0, columnspan=4, padx=6, pady=6)

    # Helper to re-draw the whole canvas and overlay selected lines from selected_stack in red
    def refresh_canvas(highlight_idx=None):
        """Redraw base plot and overlay all selected lines in red; optional highlight_idx is emphasized."""
        new_fig, new_ax = plot_layers_GUI(smsl)  # re-plot base data
        
        # Overlay previously selected lines
        for i, sidx in enumerate(selected_stack):
            if sidx is None:
                continue
            line = smsl[sidx]
            if not line:
                continue
            xs = [p[0] for p in line]
            ys = [p[1] for p in line]
            # last selected gets thicker linewidth
            lw = 2.5 if i == len(selected_stack) - 1 else 1.5
            new_ax.plot(xs, ys, linewidth=lw, linestyle='-', marker=None, zorder=5, label=f"sel_{sidx}", color='red')

        # If an extra index to highlight (e.g., immediate selection), ensure it's plotted more visible
        if highlight_idx is not None and highlight_idx not in selected_stack:
            line = smsl[highlight_idx]
            if line:
                xs = [p[0] for p in line]
                ys = [p[1] for p in line]
                new_ax.plot(xs, ys, linewidth=3.0, linestyle='-', marker=None, zorder=6, color='darkred')

        # Replace the canvas widget with new one
        nonlocal canvas
        try:
            canvas._tkcanvas.destroy()
        except Exception as e:
            pass
        canvas = FigureCanvasTkAgg(new_fig, master=app.layer_root)
        canvas.get_tk_widget().grid(row=0, column=0, columnspan=4, padx=6, pady=6)
        canvas.draw()
        # re-bind click handler to new canvas
        canvas.mpl_connect("button_press_event", on_canvas_click)

    # Utility: find candidate line indices near (x,y) in data coords
    def find_candidate_indices(x, y, tol=select_tolerance):
        # We'll look for any point in any line within tolerance (in data coords)
        candidates = set()
        # quick coarse search: iterate through lines (ok for modest number of lines)
        for idx, line in enumerate(smsl):
            for px, py in line:
                if abs(px - x) <= tol and abs(py - y) <= tol:
                    candidates.add(idx)
                    break
        return list(candidates)

    # Canvas click handler
    def on_canvas_click(event):
        # Only active while selection mode True
        if not selection_mode.get():
            return
        # event.xdata, event.ydata are None if click outside axes
        if event.xdata is None or event.ydata is None:
            messagebox.showinfo("Out of axes", "Click inside the plot area.")
            selection_mode.set(False)
            return

        x_click = round(event.xdata)
        y_click = round(event.ydata)

        candidates = find_candidate_indices(x_click, y_click, tol=select_tolerance)
        if not candidates:
            messagebox.showinfo("No layer found", "No layer found at the clicked location. Try again.")
            selection_mode.set(False)
            return

        # pick the longest line among candidates (heuristic)
        selected_idx = max(candidates, key=lambda ii: len(smsl[ii]))
        last_selected_idx[0] = selected_idx
        selected_stack.append(selected_idx)

        # update widget states
        undo_btn.config(state="normal")
        save_btn.config(state="normal")

        # redraw highlighting selected
        refresh_canvas(highlight_idx=selected_idx)

        # exit selection mode automatically
        selection_mode.set(False)

    # Button callbacks
    def activate_select_mode():
        selection_mode.set(True)
#         messagebox.showinfo("Select Mode", "Click on a layer in the canvas to select it.")  # so annoying when it is on

    def undo_selection():
        if not selected_stack:
            messagebox.showinfo("Nothing to undo", "No selection to undo.")
            return
        removed = selected_stack.pop()
        # If we removed the last selected and it was the one not saved yet, clear last_selected_idx
        if last_selected_idx[0] == removed:
            last_selected_idx[0] = None
            save_btn.config(state="disabled")
        # If stack empty then disable undo
        if not selected_stack:
            undo_btn.config(state="disabled")
        refresh_canvas()

    def save_selection():
        if last_selected_idx[0] is None:
            messagebox.showinfo("No selection", "No selected layer to save. Click 'Select Layer' first.")
            return
        idx = last_selected_idx[0]
        Layer.append(list(smsl[idx]))  # append a copy (list of tuples)  # old
        messagebox.showinfo("Saved", f"Selected layer {idx} appended to Layer (now {len(Layer)} items).")
        # after save, user may choose to continue selecting; disable save button until another selection
        save_btn.config(state="disabled")
        # keep last_selected_idx but allow further selections
        last_selected_idx[0] = None

    def export_Layer_npz():
        if not Layer:
            if not messagebox.askyesno("Empty", "Layer is empty. Export empty file?"):
                return
        # convert each line to numpy array and save as separate arrays inside single npz
        np_arrays = [np.array(line) for line in Layer]
        # name arrays arr0, arr1, ...
        kwargs = {f"arr{i}": arr for i, arr in enumerate(np_arrays)}
        np.savez("layer.npz", **kwargs)
        messagebox.showinfo("Exported", f"Saved {len(np_arrays)} lines to layer.npz")

    # Buttons and layout
    select_btn = tk.Button(app.layer_root, text="Select Layer", command=activate_select_mode)
    select_btn.grid(row=1, column=0, padx=6, pady=6)

    undo_btn = tk.Button(app.layer_root, text="Undo Selection", command=undo_selection, state="disabled")
    undo_btn.grid(row=1, column=1, padx=6, pady=6)

    save_btn = tk.Button(app.layer_root, text="Save Selection", command=save_selection, state="disabled")
    save_btn.grid(row=1, column=2, padx=6, pady=6)

    export_btn = tk.Button(app.layer_root, text="Export layer.npz", command=export_Layer_npz)
    export_btn.grid(row=1, column=3, padx=6, pady=6)


    # --- Back to Home handling ---
    def back_to_home():
        app.layer_root.destroy()
        app.root.deiconify()

    back_btn = tk.Button(app.layer_root, text="⬅ Back to Home", command=back_to_home)  # fixed by Claude, good job
    back_btn.grid(row=5, column=0, columnspan=3, pady=8)

    # initial canvas connection
    canvas.mpl_connect("button_press_event", on_canvas_click)

    # Bind window close to returning Layer (if needed by caller)
    def on_close():
        # Optionally ask for confirmation if there are unsaved picks
        if last_selected_idx[0] is not None:
            if not messagebox.askyesno("Exit", "You have a selected layer that wasn't saved — exit anyway?"):
                return
        app.layer_root.destroy()

    app.layer_root.protocol("WM_DELETE_WINDOW", on_close)

    # show initial plot
    refresh_canvas()


    # When GUI closes, return Layer so caller can use them
    return Layer
    

    
#####################################################################################
##################~~~~~~~~~~~~~~~~~~~~~~~~~ 

#manual

# matplotlib.use('TkAgg') # no need for this and later: %matplotlib inline   import matplotlib.pyplot as plt
    
def spatial_entropy(mask, radius=30):
    """Return mean local entropy over mask."""
    mask_u8 = img_as_ubyte(mask.astype(np.uint8))
    ent = entropy(mask_u8, disk(radius))
    return ent.mean()


def launch_seismic_mask_gui(app):
    app.manual_root = tk.Toplevel(app.root)
    app.manual_root.title("Manual Interpretation")
    #center_window(app.manual_root)          
    
    ### old __init__
    # Load cropped array
    try:
        data_npz = np.load('cropped.npz')
        app.arr = data_npz['data']
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load cropped.npz\n{e}")
        app.manual_root.destroy()
        return
    app.original_shape = app.arr.shape  # (H, W)
    print("test shape:", app.original_shape)    
    app.display_size = (512, 512)  # (width, height)

    # Display resize for drawing (note: PIL resize takes size=(width, height))
    app.display_arr = np.array(Image.fromarray(app.arr).resize(app.display_size[::-1], Image.BILINEAR))

    # Initialize storage
    app.layer_masks = []            # list of masks in original array shape
    app.clicks = []                 # clicks for current click-session (display coords)
    app.current_mode = 'Layer drawing'
    app.num_layers = 1

    # For "Layer drawing" mode: accumulate merged line data into this mask (original shape)
    app.current_line_mask = np.zeros(app.original_shape, dtype=np.uint8)
    # index of the layer_masks entry that holds the current_line_mask while user is merging;
    # None means not yet created in layer_masks
    app.current_saving_layer_idx = None
    
    app.layer_lines = []   # vector form (bresenham points)
    app.current_layer_lines = []  # accumulates during one layer session    

    # --------------------------------------------------------
    def _build_controls():
        frm = ttk.Frame(app.manual_root)
        frm.pack(side=tk.TOP, fill=tk.X, pady=5)

        # Number of layers
        ttk.Label(frm, text='Number of layers:').grid(row=0, column=0, sticky=tk.W)
        app.spin_layers = ttk.Spinbox(frm, from_=1, to=20, width=5)
        app.spin_layers.set(1)
        app.spin_layers.grid(row=0, column=1, padx=6)

        # Drawing mode combo (added "Layer drawing")
        ttk.Label(frm, text='Drawing mode:').grid(row=1, column=0, sticky=tk.W)
        app.combo_mode = ttk.Combobox(frm, values=['Layer drawing', 'Mask drawing'],
                                       state='readonly', width=15)
        app.combo_mode.current(0)
        app.combo_mode.bind('<<ComboboxSelected>>', _on_mode_change)
        app.combo_mode.grid(row=1, column=1, sticky='we', padx=6)

        # Layer view combo
        ttk.Label(frm, text='Layer view:').grid(row=2, column=0, sticky=tk.W)
        app.combo_layer = ttk.Combobox(frm, values=['All layers'], state='readonly', width=15)
        app.combo_layer.current(0)
        app.combo_layer.bind('<<ComboboxSelected>>', _on_layer_view_change)
        app.combo_layer.grid(row=2, column=1, sticky='we', padx=6)


        # Buttons
        # app.btn_start = ttk.Button(frm, text='Start', command=_on_start)
        # app.btn_start.grid(row=0, column=2, padx=10)
        
        app.btn_merge = ttk.Button(frm, text='Merge', command=_on_merge, state='disabled')
        app.btn_merge.grid(row=0, column=3, padx=10)

        app.btn_undo = ttk.Button(frm, text='Undo', command=_on_undo)
        app.btn_undo.grid(row=0, column=4, padx=10)
        
        app.btn_next = ttk.Button(frm, text='Finish', command=_on_next, state='disabled')  # edited Next>Finish (reason explained _on_next)
        app.btn_next.grid(row=0, column=5, padx=10)        
        
        
        # --- Back to Home handling ---
        def back_to_home():
            app.manual_root.destroy()
            app.root.deiconify()

        back_btn = ttk.Button(frm, text="⬅ Back to Home", command=back_to_home)
        back_btn.grid(row=1, column=5, padx=10)        

        # app.btn_finish = ttk.Button(frm, text='Finish', command=_on_finish, state='disabled')
        # app.btn_finish.grid(row=0, column=6, padx=10)

    # --------------------------------------------------------
    def _build_canvas():
        app.fig = plt.Figure()
        app.ax = app.fig.add_subplot(111)
        app.canvas = FigureCanvasTkAgg(app.fig, master=app.manual_root)
        app.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        app.cid_click = app.canvas.mpl_connect('button_press_event', _on_click)

    # --------------------------------------------------------
    def _refresh_display(img, overlay=None):
        """Refresh the seismic display with RdBu colormap centered on 0 amplitude."""
        app.ax.clear()

        # Dynamic symmetric range for amplitude
        vmax = np.percentile(img, 99)
        vmin = -vmax

        app.ax.imshow(
            img,
            cmap='RdBu',
            vmin=vmin,
            vmax=vmax,
            origin='upper',
            aspect='auto'
        )

        if overlay is not None:
            # overlay is expected to be in display size (512x512) and in [0,1]
            app.ax.imshow(overlay, cmap='plasma', alpha=0.8, origin='upper', aspect='auto')

        app.ax.set_xticks([])
        app.ax.set_yticks([])
        app.canvas.draw()

    # --------------------------------------------------------
    def _on_mode_change(event=None):
        app.current_mode = app.combo_mode.get()

        # adjust Merge button availability
        if app.current_mode == 'Layer drawing':
            # if there are clicks ready, enable Merge; otherwise keep disabled until click(s)
            messagebox.showinfo('Info', 'Clicks record pixels marking a layer. Pressing "Merge" converts selected points to a line.')
            app.btn_merge.config(state='normal' if len(app.clicks) > 0 else 'disabled')
        else:
            messagebox.showinfo('Info', 'Click on any number of points for choosing the mask with min and max amplitudes of selected points.')
            app.btn_merge.config(state='disabled')

    # --------------------------------------------------------
    def _on_layer_view_change(event=None):
        val = app.combo_layer.get()
        if val == 'All layers' or not app.layer_masks:
            _refresh_display(app.display_arr)
        else:
            idx = int(val.split()[-1]) - 1
            overlay = _resize_mask(app.layer_masks[idx])
            _refresh_display(app.display_arr, overlay=overlay)

    # --------------------------------------------------------
    def _on_start():
        _current_saving_layer_idx = None
        app.num_layers = int(app.spin_layers.get())
        # default to Mask drawing as previous behavior, user can change to Layer drawing
        app.combo_mode.set('Layer drawing')
        app.current_mode = 'Layer drawing'
        app.clicks = []
        app.current_line_mask = np.zeros(app.original_shape, dtype=np.uint8)
        app.btn_next.config(state='disabled')
        app.btn_merge.config(state='disabled')

    # --------------------------------------------------------
    def _on_click(event):
        # ignore clicks outside axes
        if event.xdata is None or event.ydata is None:
            return

        x, y = int(event.xdata), int(event.ydata)

        # MODE: Mask drawing (original behavior)
        if app.current_mode == 'Mask drawing':
            amp = app.display_arr[y, x]
            app.clicks.append((x, y, amp))

            # Enable "Next" if at least one click was made
            if len(app.clicks) > 0:
                app.btn_next.config(state='normal')

            # Visual feedback of clicked points
            app.ax.plot(x, y, 'yo', markersize=4)
            app.canvas.draw()
            return

        # MODE: Layer drawing (new behavior)
        if app.current_mode == 'Layer drawing':
            # record pure display pixel coords (no amplitude)
            app.clicks.append((x, y))

            # Visual feedback
            app.ax.plot(x, y, 'co', markersize=4)
            app.canvas.draw()

            # enable Merge button once there is at least one click
            app.btn_merge.config(state='normal')

            # also allow Next if we already have some merged content in current_line_mask
            if app.current_line_mask.any():
                app.btn_next.config(state='normal')

            return

        # other modes: ignore
        return

    # --------------------------------------------------------
    def _on_merge():
        """Merge current click-session into the current_line_mask using bresenham_line.

        Behavior:
        - Sort points by x-value
        - Call bresenham_line(points) which is expected to return a list of corrected (x,y) points (display coords)
        - Map those display coords back to original arr coords and set mask pixels to 1
        - Clear the click-session but keep current_line_mask so multiple merges accumulate
        - Update preview overlay
        """
        if app.current_mode != 'Layer drawing':
            return
        if not app.clicks:
            return

        # Sort clicks by x (display coords)
        pts = sorted(app.clicks, key=lambda p: p[0])

        disp_w, disp_h = app.display_size
        orig_h, orig_w = app.original_shape

        # --- Map display coords to original coords FIRST ---
        orig_pts = []
        for (dx, dy) in pts:
            ox = int(np.clip(round(dx * (orig_w / disp_w)), 0, orig_w - 1))
            oy = int(np.clip(round(dy * (orig_h / disp_h)), 0, orig_h - 1))
            orig_pts.append((ox, oy))

        # --- Apply Bresenham in original coordinate space ---
        corrected = bresenham_line(orig_pts)

        # store vector representation
        app.current_layer_lines.append(corrected)


        # --- Fill mask directly in original coordinates ---
        for (ox, oy) in corrected:
            app.current_line_mask[oy, ox] = 1


        # ---- BEGIN REPLACEMENT BLOCK ----
        # Clear clicks for next session but keep merged mask content
        app.clicks = []
        app.btn_merge.config(state='disabled')

        # If there is already a layer entry for this ongoing drawing session, update it;
        # otherwise create a new layer entry so Layer view shows the merged result immediately.
        if app.current_saving_layer_idx is None:
            # create new visible layer (placeholder) and record its index
            app.layer_masks.append(app.current_line_mask.copy().astype(np.uint8))
            app.current_saving_layer_idx = len(app.layer_masks) - 1
        else:
            # update existing saved layer so multiple merges accumulate into the same layer
            app.layer_masks[app.current_saving_layer_idx] = app.current_line_mask.copy().astype(np.uint8)

        # Update layer dropdown values and select the current layer
        vals = ['All layers'] + [f'Layer {i+1}' for i in range(len(app.layer_masks))]
        app.combo_layer['values'] = vals
        app.combo_layer.set(f'Layer {app.current_saving_layer_idx + 1}')

        # After each merge, user can still Next to finalize and move on
        app.btn_next.config(state='normal')

        # Update overlay shown on display (resize to display size)
        overlay = _resize_mask(app.current_line_mask)
        _refresh_display(app.display_arr, overlay=overlay)
        # ---- END REPLACEMENT BLOCK ----
                    
        # Clear clicks for next session but keep merged mask
        app.clicks = []
        app.btn_merge.config(state='disabled')

        # After each merge, user can still Next to save this mask
        app.btn_next.config(state='normal')

        # Update overlay shown on display (resize to display size)
        overlay = _resize_mask(app.current_line_mask)
        _refresh_display(app.display_arr, overlay=overlay)
     
        
    # --------------------------------------------------------            
    def _on_next():
        """Called when user finishes selecting clicks for current layer.

        Behavior differs by current mode:
        - Mask drawing: original amplitude-range based mask creation
        - Layer drawing: save current_line_mask (accumulated merges) as a layer
        """
        if app.current_mode == 'Mask drawing':
            if not app.clicks:
                return

            amps = [a for (_, _, a) in app.clicks]
            lower, upper = min(amps), max(amps)

            # Inclusive both bounds this time
            mask = (app.arr >= lower) & (app.arr <= upper)
            
            # --------------------------------------------------
            # ADD TOP & BOTTOM FRAME LAYERS (Bresenham-based)
            # --------------------------------------------------
            orig_h, orig_w = app.original_shape

            # --- Top frame line ---
            top_pts = [(0, 0), (orig_w - 1, 0)]
            top_line = bresenham_line(top_pts)

            top_mask = np.zeros(app.original_shape, dtype=np.uint8)
            for (ox, oy) in top_line:
                top_mask[oy, ox] = 1


            # --- Bottom frame line ---
            bottom_pts = [(0, orig_h - 1), (orig_w - 1, orig_h - 1)]
            bottom_line = bresenham_line(bottom_pts)

            bottom_mask = np.zeros(app.original_shape, dtype=np.uint8)
            for (ox, oy) in bottom_line:
                bottom_mask[oy, ox] = 1

            
            
            final_mask = mask.astype(np.uint8) | top_mask | bottom_mask
            app.layer_masks.append(final_mask)


            # Reset for next layer
            app.clicks = []
            app.btn_next.config(state='disabled')

            messagebox.showinfo(
                'Mask created',
                f'Mask {len(app.layer_masks)} created.\nAmplitude range: {lower:.2f} to {upper:.2f}'
            )
            
            # # Update layer dropdown
            # vals = ['All layers'] + [f'Layer {i+1}' for i in range(len(app.layer_masks))]
            # app.combo_layer['values'] = vals
            # if app.layer_masks:
                # app.combo_layer.set(f'Layer {len(app.layer_masks)}')

            # Show overlay of the new layer if exists
            if app.layer_masks:
                overlay = _resize_mask(app.layer_masks[-1])
                _refresh_display(app.display_arr, overlay=overlay)

            # if len(app.layer_masks) >= app.num_layers:
                # app.btn_finish.config(state='normal')
                # app.combo_mode.set('Layer drawing')
                # app.current_mode = 'Layer drawing'


            ### added _on_finish to _on_next to forbid creating extra masks/layers
            # if user has drawn a current line mask but didn't press Next, save it as final layer
            if app.current_line_mask.any():
                app.layer_masks.append(app.current_line_mask.copy().astype(np.uint8))

            if not app.layer_masks:
                messagebox.showinfo('No layers', 'No layers to save.')
                return


            masks = np.stack(app.layer_masks)

            np.savez_compressed('masks.npz', masks=masks[0])
            messagebox.showinfo('Saved', f'Saved results to masks.npz.')
            

        elif app.current_mode == 'Layer drawing':
            # If there are any pending clicks that were not merged, attempt an implicit merge.
            if app.clicks:
                try:
                    _on_merge()
                except Exception as e:
                    messagebox.showerror('Merge failed', f'Merge failed: {e}')
                    return

            # Now ensure there's something to save
            if not app.current_line_mask.any():
                messagebox.showinfo('No data', 'No drawn lines to save for this layer.')
                return

            # If merge already created/updated the visible layer, keep it; otherwise append now.
            # --------------------------------------------------
            # ADD TOP & BOTTOM FRAME LINES INTO CURRENT LAYER
            # --------------------------------------------------
            orig_h, orig_w = app.original_shape

            # Top frame
            top_pts = [(0, 0), (orig_w - 1, 0)]
            top_line = bresenham_line(top_pts)
            for (ox, oy) in top_line:
                app.current_line_mask[oy, ox] = 1

            # Bottom frame
            bottom_pts = [(0, orig_h - 1), (orig_w - 1, orig_h - 1)]
            bottom_line = bresenham_line(bottom_pts)
            for (ox, oy) in bottom_line:
                app.current_line_mask[oy, ox] = 1
            
            
            if app.current_saving_layer_idx is None:
                # (rare) case: no explicit merge happened but current_line_mask exists -> append it
                app.layer_masks.append(app.current_line_mask.copy().astype(np.uint8))
            else:
                # ensure the stored layer matches current_line_mask (update)
                app.layer_masks[app.current_saving_layer_idx] = app.current_line_mask.copy().astype(np.uint8)
                
            # for layers export
            app.current_layer_lines.append(top_line)  # frame line
            app.current_layer_lines.append(bottom_line) # frame line
            app.layer_lines = app.current_layer_lines                

            # Reset working buffers for next layer session
            app.current_line_mask = np.zeros(app.original_shape, dtype=np.uint8)
            app.current_saving_layer_idx = None
            app.clicks = []
            app.btn_next.config(state='disabled')
            # for bresenham layers
            app.current_layer_lines = []

            messagebox.showinfo('Layer created', f'Layer {len(app.layer_masks)} created (line drawing).')
            
            # # Update layer dropdown
            # vals = ['All layers'] + [f'Layer {i+1}' for i in range(len(app.layer_masks))]
            # app.combo_layer['values'] = vals
            # if app.layer_masks:
                # app.combo_layer.set(f'Layer {len(app.layer_masks)}')

            # Show overlay of the new layer if exists
            if app.layer_masks:
                overlay = _resize_mask(app.layer_masks[-1])
                _refresh_display(app.display_arr, overlay=overlay)

            # if len(app.layer_masks) >= app.num_layers:
                # app.btn_finish.config(state='normal')
                # app.combo_mode.set('Layer drawing')
                # app.current_mode = 'Layer drawing'


            ### added _on_finish to _on_next to forbid creating extra masks/layers
            # if user has drawn a current line mask but didn't press Next, save it as final layer
            if app.current_line_mask.any():
                app.layer_masks.append(app.current_line_mask.copy().astype(np.uint8))

            if not app.layer_masks:
                messagebox.showinfo('No layers', 'No layers to save.')
                return


            masks = np.stack(app.layer_masks)

            
            np_data = [np.array(line) for line in app.layer_lines]
            np.savez(f"layers.npz", *np_data)    
            messagebox.showinfo('Saved', f'Saved results to layers.npz.')
            

            

        app.manual_root.destroy()
        # Start contours GUI (user's function)
        launch_second_gui(app)

    # --------------------------------------------------------
    def _on_undo():
            # --- MODE: Layer drawing Undo ---
            if app.current_mode == 'Layer drawing':
                if app.current_layer_lines:
                    # 1. Remove the very last segment (the last merge)
                    app.current_layer_lines.pop()
                    
                    # 2. Re-build the mask from scratch based on remaining lines
                    app.current_line_mask[:] = 0
                    for line_pts in app.current_layer_lines:
                        for (ox, oy) in line_pts:
                            app.current_line_mask[oy, ox] = 1
                    
                    # 3. Update the persistent storage (layer_masks)
                    if app.current_saving_layer_idx is not None:
                        app.layer_masks[app.current_saving_layer_idx] = app.current_line_mask.copy().astype(np.uint8)
                    
                    # 4. Refresh display
                    overlay = _resize_mask(app.current_line_mask)
                    _refresh_display(app.display_arr, overlay=overlay)
                    
                    # 5. UI Housekeeping
                    if not app.current_layer_lines:
                        app.btn_next.config(state='disabled')
                    
                    return # Exit early
                else:
                    # If no lines are left in the current session, perhaps clear the clicks
                    app.clicks = []
                    app.btn_merge.config(state='disabled')
                    _refresh_display(app.display_arr)
                    messagebox.showinfo('Undo', 'No more segments to undo in this layer.')

            # --- MODE: Mask drawing (Original Behavior) ---
            else:
                app.current_saving_layer_idx = None
                app.current_line_mask[:] = 0 
                app.clicks = []
                app.btn_merge.config(state='disabled')
                app.btn_next.config(state='disabled')
                _refresh_display(app.display_arr)
                messagebox.showinfo('Undo', 'Cleared current drawing session.')
                
                
    # --------------------------------------------------------
    def _on_finish():
        # if user has drawn a current line mask but didn't press Next, save it as final layer
        if app.current_line_mask.any():
            app.layer_masks.append(app.current_line_mask.copy().astype(np.uint8))

        if not app.layer_masks:
            messagebox.showinfo('No layers', 'No layers to save.')
            return

        masks = np.stack(app.layer_masks)

        np.savez_compressed('masks.npz', masks=masks)
        messagebox.showinfo('Saved', f'Saved {len(app.layer_masks)} layers to masks.npz.')

        app.manual_root.destroy()
        # Start contours GUI (user's function)
        launch_second_gui(app)

    # --------------------------------------------------------
    def _resize_mask(mask):
        """Resize a mask (original shape) to 512x512 for overlay display."""
        img = Image.fromarray((mask * 255).astype(np.uint8))
        return np.array(img.resize(app.display_size[::-1], Image.NEAREST)) / 255.0

    # Load cropped array
    try:
        data_npz = np.load('cropped.npz')
        app.arr = data_npz['data']
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load cropped.npz\n{e}")
        app.manual_root.destroy()
        return



    _build_controls()
    _build_canvas()
    _refresh_display(app.display_arr)







    
###################################~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
#polarity mask check GUI class
###################################%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def launch_seismic_polarity_gui(app):
    # Prevent duplicate window
    if hasattr(app, "polarity_root") and app.polarity_root.winfo_exists():
        app.polarity_root.lift()
        return

    app.polarity_root = tk.Toplevel(app.root)
    app.polarity_root.title("Polarity")
    #center_window(app.polarity_root)        

    # --------------------------------------------------
    # Load cropped array
    try:
        data_npz = np.load('cropped.npz')
        app.arr = data_npz['data']
        if app.arr.ndim != 2:
            raise ValueError(f"cropped array must be 2D (got shape {app.arr.shape}).")
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load cropped.npz\n{e}")
        app.polarity_root.destroy()
        return

    # --------------------------------------------------
    # Store original and display arrays
    app.original_shape = app.arr.shape
    app.display_size = (512, 512)  # (width, height)

    pil = Image.fromarray(app.arr)
    pil_resized = pil.resize(app.display_size[::-1], Image.BILINEAR)
    app.display_arr = np.array(pil_resized).astype(np.float32)

    # --------------------------------------------------
    def _compute_amp_ranges():
        arr = app.arr

        if np.any(arr > 0):
            app.pos_min = float(np.min(arr[arr > 0]))
            app.pos_max = float(np.max(arr))
        else:
            app.pos_min, app.pos_max = 0.0, 1.0

        if np.any(arr < 0):
            app.neg_min = float(np.min(arr))
            app.neg_max = float(np.max(arr[arr < 0]))
        else:
            app.neg_min, app.neg_max = -1.0, 0.0

    _compute_amp_ranges()

    # --------------------------------------------------
    # Initial bounds
    app.pos_val = max(
        0.0,
        np.percentile(app.arr[app.arr > 0] if np.any(app.arr > 0) else [0.0], 99)
    )
    app.neg_val = min(
        0.0,
        np.percentile(app.arr[app.arr < 0] if np.any(app.arr < 0) else [0.0], 1)
    )

    app.pos_var = tk.DoubleVar(value=app.pos_val)
    app.neg_var = tk.DoubleVar(value=app.neg_val)
    app.entry_pos_var = tk.StringVar(value=f"{app.pos_val:.3f}")
    app.entry_neg_var = tk.StringVar(value=f"{app.neg_val:.3f}")

    # --------------------------------------------------
    def _refresh_display(img, overlay=None):
        app.ax.clear()

        try:
            vmax = np.percentile(np.abs(app.display_arr), 99)
            if vmax <= 0:
                vmax = np.max(np.abs(app.display_arr)) or 1.0
        except Exception:
            vmax = np.max(np.abs(app.display_arr)) or 1.0

        app.ax.imshow(
            img, cmap='RdBu', vmin=-vmax, vmax=vmax,
            origin='upper', aspect='auto'
        )

        if overlay is not None:
            app.ax.imshow(
                overlay, cmap='plasma', alpha=0.9,
                origin='upper', aspect='auto'
            )

        app.ax.set_xticks([])
        app.ax.set_yticks([])
        app.canvas.draw()

    # --------------------------------------------------
    def _apply_and_refresh(pos_value, neg_value):
        mask_bin = ((app.arr >= neg_value) & (app.arr <= pos_value)).astype(np.uint8)

        overlay = Image.fromarray((mask_bin * 255).astype(np.uint8))
        overlay_resized = np.array(
            overlay.resize(app.display_size[::-1], Image.NEAREST)
        ) / 255.0

        _refresh_display(app.display_arr, overlay=overlay_resized)
        app.last_mask = mask_bin

    # --------------------------------------------------
    def _on_slider_change(event=None):
        pos_v = float(app.pos_var.get())
        neg_v = float(app.neg_var.get())

        app.lbl_pos_value.config(text=f"{pos_v:.3f}")
        app.lbl_neg_value.config(text=f"{neg_v:.3f}")

        app.entry_pos_var.set(f"{pos_v:.6g}")
        app.entry_neg_var.set(f"{neg_v:.6g}")

        _apply_and_refresh(pos_v, neg_v)

    # --------------------------------------------------
    def _on_finish():
        if not hasattr(app, 'last_mask'):
            _apply_and_refresh(app.pos_var.get(), app.neg_var.get())

        try:
            np.savez('masks.npz', masks=app.last_mask)
            messagebox.showinfo("Saved", "Filtered mask saved to masks.npz")
        except Exception as e:
            messagebox.showerror("Save error", str(e))

        app.polarity_root.destroy()
        launch_second_gui(app)

    # --------------------------------------------------
    # Controls
    frm = ttk.Frame(app.polarity_root)
    frm.pack(side=tk.TOP, fill=tk.X, padx=6, pady=6)

    ttk.Label(frm, text='Positive upper bound:').grid(row=0, column=0, sticky=tk.W)
    app.scale_pos = tk.Scale(
        frm, from_=app.pos_min, to=app.pos_max,
        resolution=max((app.pos_max - app.pos_min) / 1000.0, 1e-6),
        orient=tk.HORIZONTAL, variable=app.pos_var,
        command=_on_slider_change, length=300
    )
    app.scale_pos.grid(row=0, column=1, padx=6)
    app.lbl_pos_value = ttk.Label(frm, text=f"{app.pos_var.get():.3f}")
    app.lbl_pos_value.grid(row=0, column=2)

    ttk.Label(frm, text='Negative lower bound:').grid(row=1, column=0, sticky=tk.W)
    app.scale_neg = tk.Scale(
        frm, from_=app.neg_min, to=app.neg_max,
        resolution=max((app.neg_max - app.neg_min) / 1000.0, 1e-6),
        orient=tk.HORIZONTAL, variable=app.neg_var,
        command=_on_slider_change, length=300
    )
    app.scale_neg.grid(row=1, column=1, padx=6)
    app.lbl_neg_value = ttk.Label(frm, text=f"{app.neg_var.get():.3f}")
    app.lbl_neg_value.grid(row=1, column=2)

    ttk.Button(frm, text='Finish (save mask)', command=_on_finish).grid(
        row=0, column=3, rowspan=2, padx=(20, 0)
    )
    
    # --- Back to Home handling ---
    def back_to_home():
        app.polarity_root.destroy()
        app.root.deiconify()

    ttk.Button(frm, text='⬅ Back to Home', command=back_to_home).grid(
        row=2, column=3, rowspan=2, padx=(20, 0)
    )

    # --------------------------------------------------
    # Canvas
    app.fig = plt.Figure(figsize=(6, 6))
    app.ax = app.fig.add_subplot(111)
    app.canvas = FigureCanvasTkAgg(app.fig, master=app.polarity_root)
    app.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    _refresh_display(app.display_arr)

    
    
    
    
#plotting
def launch_plotting_gui(app):
    # reuse existing plotting window if possible
    try:
        if hasattr(app, "plotting_root") and app.plotting_root.winfo_exists():
            app.plotting_root.lift()
            return
    except Exception:
        pass

    # create as a Toplevel of main root to allow back-to-home
    # if app.root might be None, add safe fallback
    app.plotting_root = tk.Toplevel(app.root)
    app.plotting_root.title("Plotting")
    center_window(app.plotting_root, 900, 680)        

    # build your plotting GUI inside app.plotting_root ...
    # e.g., canvas, controls, etc.


    # left controls frame, right plot frame
    controls = tk.Frame(app.plotting_root)
    controls.pack(side=tk.LEFT, fill=tk.Y, padx=6, pady=6)

    plot_frame = tk.Frame(app.plotting_root)
    plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    # --- Variables ---
    app.plot_type_var = tk.StringVar(value="cropped")  # default
    app.axis_off_var = tk.BooleanVar(value=False)      # False => axis on
    app.svg_var = tk.BooleanVar(value=False)
    app.plot_width_var = tk.IntVar(value=512)
    app.plot_height_var = tk.IntVar(value=512)

    # NEW: overlay option and line color
    app.overlay_var = tk.StringVar(value="none")  # 'none' | 'original' | 'marked'
    app.line_color_var = tk.StringVar(value="lime")  # default color

    # radio group: file types
    tk.Label(controls, text="Source (npz):").pack(anchor='w')
    options = [("cropped.npz", "cropped"),
               ("contours.npz", "contours"),
               ("tsp.npz", "tsp"),
               ("onlap.npz", "onlap"),
               ("smooth.npz", "smooth"),
               ("feature.npz", "feature"),
               ("manual.npz", "manual")]
    for text, val in options:
        tk.Radiobutton(controls, text=text, variable=app.plot_type_var, value=val).pack(anchor='w')

    # overlay option radio group
    tk.Label(controls, text="Overlay option:").pack(anchor='w', pady=(8,0))
    tk.Radiobutton(controls, text="None", variable=app.overlay_var, value="none").pack(anchor='w')
    tk.Radiobutton(controls, text="Original", variable=app.overlay_var, value="original").pack(anchor='w')
    tk.Radiobutton(controls, text="Marked Layers", variable=app.overlay_var, value="marked").pack(anchor='w')

    # color combobox for line_color (common names)
    tk.Label(controls, text="Line color:").pack(anchor='w', pady=(8,0))
    # a compact list of useful color names; add/remove as you like
    color_choices = ['lime', 'red', 'blue', 'green', 'yellow', 'cyan', 'magenta', 'black', 'white', 'orange']
    color_combo = ttk.Combobox(controls, values=color_choices, textvariable=app.line_color_var, state="readonly", width=12)
    color_combo.pack(anchor='w')
    color_combo.set(app.line_color_var.get())

    # axis on/off
    tk.Label(controls, text="Axis:").pack(anchor='w', pady=(8,0))
    tk.Radiobutton(controls, text="On", variable=app.axis_off_var, value=False).pack(anchor='w')
    tk.Radiobutton(controls, text="Off", variable=app.axis_off_var, value=True).pack(anchor='w')

    # SVG file save on/off
    tk.Label(controls, text="Save SVG:").pack(anchor='w', pady=(8,0))
    tk.Radiobutton(controls, text="No", variable=app.svg_var, value=False).pack(anchor='w')
    tk.Radiobutton(controls, text="Yes", variable=app.svg_var, value=True).pack(anchor='w')    

    # dimensions
    tk.Label(controls, text="Dimensions (px):").pack(anchor='w', pady=(8,0))
    dim_frame = tk.Frame(controls)
    dim_frame.pack(anchor='w')
    tk.Label(dim_frame, text="W:").pack(side=tk.LEFT)
    tk.Entry(dim_frame, textvariable=app.plot_width_var, width=6).pack(side=tk.LEFT, padx=(2,6))
    tk.Label(dim_frame, text="H:").pack(side=tk.LEFT)
    tk.Entry(dim_frame, textvariable=app.plot_height_var, width=6).pack(side=tk.LEFT, padx=(2,0))

    # Plot & Save button
    def on_plot_and_save():
        # read controls
        plot_type = app.plot_type_var.get()
        overlay_opt = app.overlay_var.get()
        axis_off = bool(app.axis_off_var.get())
        svg_option = bool(app.svg_var.get())
        
        # if svg_option:
            # messagebox.showinfo("About saving .svg images", "Do not save .svg images on .exe version of this program because it will encounter an error.\n You can export plots as .svg images in the .py version.")
        
        width = int(app.plot_width_var.get())
        height = int(app.plot_height_var.get())
        pixels = (width, height)
        line_color = str(app.line_color_var.get())

        # load cropped.npz (required)
        try:
            cropped_npz = np.load("cropped.npz")
            if 'data' not in cropped_npz:
                messagebox.showerror("File error", "cropped.npz does not contain key 'data'.")
                return
            arr = (cropped_npz['data']).T   # shape (n_samples, n_traces)
        except FileNotFoundError:
            messagebox.showerror("File missing", "cropped.npz not found in current dir.")
            return
        except Exception as e:
            messagebox.showerror("Error loading cropped.npz", str(e))
            return

        # pull inline if present
        inline = int(cropped_npz["inline"][0]) if "inline" in cropped_npz else -999

        # required plotting params (from your object)
        x_start = getattr(app, "x_start", 0)
        y_start = getattr(app, "y_start", 0)
        crossline = getattr(app, "crossline_val", 0)
        twt_spacing = getattr(app, "twt_spacing", 4)
        initial_inline = getattr(app, "initial_inline", 0)
        #debug  
        print(x_start := getattr(app, "x_start", 0), y_start := getattr(app, "y_start", 0), crossline := getattr(app, "crossline_val", 0), twt_spacing := getattr(app, "twt_spacing", 4), initial_inline := getattr(app, "initial_inline", 0))
        try:
            # --- CASE: cropped.npz selected ---
            if plot_type == "cropped":
                if overlay_opt == "none":
                    # plain seismic
                    fig = plot_seismic_CROPPED(
                        seismic_profile=arr,
                        pixels=pixels,
                        marked_points=False,
                        axis_off=axis_off,
                        seismic_relief=99,
                        inline=(inline or 0),
                        twt_spacing=twt_spacing,
                        x_start=x_start,
                        y_start=y_start,
                        initial_inline=initial_inline,
                        crossline=crossline,
                        save_option=True,
                        close_option=True,
                        fixed_name="cropped",
                        return_fig=True,
                        svg_option=svg_option
                    )

                elif overlay_opt == "original":
                    # Seismic overlay on cropped does not make sense because cropped is the seismic profile itapp.
                    # Inform the user and show the plain cropped plot instead.
                    messagebox.showinfo("Overlay not available", "Original overlay is not applicable to cropped.npz — showing the plain cropped profile.")
                    fig = plot_seismic_CROPPED(
                        seismic_profile=arr,
                        pixels=pixels,
                        marked_points=False,
                        axis_off=axis_off,
                        seismic_relief=99,
                        inline=(inline or 0),
                        twt_spacing=twt_spacing,
                        x_start=x_start,
                        y_start=y_start,
                        initial_inline=initial_inline,
                        crossline=crossline,
                        save_option=True,
                        close_option=True,
                        fixed_name="cropped",
                        return_fig=True,
                        svg_option=svg_option
                    )

                elif overlay_opt == "marked":
                    # overlay marked layers on the cropped seismic (read layer.npz)
                    try:
                        np_data = np.load("layer.npz")
                        marked_layers = [list(map(tuple, np_data[key])) for key in np_data.files]
                        if not marked_layers:
                            raise ValueError("layer.npz contains no layers.")
                    except FileNotFoundError:
                        messagebox.showerror("File missing", "layer.npz not found in current dir.")
                        return
                    except Exception as e:
                        messagebox.showerror("Error loading layer.npz", str(e))
                        return

                    # map amplitude to the marked layers (overlay)
                    print(f"successfully reached map: {len(marked_layers)}")
                    fig = map_amplitude_to_lines(
                        seismic_profile=arr,
                        lines=marked_layers,
                        pixels=pixels,
                        seismic_relief=99,
                        inline=(inline or 0),
                        overlay_debug=False,
                        twt_spacing=twt_spacing,
                        axis_off=axis_off,
                        x_start=x_start,
                        y_start=y_start,
                        initial_inline=initial_inline,
                        crossline=crossline,
                        line_color=line_color,
                        return_fig=True,
                        close_option=False,
                        svg_option=svg_option
                    )
                else:
                    raise ValueError("Unknown overlay option")

            # --- CASE: any other .npz (line-based) ---
            else:
                fname = f"{plot_type}.npz"
                try:
                    npz = np.load(fname)
                except FileNotFoundError:
                    messagebox.showerror("File missing", f"{fname} not found.")
                    return
                except Exception as e:
                    messagebox.showerror("Error loading " + fname, str(e))
                    return

                # convert npz arrays into list of tuple lists (interpreted_layers)
                try:
                    interpreted_layers = []
                    for key in npz.files:
                        if key != "flag":
                            arr_layer = npz[key]
                            layer = [tuple(pt) for pt in arr_layer]
                            interpreted_layers.append(layer)
                except Exception as e:
                    messagebox.showerror("NPZ parse error", f"Couldn't parse {fname}: {e}")
                    return

                seismic_profile_shape = np.shape(arr)

                if overlay_opt == "none":
                    # existing behavior: plot layers only
                    fig = plot_layers(
                        layers=interpreted_layers,
                        pixels=pixels,
                        save_option=True,
                        marked_layers=False,
                        axis_off=axis_off,
                        close_option=False,
                        fixed_name=plot_type,
                        inline=inline,
                        initial_inline=initial_inline,
                        crossline=crossline,
                        x_start=x_start,
                        y_start=y_start,
                        twt_spacing=twt_spacing,
                        seismic_profile_shape=seismic_profile_shape,
                        return_fig=True,
                        svg_option=svg_option
                    )

                elif overlay_opt == "original":
                    # Use map_amplitude_to_lines: overlay the selected interpreted_layers onto seismic_profile
                    # (Seismic overlay only makes sense for non-cropped types)
                    fig = map_amplitude_to_lines(
                        seismic_profile=arr,
                        lines=interpreted_layers,
                        pixels=pixels,
                        seismic_relief=99,
                        inline=(inline or 0),
                        overlay_debug=False,
                        twt_spacing=twt_spacing,
                        axis_off=axis_off,
                        x_start=x_start,
                        y_start=y_start,
                        initial_inline=initial_inline,
                        crossline=crossline,
                        line_color=line_color,
                        return_fig=True,
                        close_option=False,
                        svg_option=svg_option
                    )

                elif overlay_opt == "marked":
                    # Read marked layers and pass them to plot_layers together with interpreted_layers
                    try:
                        np_data = np.load("layer.npz")
                        marked_layers = [list(map(tuple, np_data[key])) for key in np_data.files]
                        if not marked_layers:
                            raise ValueError("layer.npz contains no layers.")
                    except FileNotFoundError:
                        messagebox.showerror("File missing", "layer.npz not found in current dir.")
                        return
                    except Exception as e:
                        messagebox.showerror("Error loading layer.npz", str(e))
                        return

                    # plot_layers enhanced with line_color and marked_layers
                    fig = plot_layers(
                        layers=interpreted_layers,
                        pixels=pixels,
                        save_option=True,
                        marked_layers=marked_layers,
                        axis_off=axis_off,
                        close_option=False,
                        fixed_name=f"{plot_type}_marked",
                        inline=inline,
                        initial_inline=initial_inline,
                        crossline=crossline,
                        x_start=x_start,
                        y_start=y_start,
                        twt_spacing=twt_spacing,
                        seismic_profile_shape=seismic_profile_shape,
                        line_color=line_color,
                        return_fig=True,
                        svg_option=svg_option
                    )
                else:
                    raise ValueError("Unknown overlay option")

        except Exception as e:
            messagebox.showerror("Plot error", str(e))
            return

        # embed figure into Tk
        try:
            # clear previous canvas if exists
            if hasattr(app, "plot_canvas") and app.plot_canvas:
                try:
                    app.plot_canvas.get_tk_widget().pack_forget()
                    app.plot_canvas.get_tk_widget().destroy()
                except Exception as e:
                    pass

            app.plot_canvas = FigureCanvasTkAgg(fig, master=plot_frame)
            widget = app.plot_canvas.get_tk_widget()
            widget.pack(fill=tk.BOTH, expand=True)
            app.plot_canvas.draw()
        except Exception as e:
            messagebox.showerror("Canvas error", str(e))
            return


    # --- Back to Home handling ---
    def back_to_home():
        app.plotting_root.destroy()
        app.root.deiconify()

    tk.Button(controls, text="Plot and Save", command=on_plot_and_save).pack(pady=(12,0))
    tk.Button(controls, text="⬅ Back to Home", command=back_to_home).pack(pady=(12,0))

    # small info label
    tk.Label(controls, text="Note: cropped.npz is required for all plots.").pack(anchor='w', pady=(8,0))

    # make plot_frame scrollable / large
    plot_frame.config(width=512, height=512)
    plot_frame.pack_propagate(True)

    # keep reference to avoid garbage collection
    app.plotting_root = app.plotting_root






def center_window(window, width=None, height=None):
    """
    Centers a tkinter window. If width/height are None, 
    it calculates size based on widget packing.
    """
    if width is None or height is None:
        # Step A: Let packing decide the size
        window.update_idletasks()
        width = window.winfo_width()
        height = window.winfo_height()
    
    # Step B: Get screen dimensions
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()

    # Step C: Calculate center coordinates
    x = int((screen_width - width) / 2)
    y = int((screen_height - height) / 2)

    # Step D: Apply geometry
    window.geometry(f'{width}x{height}+{x}+{y}')
    
    

###################################%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        

# --- Main Entry Point ---
def main():
    root = tk.Tk()
    app = SeismicApp(root)

    def force_close():
        root.destroy()
        # Option A: Polite exit 
        sys.exit(0) 
        # Option B: Nuclear exit
        # os._exit(0) 

    root.protocol("WM_DELETE_WINDOW", force_close)
    root.mainloop()

if __name__ == "__main__":
    multiprocessing.freeze_support() # This line stops the subprocesses from re-running main() recursively
    

    main()
