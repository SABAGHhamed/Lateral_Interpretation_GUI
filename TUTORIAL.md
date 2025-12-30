# Quick tutorial — Interpret XL480 (example)

This short, hands-on tutorial shows how to interpret crossline **XL480** from the example subvolume
using the **Lateral Interpretation GUI**. It is intended for quick testing and demonstration;
for full methodological details please read the paper and the GUI *Help*.

**Prerequisites**
- You have the repository files on your system and have installed dependencies:
  ```bash
  python -m pip install -r requirements.txt
````

* Example SEGY (small demo) is available at:
  `example_data/F3_Dip_steered_median_subvolume_IL230-430_XL475-675_T1600-1800-CDP170-370.sgy`
* Tutorial images/screenshots are in `tutorial_images/`

---

## Overview of this example

We will:

1. Open the GUI.
2. Load the example SEGY and select crossline XL480.
3. Crop the profile, run the automated Lateral Interpretation pipeline,
   inspect automatic contours, run TSP/onlap/smoothing stages,
   and finally demonstrate a short manual edit step.
4. Save example `.npz` outputs and plots.

---

## Steps

### 1) Prepare and launch

1. From the repository root, check the environment and install if needed:

   ```bash
   python check_env.py
   # then, if needed
   pip install -r requirements.txt
   ```
2. Launch the GUI:

   ```bash
   python run_gui.py
   ```

---

### 2) Load the example SEGY

1. In the GUI, click **Browse** and navigate to:

   ```
   example_data/F3_Dip_steered_median_subvolume_IL230-430_XL475-675_T1600-1800-CDP170-370.sgy
   ```
2. Set seismic attributes in the main GUI panel before loading:

   * **Initial inline/crossline attribute**: `475`
     (the example subvolume crosslines start at XL475; we will study XL480)
   * **F3 Time Sample Interval**: keep default `4`
   * **x_start**: `170`
   * **y_start**: `1600`
   * Select the **Crossline** radio button (not Inline)
   * Set **Manual interpretation**: `No` (we will do manual steps later)
3. Click **Load Selected File**.
   In the file selection window that opens, choose `480` (the GUI will list available inline/crossline indices) — this loads XL480.

---

### 3) Crop the profile for the study area

On the plotted seismic profile, **click and drag** to select the study window (or press and drag to select entire profile).
Select the *whole* visible profile for this tutorial.

![Crop selection animation](tutorial_images/01.gif)
*Figure: click-and-drag to crop the study window.*

---

### 4) Run Lateral Interpretation — contour extraction

1. Press **Proceed to Lateral Interpretation**.
2. In the new LI window:

   * Change **Line Length Threshold** to `10` (50 is too large for this zoomed view).
   * Leave other parameters as default.
3. Click **Run contour extraction & save contours.npz**.

A small status window will pop up reporting the number of initial layers detected.

> **Note:** `contours.npz` (and other intermediate `.npz` outputs) are saved to the current working directory (the folder from which you launched `run_gui.py`).

---

### 5) Optional manual merge step (skip for now)

When the contour-extraction canvas opens, you can manually merge or edit layers. For this quick tutorial press **Finish Merging** to continue with automated steps.

---

### 6) TSP, Onlap/Downlap, Vertical smoothing

For each stage (TSP connection, onlap/downlap detection, vertical smoothing):

1. Use the default parameter values (unless noted otherwise).
2. Click the corresponding **Run ...** button to execute the stage.
3. The GUI saves outputs after each stage (for example `tsp.npz`, `onlap.npz`, `smooth.npz`).

---

### 7) Plot & save interpreted results

1. Go to **More → Plotting**.
2. For a simple illustrative plot:

   * **Source (npz)**: select `feature.npz`
   * **Overlay**: `Original`
3. Click **Plot and Save**.

This will produce a plot image (saved in the working directory). Example of the automatic interpretation for XL480:

![Automatic interpretation example](tutorial_images/XL480_interpreted.png)
*Figure: automatic LI result for XL480 (feature overlay on original).*

---

### 8) Manual interpretation (short example)

1. Go to **More → Manual Interpretation**.
2. Click on a red layer outline in the canvas to select it, then press **Merge** (or use provided edit tools).
![Manual edit example](tutorial_images/manual.gif)
*Figure: a short manual edit.*
3. When done, press **Finish**.
4. A **Contour Detection** dialog will appear. Use these settings:

   * **Use Manual Interpretation**: `Yes`
   * **Manual Interpretation Type**: `Layer`
   * **Lateral Interpretation Setting**: `None`
5. Click **Run**.

After the manual step, return to **More → Plotting**, choose `manual.npz`, set **Overlay** to `None`, then **Plot and Save**.

![Manual edit example](tutorial_images/XL480_manual.png)
*Figure: XL480 after a short manual edit (manual overlay on original).*

---

## Expected files produced (examples)

* `contours.npz` — extracted contours
* `tsp.npz` — after TSP connection
* `onlap.npz` — onlap/downlap result
* `feature.npz` — features used for plotting
* `manual.npz` — after manual edits
* plot images (`*.png`) saved by the plotting dialog

> These files are saved to the directory where `run_gui.py` was launched (unless the GUI save path is changed).

---

## Troubleshooting & tips

* If **no layers** are detected at contour extraction:

  * Try lowering **Line Length Threshold**
  * Ensure you selected the correct crossline/inline index
* If contours look noisy:

  * Increase smoothing parameter in the smoothing stage
* If GUI cannot read the SEGY file:

  * Make sure `cigsegy` is installed and correct SEGY path is supplied
  * Check `check_env.py` for missing dependencies

---

## Where to find more detail

* For full method description and validation, **read the paper**: *Lateral Interpretation: A New Method for High Resolution Seismic Stratigraphy*.
* For GUI-specific usage details and credits, open the **Help** menu inside the Lateral Interpretation GUI.

````

---