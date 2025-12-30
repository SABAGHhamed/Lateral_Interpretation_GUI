# Lateral Interpretation GUI

This repository contains the Python source code for the **Lateral Interpretation (LI)**
method described in:

> Sabagh, H., Naruse, H., *Lateral Interpretation: A New Method for High Resolution Seismic Stratigraphy*, Computers & Geosciences (under review).

## Requirements
- Python ≥ 3.8
- Windows / Linux / macOS
- No admin privileges required

### Primary dataset (used in the paper)

The main experiments in the accompanying manuscript are based on the
raw F3 seismic volume provided in the *F3 Demo 2023* release.

Source:
https://terranubis.com/datainfo/F3-Demo-2023

After downloading the dataset, the seismic file used in this study is
located at:

F3_Demo_2023/Rawdata/Seismic_data.sgy

Due to the large file size (≈680 MB), this seismic volume is **not included**
in the repository. Users and reviewers are encouraged to obtain the data
directly from the official source.

**Dutch F3 Offshore Dataset** (offshore Netherlands) is publicly available for research and
educational use under the **CC-BY-SA license**.

---

### Small example dataset (for quick testing)

For demonstration and quick testing purposes, this project is also
compatible with a small subvolume extracted from the Dutch F3 dataset:

F3_Dip_steered_median_subvolume_IL230-430_XL475-675_T1600-1800-CDP170-370.sgy

This file is provided via the **segyio-notebooks** repository maintained by
**Equinor**:

https://github.com/equinor/segyio-notebooks/tree/master/data/basic

The limited spatial extent of the sample data makes it overly simplified for demonstrating
the full capabilities of the Lateral Interpretation method. Please use F3 raw dataset for an extensive test.


## Third-Party Libraries

This software relies on external open-source libraries. In particular, SEG-Y
file reading is handled using the **CIGSEGY** library:

Li, J., *CIGSEGY: A tool for exchanging data between SEG-Y format and NumPy array
inside Python environment*.

Project URL:
https://github.com/JintaoLee-Roger/cigsegy

Documentation:
https://cigsegy.readthedocs.io/en/latest/


## Quick Start

1. Unzip the archive
2. Open terminal in this folder
3. Run:

   python check_env.py

4. If needed:

   pip install -r requirements.txt

5. Launch GUI:

   python run_gui.py

## Tutorial

- 👉 [Quick tutorial: Interpret XL480 (example)](TUTORIAL.md)

## Notes
- Further details are provided within the paper and GUI Help section.
- New .npz files and plotting results might overwrite the existing ones in the main directory during each interpretation stage.
- No system-wide modifications are performed.

## License

This software is released under the MIT License. See the `LICENSE` file for details.
