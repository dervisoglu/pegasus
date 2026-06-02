# PEGASUS: Parametric Echelle Graphical Assistant for Spectroscopic Unification and Spline-fitting

**PEGASUS** is a modern, high-performance, interactive visual tool written in Python and PyQt5 designed for reducing, normalizing, spline-fitting, and merging stellar echelle spectra. 

Inspired by the spectral reduction methodologies detailed in *"Spectroscopically resolving the Algol triple system"* (Kolbas et al.), PEGASUS consolidates legacy multi-window command workflows into a unified, highly interactive graphical dashboard.

---

## Table of Contents
1. [Acronym Expansion](#acronym-expansion)
2. [Prerequisites & Installation](#prerequisites--installation)
3. [Quick Start](#quick-start)
4. [Extensive Usage Guide](#extensive-usage-guide)
   - [File Loading and Sorting](#1-file-loading-and-sorting)
   - [Continuum Anchor Placement and Fitting](#2-continuum-anchor-placement-and-fitting)
   - [Advanced Continuum Manipulation (Mouse Wheel)](#3-advanced-continuum-manipulation-mouse-wheel)
   - [Continuum Snap-to-Median snapping](#4-continuum-snap-to-median-snapping)
   - [Blaze Continuity (Copying and Interpolation)](#5-blaze-continuity-copying-and-interpolation)
   - [Reference Spectrum and Doppler Shift Adjustments](#6-reference-spectrum-and-doppler-shift-adjustments)
   - [Saving and Loading Blaze Configurations](#7-saving-and-loading-blaze-configurations)
   - [Visual Trimming and Scientific Merging](#8-visual-trimming-and-scientific-merging)
5. [Keyboard and Mouse Reference](#keyboard-and-mouse-reference)
6. [Data Interoperability Details](#data-interoperability-details)
7. [License & Credits](#license--credits)

---

## Acronym Expansion

- **P**arametric
- **E**chelle
- **G**raphical
- **A**ssistant for
- **S**pectroscopic
- **U**nification and
- **S**pline-fitting

---

## Prerequisites & Installation

### Core Dependencies
PEGASUS is fully compatible with Python 3.8 to 3.12. Ensure you have the standard scientific Python libraries and PyQt5 installed. 

Install the required packages in your active environment:
```bash
pip install pyqt5 matplotlib numpy scipy
```

### Running PEGASUS
Navigate to your repository and execute:
```bash
python pegasus.py
```

---

## Quick Start

1. **Launch** `pegasus.py`.
2. Click **Load Spectra** and select your raw echelle orders.
3. Left-click on the top plot to place continuum anchor points.
4. Scale all points together using the **Mouse Wheel** to line up with the spectrum.
5. Navigate through the orders using **Next** and **Previous** buttons.
6. Click **Merge Spectra** to trim overlapping boundaries and combine them into a single 1D spectrum.

---

## Extensive Usage Guide

### 1. File Loading and Sorting
PEGASUS parses plain-text ASCII files representing echelle orders. 
- **Format Requirements**: The input files should contain two columns separated by spaces or tabs: Wavelength (in Å) and Intensity/Flux. Text headers or invalid line rows are automatically skipped.
- **Sequential Sorting**: When you click **Load Spectra** and select multiple files, PEGASUS automatically sorts them lexicographically by filename. This guarantees that spectral orders are arranged sequentially in your memory, facilitating smooth previous/next navigation and edge merging.

### 2. Continuum Anchor Placement and Fitting
A crucial step in normalization is fitting the *blaze function* (continuum profile). PEGASUS offers two robust algorithms for this task:
- **Cubic Spline (Default)**: Best for complex instrument profiles that exhibit steep or non-polynomial instrument slopes. Requires at least 4 anchor points. If fewer points are placed, it falls back to Quadratic (3 points) or Linear (2 points) interpolation.
- **Polynomial**: Fits a polynomial curve to the anchors. Ideal for smooth, slowly-varying profiles.
  - *Polynomial Degree*: You can select degrees from 1 to 9.
  - *Mathematical Stability*: PEGASUS utilizes a centered wavelength model (subtracting the order's median wavelength $\lambda_{\text{mid}}$) to prevent computational overflow during high-order polynomial fits (an issue commonly found in raw fitting tools).
  - *Interactive Degree Tuning*: When using the Polynomial method, adding or removing points automatically scales the default polynomial degree to $N-1$ (up to degree 9) for rapid drafting.

### 3. Advanced Continuum Manipulation (Mouse Wheel)
Placing dozens of points on every order can be tedious. PEGASUS implements a **Mouse Scroll Wheel Scaling** feature that adjusts all continuum points on the active order simultaneously:
- **Multiply/Divide Mode (Default)**: Scrolling up multiplies all anchor $Y$-coordinates by the step size (e.g., `1.01` or `1.001`); scrolling down divides them. This acts as a scale multiplier, raising or lowering the entire continuum curve while preserving its fractional curvature.
- **Add/Remove Mode**: Scrolling up/down adds/subtracts a fixed constant value (e.g., `0.01` or `0.1`) to all anchors. This shifts the curve vertically without rescaling.
- *Tip*: Combine this with small step sizes (`1.001` or `0.01`) for ultra-fine adjustments to match your normalized continuum to $1.0$ exactly!

### 4. Continuum Snap-to-Median Snapping
In noisy spectra or orders with closely-packed absorption lines, manual cursor clicks can easily land in line cores, causing the continuum curve to sag. 
- To prevent this, toggle the **"Snap to 1Å Local Median"** checkbox.
- When enabled, any point you left-click to add, or any point you drag, will automatically snap its $Y$-coordinate to the **median intensity** of all raw spectrum pixels located inside a $\pm 0.5\text{ Å}$ wavelength window centered at the mouse cursor.
- This allows you to rapidly place points across absorption lines and remain confident that they align with the local envelope.

### 5. Blaze Continuity (Copying and Interpolation)
When processing sequential echelle orders, the instrument blaze function changes slowly. PEGASUS provides two shortcuts to copy continuum shapes across orders:
- **Copy from Order**: Select a source order in the dropdown. PEGASUS copies the point configuration, shifts it horizontally by the difference in starting wavelength $\Delta\lambda = \lambda_{\text{start, dest}} - \lambda_{\text{start, source}}$, and fits it to the new order.
- **Interpolate between Orders**: In cases where a single echelle order is heavily contaminated by a broad feature (such as the $H\alpha$ absorption profile) making continuum placement impossible, you can interpolate. Choose the previous order ($c1$) and the next order ($c2$). Click **Interpolate** to average their fitted continuum shapes and project the interpolated curve onto the active order.

### 6. Reference Spectrum and Doppler Shift Adjustments
To verify the quality of your normalization:
- Click **Load Ref Spec** to load a synthetic stellar model or high-SNR atlas (Vega, Arcturus, etc.) in the bottom overlap panel. It will be rendered as a dark gray reference curve.
- Enter a velocity value in the **Doppler Shift (Å)** text field and click **Shift** to translate the reference spectrum horizontally in real-time. This allows you to align absorption lines and ensure your continuum normalization is perfect.

### 7. Saving and Loading Blaze Configurations
Your work is fully restorable.
- **Save Blazes**: Exports all your picked continuum anchor coordinates into a plain-text database. The format matches legacy Java logs perfectly:
  1. Total number of orders.
  2. A list of anchor point counts for each order.
  3. All space-separated $(X, Y)$ anchor coordinates in sequence.
- **Load Blazes**: Re-imports a saved blaze database, immediately re-populating all anchor points and re-generating their spline/polynomial fits.

### 8. Visual Trimming and Scientific Merging
Echelle spectra suffer from low signal-to-noise ratios (SNR) and severe instrument roll-off at the overlapping edges of each order. Directly merging raw orders creates jagged, overlapping spikes. PEGASUS solves this with an interactive visual trimmer:
- Click **Merge Spectra** to launch the resizable trimmer workspace dialog.
- **Interactive Boundary Dragging**: Each order is displayed individually. Hovering over the royal blue dashed vertical lines changes your mouse cursor to a double-sided horizontal arrow (`Qt.SizeHorCursor`). Left-click and drag these lines to narrow the active wavelength range of the order.
- **Dynamic Masking**: Discarded edge data is visually shaded with a translucent gray overlay, showing exactly which wavelength ranges will be ignored.
- **Scientific Re-Gridding and Combining**:
  1. Once all orders have been trimmed, click **Merge & Save 1D Spectrum**.
  2. PEGASUS automatically calculates the optimal unified pixel bin spacing $\Delta\lambda$ based on the average resolution of your spectral segments.
  3. It generates a uniform, target wavelength grid.
  4. For every wavelength bin, PEGASUS performs a **robust median combine** of all overlapping pixel intensities. Taking the median successfully rejects outliers like cosmic ray hits, telluric spikes, and remnant edge roll-offs.
  5. Any rare physical gaps are safely reconstructed using linear interpolation.
  6. The finalized merged 1D spectrum is exported as a clean two-column ASCII file.

---

## Keyboard and Mouse Reference

### Top Fitting Canvas
- **Left-Click (Empty space)**: Add a new continuum anchor point.
- **Left-Click & Drag (On point)**: Smoothly drag and reposition the anchor point.
- **Right or Middle-Click (On point)**: Delete the anchor point.
- **Scroll Wheel**: Vertically scale or shift all anchors in the active order.

### Visual Trimmer Canvas (Popup Dialog)
- **Hover near boundary line**: Cursor changes to `SizeHorCursor`.
- **Left-Click & Drag boundary line**: Interactively adjust the `Trim Min` and `Trim Max` wavelength cutoffs.

---

## Data Interoperability Details

- **Normalized Spectrum Save**: Click **Save Normed** to select an export folder. PEGASUS saves each fitted order under `[original_name]-norm.[ext]`. The output contains two tab-separated columns:
  $$\lambda \quad \left(\frac{I_{\text{raw}}}{I_{\text{fit}}}\right)$$
- **Output Compatibility**: The exported files are fully compatible with standard astronomical software such as IRAF, PyRAF, and spectroscopy analysis routines.

---

## License & Credits
PEGASUS is designed and maintained as an open-source tool for astronomical echelle spectroscopy. It is inspired by the spectral reduction algorithms of the Kolbas et al. Algol triple system paper.
