#!/usr/bin/env python3
# =========================================================================
# PEGASUS: Parametric Echelle Graphical Assistant for Spectroscopic Unification and Spline-fitting
#
# A modern, high-performance interactive visual tool for reducing,
# normalizing, spline-fitting, and merging echelle stellar spectra.
#
# Inspired by "Spectroscopically resolving the Algol triple system" (Kolbas et al.)
# =========================================================================

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QGridLayout, QPushButton, QComboBox, QLabel, QLineEdit, 
    QFileDialog, QSplitter, QFrame, QMessageBox, QSlider, QCheckBox,
    QSizePolicy, QStatusBar, QGroupBox, QDialog, QListWidget, QListWidgetItem, QAbstractItemView
)
from PyQt5.QtCore import Qt, QSize
from PyQt5.QtGui import QFont, QIcon

# ---------------------------------------------------------
# SpectrumOrder Class
# Manages data for a single echelle order
# ---------------------------------------------------------
class SpectrumOrder:
    def __init__(self, filepath, wavelength=None, intensity=None):
        self.filepath = filepath
        self.filename = os.path.basename(filepath)
        
        # Load spectrum data from file or use pre-loaded arrays
        if wavelength is not None and intensity is not None:
            self.wavelength = np.array(wavelength)
            self.intensity = np.array(intensity)
        else:
            wavelength_list = []
            intensity_list = []
            with open(filepath, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            wavelength_list.append(float(parts[0]))
                            intensity_list.append(float(parts[1]))
                        except ValueError:
                            continue
                            
            self.wavelength = np.array(wavelength_list)
            self.intensity = np.array(intensity_list)
        
        # Intermediate/Result values
        self.pts_x = []         # Picked continuum wavelengths
        self.pts_y = []         # Picked continuum intensities
        self.fit_y = None       # Evaluated continuum curve
        self.norm_y = None      # Normalized spectrum (intensity / fit_y)
        
        self.degree = 3         # Fit polynomial degree (if Polynomial)
        self.fit_method = "Cubic Spline" # Default fit method is Spline
        self.is_fitted = False
        
        # Center coordinate (pomak) for polynomial stability & active trim limits
        if len(self.wavelength) > 0:
            self.wavelength_mid = (self.wavelength[0] + self.wavelength[-1]) / 2.0
            self.trim_min = self.wavelength[0]
            self.trim_max = self.wavelength[-1]
        else:
            self.wavelength_mid = 0.0
            self.trim_min = 0.0
            self.trim_max = 0.0

    def sort_points(self):
        if len(self.pts_x) > 1:
            x = np.array(self.pts_x)
            y = np.array(self.pts_y)
            sort_idx = np.argsort(x)
            self.pts_x = list(x[sort_idx])
            self.pts_y = list(y[sort_idx])

    def fit(self):
        if len(self.pts_x) == 0:
            self.fit_y = None
            self.norm_y = None
            self.is_fitted = False
            return
            
        self.sort_points()
        x_pts = np.array(self.pts_x)
        y_pts = np.array(self.pts_y)
        
        if self.fit_method == "Polynomial":
            deg = min(self.degree, len(x_pts) - 1)
            if deg < 0:
                deg = 0
            if deg == 0:
                mean_y = np.mean(y_pts)
                self.fit_y = np.full_like(self.wavelength, mean_y)
            else:
                # Centered fitting to prevent mathematical overflow (exactly like Java pomak!)
                x_shifted = x_pts - self.wavelength_mid
                coeffs = np.polyfit(x_shifted, y_pts, deg)
                self.fit_y = np.polyval(coeffs, self.wavelength - self.wavelength_mid)
        else:  # Spline fit
            if len(x_pts) < 2:
                mean_y = np.mean(y_pts)
                self.fit_y = np.full_like(self.wavelength, mean_y)
            elif len(x_pts) == 2:
                # Linear interpolator
                f = interp1d(x_pts, y_pts, kind='linear', fill_value='extrapolate')
                self.fit_y = f(self.wavelength)
            elif len(x_pts) == 3:
                # Quadratic interpolator
                f = interp1d(x_pts, y_pts, kind='quadratic', fill_value='extrapolate')
                self.fit_y = f(self.wavelength)
            else:
                # Cubic Spline
                try:
                    cs = CubicSpline(x_pts, y_pts, bc_type='natural')
                    self.fit_y = cs(self.wavelength)
                except Exception:
                    # Fallback to linear in case of scipy error
                    f = interp1d(x_pts, y_pts, kind='linear', fill_value='extrapolate')
                    self.fit_y = f(self.wavelength)
                    
        # Avoid division by zero
        self.fit_y = np.where(self.fit_y == 0, 1e-10, self.fit_y)
        self.norm_y = self.intensity / self.fit_y
        self.is_fitted = True

    def find_closest_point(self, mouse_x, mouse_y, ax):
        if not self.pts_x:
            return None, None
            
        # Convert picked points to display pixel coordinates for scale-independent distance
        points_pixels = ax.transData.transform(np.column_stack((self.pts_x, self.pts_y)))
        mouse_pixel = np.array([mouse_x, mouse_y])
        
        distances = np.linalg.norm(points_pixels - mouse_pixel, axis=1)
        closest_idx = np.argmin(distances)
        return closest_idx, distances[closest_idx]

    def snap_to_median(self, x):
        # Local snap to median intensity in a 1.0 Å window around x
        window = 0.5
        mask = (self.wavelength >= x - window) & (self.wavelength <= x + window)
        if np.any(mask):
            return np.median(self.intensity[mask])
        return None

    def auto_find_continuum(self):
        total_len = len(self.wavelength)
        if total_len == 0:
            return
            
        self.pts_x = []
        self.pts_y = []
        
        # Median filter size (e.g., 21 pixels to smooth out noise and absorption lines)
        from scipy.ndimage import median_filter
        smoothed_intensity = median_filter(self.intensity, size=21)
        
        # Select up to 15 unique indices
        indices = set()
        
        # 1. First point in the first 10 pixels (indices 0 to min(9, total_len-1))
        limit_first = min(10, total_len)
        idx_first = np.argmax(smoothed_intensity[0:limit_first])
        indices.add(idx_first)
        
        # 2. Last point in the last 10 pixels (indices max(0, total_len-10) to total_len-1)
        start_last = max(0, total_len - 10)
        idx_last = start_last + np.argmax(smoothed_intensity[start_last:total_len])
        indices.add(idx_last)
        
        # 3. Middle 13 points distributed between idx_first and idx_last
        needed = 15 - len(indices)
        if needed > 0 and idx_last > idx_first + 1:
            bin_edges = np.linspace(idx_first + 1, idx_last, needed + 1, dtype=int)
            for i in range(needed):
                s_idx = bin_edges[i]
                e_idx = bin_edges[i+1]
                if s_idx < e_idx:
                    bin_max = s_idx + np.argmax(smoothed_intensity[s_idx:e_idx])
                    indices.add(bin_max)
                    
        # If we still don't have enough points (e.g. due to overlap or small size),
        # add more points from the remaining range to make it up to 15 (or total_len)
        if len(indices) < 15 and total_len > len(indices):
            for i in np.linspace(0, total_len - 1, 15, dtype=int):
                indices.add(i)
                if len(indices) >= min(15, total_len):
                    break
                    
        # Convert to a sorted list of unique indices
        sorted_indices = sorted(list(indices))
        
        # Map indices to wavelength and intensity values
        for idx in sorted_indices:
            self.pts_x.append(self.wavelength[idx])
            self.pts_y.append(smoothed_intensity[idx])
            
        self.sort_points()
        self.fit_method = "Cubic Spline"
        self.fit()

    def copy_blaze(self, source_order):
        self.pts_x = []
        self.pts_y = []
        # Calculate start wavelength difference to shift points horizontally
        dx = self.wavelength[0] - source_order.wavelength[0]
        for sx, sy in zip(source_order.pts_x, source_order.pts_y):
            self.pts_x.append(sx + dx)
            self.pts_y.append(sy)
        self.degree = source_order.degree
        self.fit_method = source_order.fit_method
        self.fit()

    def interpolate_blazes(self, order1, order2):
        if not order1.is_fitted or not order2.is_fitted:
            return
            
        # Interpolate by normalized pixel coordinate (0 to 1) instead of absolute wavelength,
        # because the blaze profile shape is physically tied to the detector pixel coordinate
        # and adjacent orders are shifted in wavelength.
        x1 = np.linspace(0.0, 1.0, len(order1.wavelength))
        x2 = np.linspace(0.0, 1.0, len(order2.wavelength))
        x_self = np.linspace(0.0, 1.0, len(self.wavelength))
        
        fit1 = np.interp(x_self, x1, order1.fit_y)
        fit2 = np.interp(x_self, x2, order2.fit_y)
        
        self.fit_y = (fit1 + fit2) / 2.0
        
        # Sample points from the interpolated curve to show picked points
        self.pts_x = []
        self.pts_y = []
        n_pts = 15
        indices = np.linspace(0, len(self.wavelength) - 1, n_pts, dtype=int)
        for idx in indices:
            self.pts_x.append(self.wavelength[idx])
            self.pts_y.append(self.fit_y[idx])
            
        self.degree = 3
        self.fit_method = "Cubic Spline"
        self.fit()

# ---------------------------------------------------------
# InteractiveFittingCanvas Class
# Main plotting area with click, drag, scroll to fit blazes
# ---------------------------------------------------------
class InteractiveFittingCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(7, 4.5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
        self.main_window = parent
        
        self.fig.patch.set_facecolor('#ffffff')
        self.ax.set_facecolor('#ffffff')
        self.ax.tick_params(colors='#000000', grid_color='#e0e0e0')
        for spine in self.ax.spines.values():
            spine.set_color('#cccccc')
            
        self.order = None
        self.active_drag_idx = None
        
        # Connect Matplotlib events
        self.mpl_connect('button_press_event', self.on_button_press)
        self.mpl_connect('button_release_event', self.on_button_release)
        self.mpl_connect('motion_notify_event', self.on_motion_notify)
        self.mpl_connect('scroll_event', self.on_scroll)
        
        self.setCursor(Qt.CrossCursor)
        
        # Custom callbacks to notify UI window
        self.on_change_callback = None
        self.status_callback = None

    def set_order(self, order):
        self.order = order
        self.active_drag_idx = None
        self.draw_plot(reset_zoom=True)

    def draw_plot(self, reset_zoom=False):
        # Save current view limits if not resetting zoom
        xlim = self.ax.get_xlim() if not reset_zoom and self.order else None
        ylim = self.ax.get_ylim() if not reset_zoom and self.order else None
        
        self.ax.clear()
        self.ax.grid(True, color='#e0e0e0', linestyle='--')
        
        if self.order is None or len(self.order.wavelength) == 0:
            self.ax.text(0.5, 0.5, "No Spectrum Loaded", color='#888888',
                         ha='center', va='center', transform=self.ax.transAxes, fontsize=14)
            self.draw()
            return
            
        # Draw raw spectrum (dark steel blue)
        self.ax.plot(self.order.wavelength, self.order.intensity, color='#34495e', label='Spectrum', alpha=0.85, linewidth=1.2)
        
        # Draw fitted continuum (crimson red)
        if self.order.is_fitted and self.order.fit_y is not None:
            self.ax.plot(self.order.wavelength, self.order.fit_y, color='#e74c3c', label='Continuum Fit', linewidth=2.0)
            
        # Draw selected continuum points (sky blue diamonds/crosses)
        if self.order.pts_x:
            self.ax.scatter(self.order.pts_x, self.order.pts_y, color='#3498db', edgecolor='#2c3e50',
                            marker='D', s=45, zorder=5, label='Continuum Points')
            
        self.ax.set_title(f"Active Order: {self.order.filename}", color='#000000', fontsize=12)
        self.ax.set_xlabel("Wavelength (Å)", color='#000000')
        self.ax.set_ylabel("Intensity", color='#000000')
        self.ax.legend(facecolor='#ffffff', edgecolor='#cccccc', labelcolor='#000000')
        
        if not reset_zoom and xlim is not None:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        else:
            # Add padding
            dx = (self.order.wavelength[-1] - self.order.wavelength[0]) * 0.02
            self.ax.set_xlim(self.order.wavelength[0] - dx, self.order.wavelength[-1] + dx)
            ymin, ymax = np.min(self.order.intensity), np.max(self.order.intensity)
            dy = (ymax - ymin) * 0.05 if ymax != ymin else 1.0
            self.ax.set_ylim(ymin - dy, ymax + dy)
            
        self.draw()

    def on_button_press(self, event):
        if event.inaxes != self.ax or self.order is None:
            return
            
        # Check toolbar mode (avoid adding points while zooming/panning)
        try:
            if self.main_window.toolbar_fitting.mode != '':
                return
        except AttributeError:
            pass
            
        # Left click: Add or select to drag
        if event.button == 1:
            idx, dist_px = self.order.find_closest_point(event.x, event.y, self.ax)
            
            # Drag if click is close (within 30 pixels for high-DPI compatibility)
            if idx is not None and dist_px < 30:
                self.active_drag_idx = idx
                self.setCursor(Qt.ClosedHandCursor)
            else:
                # Add a new point (only if we clicked inside the plot axes)
                if event.xdata is not None and event.ydata is not None:
                    y_val = event.ydata
                    if self.main_window.chk_median_snap.isChecked():
                        median_y = self.order.snap_to_median(event.xdata)
                        if median_y is not None:
                            y_val = median_y
                            
                    self.order.pts_x.append(event.xdata)
                    self.order.pts_y.append(y_val)
                    self.order.sort_points()
                    
                    # Auto-adjust polynomial degree exactly like Java
                    n_pts = len(self.order.pts_x)
                    if self.order.fit_method == "Polynomial" and 1 < n_pts < 11:
                        self.order.degree = n_pts - 1
                        
                    self.order.fit()
                    self.draw_plot()
                    if self.on_change_callback:
                        self.on_change_callback()
                    
        # Right or Middle click: Delete closest point
        elif event.button in [2, 3]:
            idx, dist_px = self.order.find_closest_point(event.x, event.y, self.ax)
            if idx is not None and dist_px < 30:
                self.order.pts_x.pop(idx)
                self.order.pts_y.pop(idx)
                
                # Auto-adjust degree downwards
                n_pts = len(self.order.pts_x)
                if self.order.fit_method == "Polynomial" and 1 < n_pts < 11:
                    self.order.degree = n_pts - 1
                    
                self.order.fit()
                self.draw_plot()
                if self.on_change_callback:
                    self.on_change_callback()

    def on_motion_notify(self, event):
        if self.order is None:
            return
            
        # Hover coordinate updates
        if event.inaxes == self.ax:
            if self.status_callback:
                self.status_callback(f"Wavelength: {event.xdata:.2f} Å  |  Intensity: {event.ydata:.2f}")
                
            # Handle point dragging
            if self.active_drag_idx is not None:
                if event.xdata is not None and event.ydata is not None:
                    y_val = event.ydata
                    if self.main_window.chk_median_snap.isChecked():
                        median_y = self.order.snap_to_median(event.xdata)
                        if median_y is not None:
                            y_val = median_y
                            
                    self.order.pts_x[self.active_drag_idx] = event.xdata
                    self.order.pts_y[self.active_drag_idx] = y_val
                    self.order.fit()
                    self.draw_plot()
                    if self.on_change_callback:
                        self.on_change_callback()
            else:
                # Hover cursor feedback (Open hand when near a point, crosshair otherwise)
                idx, dist_px = self.order.find_closest_point(event.x, event.y, self.ax)
                if idx is not None and dist_px < 30:
                    self.setCursor(Qt.OpenHandCursor)
                else:
                    self.setCursor(Qt.CrossCursor)

    def on_button_release(self, event):
        if self.active_drag_idx is not None:
            self.active_drag_idx = None
            self.order.sort_points()
            self.order.fit()
            self.draw_plot()
            if self.on_change_callback:
                self.on_change_callback()
        self.setCursor(Qt.CrossCursor)

    def on_scroll(self, event):
        if self.order is None or not self.order.pts_y:
            return
            
        # Get adjust values from parent
        try:
            method = self.main_window.method_selector.currentText()
            step = float(self.main_window.by_how_much.currentText())
        except Exception:
            method = "Multiply/Divide"
            step = 1.01
            
        # Scale/adjust ALL picked points coordinates, exactly like Java mouse wheel!
        if event.button == 'up':
            if method == "Multiply/Divide":
                self.order.pts_y = [y * step for y in self.order.pts_y]
            else:
                self.order.pts_y = [y + step for y in self.order.pts_y]
        elif event.button == 'down':
            if method == "Multiply/Divide":
                self.order.pts_y = [y / step for y in self.order.pts_y]
            else:
                self.order.pts_y = [y - step for y in self.order.pts_y]
                
        self.order.fit()
        self.draw_plot()
        if self.on_change_callback:
            self.on_change_callback()

# ---------------------------------------------------------
# NormalizedOverlapCanvas Class
# Bottom plotting canvas for visualizing order overlaps
# ---------------------------------------------------------
class NormalizedOverlapCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(7, 4.0), dpi=100)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
        
        self.fig.patch.set_facecolor('#ffffff')
        self.ax.set_facecolor('#ffffff')
        self.ax.tick_params(colors='#000000', grid_color='#e0e0e0')
        for spine in self.ax.spines.values():
            spine.set_color('#cccccc')

    def draw_overlap(self, orders, active_idx, draw_adjacent=True, ref_wavelength=None, ref_intensity=None, from_y=0.0, reset_axes=False):
        xlim = self.ax.get_xlim() if not reset_axes else None
        ylim = self.ax.get_ylim() if not reset_axes else None

        self.ax.clear()
        self.ax.grid(True, color='#e0e0e0', linestyle='--')
        
        if not orders or active_idx < 0 or active_idx >= len(orders):
            self.ax.text(0.5, 0.5, "No Normalized Spectrum to Plot", color='#888888',
                         ha='center', va='center', transform=self.ax.transAxes, fontsize=14)
            self.draw()
            return
            
        active_order = orders[active_idx]
        
        # Continuum line at 1.0 (crimson red)
        self.ax.axhline(1.0, color='#e74c3c', linestyle='--', linewidth=1.5, alpha=0.6, label='Continuum (1.0)')
        
        # Draw reference spectrum in dark gray if loaded
        if ref_wavelength is not None and ref_intensity is not None:
            self.ax.plot(ref_wavelength, ref_intensity, color='#7f8c8d', label='Reference', linewidth=1.0, alpha=0.7)
            
        # Draw previous order (cyan/blue)
        if draw_adjacent and active_idx > 0:
            prev_order = orders[active_idx - 1]
            if prev_order.is_fitted and prev_order.norm_y is not None:
                self.ax.plot(prev_order.wavelength, prev_order.norm_y, color='#3498db', 
                             label=f'Prev: {prev_order.filename}', alpha=0.6, linewidth=1.0)
                             
        # Draw next order (green)
        if draw_adjacent and active_idx < len(orders) - 1:
            next_order = orders[active_idx + 1]
            if next_order.is_fitted and next_order.norm_y is not None:
                self.ax.plot(next_order.wavelength, next_order.norm_y, color='#27ae60', 
                             label=f'Next: {next_order.filename}', alpha=0.6, linewidth=1.0)
                             
        # Draw current active normalized order (slate blue)
        if active_order.is_fitted and active_order.norm_y is not None:
            self.ax.plot(active_order.wavelength, active_order.norm_y, color='#2c3e50', 
                         label=f'Current: {active_order.filename}', linewidth=1.5)
        else:
            self.ax.text(0.5, 0.3, "(Current order not fitted yet)", color='#e74c3c',
                         ha='center', va='center', transform=self.ax.transAxes, fontsize=11)
            
        self.ax.set_title("Normalized Spectrum & Adjacent Order Overlaps", color='#000000', fontsize=12)
        self.ax.set_xlabel("Wavelength (Å)", color='#000000')
        self.ax.set_ylabel("Normalized Intensity", color='#000000')
        
        # Preserve zoom bounds unless resetting
        if not reset_axes and xlim is not None and ylim is not None and xlim[0] < xlim[1]:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        else:
            # Set limits automatically based on current, previous, and next orders if available
            w_min = active_order.wavelength[0]
            w_max = active_order.wavelength[-1]
            if draw_adjacent and active_idx > 0:
                w_min = min(w_min, orders[active_idx - 1].wavelength[0])
            if draw_adjacent and active_idx < len(orders) - 1:
                w_max = max(w_max, orders[active_idx + 1].wavelength[-1])
                
            self.ax.set_xlim(w_min, w_max)
            self.ax.set_ylim(from_y, 1.15)
            
        # Limit control
        self.ax.legend(facecolor='#ffffff', edgecolor='#cccccc', labelcolor='#000000', loc='lower left')
        self.draw()

# ---------------------------------------------------------
# PegasusWindow Class
# Unified dashboard interface
# ---------------------------------------------------------
class PegasusWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PEGASUS - Stellar Echelle Spectra Normalizer & Merger")
        self.resize(1400, 900)
        
        # State variables
        self.orders = []
        self.active_idx = -1
        
        # Reference spectrum original copy (for Doppler shifts)
        self.ref_wavelength_orig = None
        self.ref_wavelength = None
        self.ref_intensity = None
        
        # Build UI layout
        self.init_ui()

    def init_ui(self):
        # Main layout splitter
        main_splitter = QSplitter(Qt.Horizontal)
        self.setCentralWidget(main_splitter)
        
        # Left Panel (Controls Sidebar)
        sidebar = QWidget()
        sidebar.setMaximumWidth(380)
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setContentsMargins(10, 10, 10, 10)
        sidebar_layout.setSpacing(12)
        
        # Group 1: Files Controls
        files_group = QGroupBox("File Actions")
        files_layout = QGridLayout(files_group)
        self.btn_load = QPushButton("Load Spectra")
        self.btn_save = QPushButton("Save Normed")
        self.btn_load_blazes = QPushButton("Load Blazes")
        self.btn_save_blazes = QPushButton("Save Blazes")
        self.btn_merge = QPushButton("Merge Spectra")
        self.btn_load_ref = QPushButton("Load Ref Spec")
        self.btn_load_fits = QPushButton("Load IRAF FITS")
        
        self.btn_load.clicked.connect(self.load_spectra)
        self.btn_save.clicked.connect(self.save_norm_spectra)
        self.btn_load_blazes.clicked.connect(self.load_blazes)
        self.btn_save_blazes.clicked.connect(self.save_blazes)
        self.btn_merge.clicked.connect(self.open_spectra_merger)
        self.btn_load_ref.clicked.connect(self.load_ref_spectrum)
        self.btn_load_fits.clicked.connect(self.load_iraf_fits)
        
        files_layout.addWidget(self.btn_load, 0, 0)
        files_layout.addWidget(self.btn_save, 0, 1)
        files_layout.addWidget(self.btn_load_blazes, 1, 0)
        files_layout.addWidget(self.btn_save_blazes, 1, 1)
        files_layout.addWidget(self.btn_merge, 2, 0)
        files_layout.addWidget(self.btn_load_ref, 2, 1)
        files_layout.addWidget(self.btn_load_fits, 3, 0, 1, 2)
        sidebar_layout.addWidget(files_group)
        
        # Group 2: Order Navigation
        nav_group = QGroupBox("Order Navigation")
        nav_layout = QHBoxLayout(nav_group)
        self.btn_prev = QPushButton("◀ Previous")
        self.btn_next = QPushButton("Next ▶")
        self.lbl_selected_order = QLabel("No spectra loaded")
        self.lbl_selected_order.setAlignment(Qt.AlignCenter)
        self.lbl_selected_order.setStyleSheet("font-weight: bold; color: #2c3e50;")
        
        self.btn_prev.clicked.connect(self.select_previous)
        self.btn_next.clicked.connect(self.select_next)
        self.btn_prev.setEnabled(False)
        self.btn_next.setEnabled(False)
        
        nav_layout.addWidget(self.btn_prev)
        nav_layout.addWidget(self.btn_next)
        sidebar_layout.addWidget(nav_group)
        sidebar_layout.addWidget(self.lbl_selected_order)
        
        # Group 3: Fitting Settings
        fit_group = QGroupBox("Fitting Controls")
        fit_layout = QGridLayout(fit_group)
        
        fit_layout.addWidget(QLabel("Method:"), 0, 0)
        self.fit_method_selector = QComboBox()
        self.fit_method_selector.addItems(["Polynomial", "Cubic Spline"])
        self.fit_method_selector.currentIndexChanged.connect(self.fit_method_changed)
        fit_layout.addWidget(self.fit_method_selector, 0, 1)
        
        self.lbl_degree = QLabel("Poly Degree:")
        fit_layout.addWidget(self.lbl_degree, 1, 0)
        self.deg_selector = QComboBox()
        self.deg_selector.addItems([str(i) for i in range(1, 10)])
        self.deg_selector.currentIndexChanged.connect(self.degree_changed)
        fit_layout.addWidget(self.deg_selector, 1, 1)
        
        self.chk_median_snap = QCheckBox("Snap to 1Å Local Median")
        self.chk_median_snap.setChecked(False)
        fit_layout.addWidget(self.chk_median_snap, 2, 0, 1, 2)
        
        self.chk_auto_continuum = QCheckBox("Auto-Find Continuum (15 pts)")
        self.chk_auto_continuum.setChecked(False)
        self.chk_auto_continuum.stateChanged.connect(self.auto_continuum_toggled)
        fit_layout.addWidget(self.chk_auto_continuum, 3, 0, 1, 2)
        
        self.btn_fit = QPushButton("Re-Fit")
        self.btn_fit.clicked.connect(self.re_fit_current)
        fit_layout.addWidget(self.btn_fit, 4, 0, 1, 2)
        sidebar_layout.addWidget(fit_group)
        
        # Group 4: Point Manipulation (Mouse wheel settings)
        manip_group = QGroupBox("Point Scale (Wheel)")
        manip_layout = QGridLayout(manip_group)
        
        manip_layout.addWidget(QLabel("Wheel Action:"), 0, 0)
        self.method_selector = QComboBox()
        self.method_selector.addItems(["Multiply/Divide", "Add/Remove"])
        manip_layout.addWidget(self.method_selector, 0, 1)
        
        manip_layout.addWidget(QLabel("Step Size:"), 1, 0)
        self.by_how_much = QComboBox()
        self.by_how_much.addItems(["1.001", "1.01", "1.1", "0.01", "0.1", "1", "10", "100"])
        manip_layout.addWidget(self.by_how_much, 1, 1)
        sidebar_layout.addWidget(manip_group)
        
        # Group 5: Copy / Interpolate Continuum (Blaze)
        blaze_group = QGroupBox("Blaze Continuity")
        blaze_layout = QGridLayout(blaze_group)
        
        blaze_layout.addWidget(QLabel("Copy from Order:"), 0, 0)
        self.blaze_selector = QComboBox()
        self.blaze_selector.currentIndexChanged.connect(self.copy_blaze_selected)
        blaze_layout.addWidget(self.blaze_selector, 0, 1)
        
        blaze_layout.addWidget(QLabel("Interpolate between:"), 1, 0, 1, 2)
        self.c1_selector = QComboBox()
        self.c2_selector = QComboBox()
        self.btn_interpolate = QPushButton("Interpolate")
        self.btn_interpolate.clicked.connect(self.interpolate_blazes)
        
        blaze_layout.addWidget(self.c1_selector, 2, 0)
        blaze_layout.addWidget(self.c2_selector, 2, 1)
        blaze_layout.addWidget(self.btn_interpolate, 3, 0, 1, 2)
        sidebar_layout.addWidget(blaze_group)
        
        # Group 6: Visualizations & Doppler Shift
        vis_group = QGroupBox("Visualization & Doppler")
        vis_layout = QGridLayout(vis_group)
        
        self.chk_draw_adjacent = QCheckBox("Draw Adjacent Orders")
        self.chk_draw_adjacent.setChecked(True)
        self.chk_draw_adjacent.stateChanged.connect(self.update_overlap_plot)
        vis_layout.addWidget(self.chk_draw_adjacent, 0, 0, 1, 2)
        
        vis_layout.addWidget(QLabel("Y-Limit (Norm Plot):"), 1, 0)
        self.from_selector = QComboBox()
        self.from_selector.addItems([f"{i/10:.1f}" for i in range(10)]) # 0.0 to 0.9
        self.from_selector.setCurrentIndex(0)
        self.from_selector.currentIndexChanged.connect(self.update_overlap_plot)
        vis_layout.addWidget(self.from_selector, 1, 1)
        
        vis_layout.addWidget(QLabel("Doppler Shift (Å):"), 2, 0)
        self.doppler_field = QLineEdit("0.0")
        self.btn_doppler_move = QPushButton("Shift")
        self.btn_doppler_move.clicked.connect(self.doppler_shift)
        
        vis_layout.addWidget(self.doppler_field, 2, 1)
        vis_layout.addWidget(self.btn_doppler_move, 3, 0, 1, 2)
        sidebar_layout.addWidget(vis_group)
        
        sidebar_layout.addStretch()
        main_splitter.addWidget(sidebar)
        
        # Right Panel (Plots Area)
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)
        
        plot_splitter = QSplitter(Qt.Vertical)
        right_layout.addWidget(plot_splitter)
        
        # Top Plot Widget (Fitting Plot)
        top_plot_widget = QWidget()
        top_plot_layout = QVBoxLayout(top_plot_widget)
        top_plot_layout.setContentsMargins(5, 5, 5, 5)
        
        self.canvas_fitting = InteractiveFittingCanvas(self)
        self.canvas_fitting.on_change_callback = self.on_fit_changed
        self.canvas_fitting.status_callback = self.update_status_coordinates
        self.toolbar_fitting = NavigationToolbar(self.canvas_fitting, self)
        
        # Add zoom sliders directly in top plotting area to emulate ZoomWindow
        zoom_controls = QHBoxLayout()
        zoom_controls.addWidget(QLabel("Wavelength Zoom: "))
        self.sld_zoom_min = QSlider(Qt.Horizontal)
        self.sld_zoom_max = QSlider(Qt.Horizontal)
        self.sld_zoom_min.valueChanged.connect(self.on_zoom_slider_changed)
        self.sld_zoom_max.valueChanged.connect(self.on_zoom_slider_changed)
        zoom_controls.addWidget(self.sld_zoom_min)
        zoom_controls.addWidget(self.sld_zoom_max)
        
        top_plot_layout.addWidget(self.toolbar_fitting)
        top_plot_layout.addWidget(self.canvas_fitting)
        top_plot_layout.addLayout(zoom_controls)
        
        plot_splitter.addWidget(top_plot_widget)
        
        # Bottom Plot Widget (Overlap Plot)
        bottom_plot_widget = QWidget()
        bottom_plot_layout = QVBoxLayout(bottom_plot_widget)
        bottom_plot_layout.setContentsMargins(5, 5, 5, 5)
        
        self.canvas_overlap = NormalizedOverlapCanvas(self)
        self.toolbar_overlap = NavigationToolbar(self.canvas_overlap, self)
        
        # Add overlap range slider directly below overlap plot to emulate ResultZoomWindow
        overlap_controls = QHBoxLayout()
        overlap_controls.addWidget(QLabel("Overlap Window: "))
        self.sld_overlap_min = QSlider(Qt.Horizontal)
        self.sld_overlap_max = QSlider(Qt.Horizontal)
        self.sld_overlap_min.valueChanged.connect(self.on_overlap_slider_changed)
        self.sld_overlap_max.valueChanged.connect(self.on_overlap_slider_changed)
        overlap_controls.addWidget(self.sld_overlap_min)
        overlap_controls.addWidget(self.sld_overlap_max)
        
        bottom_plot_layout.addWidget(self.toolbar_overlap)
        bottom_plot_layout.addWidget(self.canvas_overlap)
        bottom_plot_layout.addLayout(overlap_controls)
        
        plot_splitter.addWidget(bottom_plot_widget)
        main_splitter.addWidget(right_panel)
        
        # Status Bar
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.statusBar.showMessage("Welcome to PEGASUS. Load spectrum files to start.")

    # ---------------------------------------------------------
    # Core Controller Logic
    # ---------------------------------------------------------
    def load_spectra(self):
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select Echelle Spectra Orders", "", "Spectra Files (*.dat *.order* *.*)"
        )
        if not files:
            return
            
        # If the user selected a single FITS file, route it to the FITS loader automatically
        if len(files) == 1 and files[0].lower().endswith(('.fits', '.fit')):
            self.load_iraf_fits_file(files[0])
            return
            
        # Clean loading
        self.orders = []
        # Sort files by name so echelle orders are in sequence
        files.sort()
        
        for filepath in files:
            if filepath.lower().endswith(('.fits', '.fit')):
                QMessageBox.warning(
                    self, "Invalid File Type", 
                    f"FITS files ({os.path.basename(filepath)}) cannot be loaded as ASCII text. Please use 'Load IRAF FITS' instead."
                )
                continue
            try:
                order = SpectrumOrder(filepath)
                if len(order.wavelength) > 0:
                    self.orders.append(order)
            except Exception as e:
                QMessageBox.warning(self, "Load Error", f"Failed to load {filepath}: {str(e)}")
                
        if not self.orders:
            self.statusBar.showMessage("No valid spectra files loaded.")
            return
            
        # Sort echelle orders by increasing wavelength
        self.orders.sort(key=lambda o: o.wavelength[0])
            
        self.active_idx = 0
        self.update_ui_state()
        
        # Setup Zoom and Overlap sliders
        self.setup_sliders()
        
        if self.chk_auto_continuum.isChecked():
            for order in self.orders:
                order.auto_find_continuum()
            self.update_ui_state()
            
        self.canvas_fitting.set_order(self.orders[self.active_idx])
        self.update_overlap_plot(reset_axes=True)
        
        self.statusBar.showMessage(f"Successfully loaded {len(self.orders)} spectral orders.")

    def load_iraf_fits(self):
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Select IRAF Echelle FITS File", "", "FITS Files (*.fits *.fit)"
        )
        if not filepath:
            return
        self.load_iraf_fits_file(filepath)

    def load_iraf_fits_file(self, filepath):
        try:
            import astropy.io.fits as pyfits
            import re
            
            with pyfits.open(filepath) as hdul:
                hdul.verify('silentfix')
                header = hdul[0].header
                data = hdul[0].data
                
                if data is None or len(data.shape) not in [2, 3]:
                    raise ValueError("FITS primary data must be a 2D or 3D array (orders x pixels or bands x orders x pixels).")
                    
                # Re-construct WAT2 string
                wat_str = ''.join([header[k].ljust(68) for k in sorted(header.keys()) if k.startswith('WAT2_')])
                if not wat_str:
                    raise ValueError("No WAT2_ WCS headers found. This does not appear to be an IRAF multispec FITS file.")
                    
                specs = re.findall(r'spec(\d+)\s*=\s*\"([^\"]+)\"', wat_str)
                if not specs:
                    raise ValueError("Failed to parse spec lines from WAT2 headers.")
                    
                # Sort specs by spec_id
                specs.sort(key=lambda x: int(x[0]))
                
                parsed_orders = []
                for spec_id_str, spec_val in specs:
                    parts = spec_val.split()
                    if len(parts) < 6:
                        continue
                    spec_id = int(parts[0])
                    order_id = int(parts[1])
                    disp_type = int(parts[2])
                    w_start = float(parts[3])
                    w_delta = float(parts[4])
                    n_pixels = int(parts[5])
                    
                    if disp_type != 0:
                        raise ValueError(f"Non-linear dispersion (type {disp_type}) in order {order_id} is not supported.")
                        
                    idx = spec_id - 1
                    
                    # Support both 2D and 3D data arrays
                    if len(data.shape) == 3:
                        if idx < 0 or idx >= data.shape[1]:
                            continue
                        intensity = data[0, idx, :n_pixels]
                    else:
                        if idx < 0 or idx >= data.shape[0]:
                            continue
                        intensity = data[idx, :n_pixels]
                        
                    wavelength = w_start + np.arange(len(intensity)) * w_delta
                    
                    basename = os.path.basename(filepath)
                    order_filename = f"{basename}_order_{order_id:03d}.dat"
                    
                    # Create the SpectrumOrder object using pre-loaded data
                    order = SpectrumOrder(filepath, wavelength=wavelength, intensity=intensity)
                    order.filename = order_filename
                    parsed_orders.append(order)
                    
            if not parsed_orders:
                raise ValueError("No valid orders parsed from FITS file.")
                
            # Sort echelle orders by increasing wavelength
            parsed_orders.sort(key=lambda o: o.wavelength[0])
            
            # Clean load
            self.orders = parsed_orders
            self.active_idx = 0
            self.update_ui_state()
            
            # Setup Zoom and Overlap sliders
            self.setup_sliders()
            
            if self.chk_auto_continuum.isChecked():
                for order in self.orders:
                    order.auto_find_continuum()
                self.update_ui_state()
                
            self.canvas_fitting.set_order(self.orders[self.active_idx])
            self.update_overlap_plot(reset_axes=True)
            
            self.statusBar.showMessage(f"Successfully loaded {len(self.orders)} spectral orders from FITS.")
            
        except Exception as e:
            QMessageBox.critical(self, "FITS Load Error", f"Failed to load FITS file: {str(e)}")

    def setup_sliders(self):
        if not self.orders or self.active_idx < 0:
            return
            
        order = self.orders[self.active_idx]
        w_min, w_max = order.wavelength[0], order.wavelength[-1]
        
        # Block signals to prevent infinite loops during configuration
        self.sld_zoom_min.blockSignals(True)
        self.sld_zoom_max.blockSignals(True)
        
        # Top Plot Zoom range: 0.1 Å resolution
        self.sld_zoom_min.setRange(int(w_min * 10), int(w_max * 10))
        self.sld_zoom_max.setRange(int(w_min * 10), int(w_max * 10))
        self.sld_zoom_min.setValue(int(w_min * 10))
        self.sld_zoom_max.setValue(int(w_max * 10))
        
        self.sld_zoom_min.blockSignals(False)
        self.sld_zoom_max.blockSignals(False)
        
        # Overlap plot zoom range (ResultZoomWindow logic)
        self.setup_overlap_sliders()

    def setup_overlap_sliders(self):
        if not self.orders or self.active_idx < 0:
            return
            
        # The overlap window covers range from previous order start to next order end
        start_idx = max(0, self.active_idx - 1)
        end_idx = min(len(self.orders) - 1, self.active_idx + 1)
        
        total_w_min = self.orders[start_idx].wavelength[0]
        total_w_max = self.orders[end_idx].wavelength[-1]
        
        self.sld_overlap_min.blockSignals(True)
        self.sld_overlap_max.blockSignals(True)
        
        self.sld_overlap_min.setRange(int(total_w_min * 10), int(total_w_max * 10))
        self.sld_overlap_max.setRange(int(total_w_min * 10), int(total_w_max * 10))
        
        # Default overlap view is current order range
        curr_order = self.orders[self.active_idx]
        self.sld_overlap_min.setValue(int(curr_order.wavelength[0] * 10))
        self.sld_overlap_max.setValue(int(curr_order.wavelength[-1] * 10))
        
        self.sld_overlap_min.blockSignals(False)
        self.sld_overlap_max.blockSignals(False)

    def update_ui_state(self):
        if not self.orders:
            return
            
        # Navigation buttons state
        self.btn_prev.setEnabled(self.active_idx > 0)
        self.btn_next.setEnabled(self.active_idx < len(self.orders) - 1)
        
        # Labels and active selections
        self.lbl_selected_order.setText(f"Order {self.active_idx + 1} / {len(self.orders)}")
        
        # Update dropdowns
        self.blaze_selector.blockSignals(True)
        self.c1_selector.blockSignals(True)
        self.c2_selector.blockSignals(True)
        
        self.blaze_selector.clear()
        self.c1_selector.clear()
        self.c2_selector.clear()
        
        for i in range(len(self.orders)):
            item_text = f"Order {i+1}"
            self.blaze_selector.addItem(item_text)
            self.c1_selector.addItem(item_text)
            self.c2_selector.addItem(item_text)
            
        self.blaze_selector.setCurrentIndex(self.active_idx)
        self.c1_selector.setCurrentIndex(max(0, self.active_idx - 1))
        self.c2_selector.setCurrentIndex(min(len(self.orders) - 1, self.active_idx + 1))
        
        self.blaze_selector.blockSignals(False)
        self.c1_selector.blockSignals(False)
        self.c2_selector.blockSignals(False)
        
        # Fit controls update
        order = self.orders[self.active_idx]
        self.fit_method_selector.blockSignals(True)
        self.deg_selector.blockSignals(True)
        
        self.fit_method_selector.setCurrentText(order.fit_method)
        self.deg_selector.setCurrentText(str(order.degree))
        self.lbl_degree.setEnabled(order.fit_method == "Polynomial")
        self.deg_selector.setEnabled(order.fit_method == "Polynomial")
        
        self.fit_method_selector.blockSignals(False)
        self.deg_selector.blockSignals(False)

    def select_previous(self):
        if self.active_idx > 0:
            self.active_idx -= 1
            self.update_ui_state()
            self.setup_sliders()
            if self.chk_auto_continuum.isChecked() and not self.orders[self.active_idx].pts_x:
                self.orders[self.active_idx].auto_find_continuum()
                self.update_ui_state()
            self.canvas_fitting.set_order(self.orders[self.active_idx])
            self.update_overlap_plot(reset_axes=True)

    def select_next(self):
        if self.active_idx < len(self.orders) - 1:
            self.active_idx += 1
            self.update_ui_state()
            self.setup_sliders()
            if self.chk_auto_continuum.isChecked() and not self.orders[self.active_idx].pts_x:
                self.orders[self.active_idx].auto_find_continuum()
                self.update_ui_state()
            self.canvas_fitting.set_order(self.orders[self.active_idx])
            self.update_overlap_plot(reset_axes=True)

    def fit_method_changed(self):
        if self.active_idx < 0:
            return
        method = self.fit_method_selector.currentText()
        self.orders[self.active_idx].fit_method = method
        self.lbl_degree.setEnabled(method == "Polynomial")
        self.deg_selector.setEnabled(method == "Polynomial")
        
        self.orders[self.active_idx].fit()
        self.canvas_fitting.draw_plot(reset_zoom=False)
        self.update_overlap_plot()

    def degree_changed(self):
        if self.active_idx < 0:
            return
        deg = int(self.deg_selector.currentText())
        self.orders[self.active_idx].degree = deg
        self.orders[self.active_idx].fit()
        self.canvas_fitting.draw_plot(reset_zoom=False)
        self.update_overlap_plot()

    def re_fit_current(self):
        if self.active_idx < 0:
            return
        self.orders[self.active_idx].fit()
        self.canvas_fitting.draw_plot(reset_zoom=False)
        self.update_overlap_plot()

    def auto_continuum_toggled(self, state):
        if state == Qt.Checked:
            if self.orders:
                for order in self.orders:
                    if not order.pts_x:
                        order.auto_find_continuum()
                self.update_ui_state()
                self.canvas_fitting.draw_plot(reset_zoom=False)
                self.update_overlap_plot()

    def on_fit_changed(self):
        # Called automatically when canvas points are added/dragged/removed/scrolled
        self.update_overlap_plot()
        # Ensure degree matches automatically in UI (if changed via automatic fitting)
        self.deg_selector.blockSignals(True)
        self.deg_selector.setCurrentText(str(self.orders[self.active_idx].degree))
        self.deg_selector.blockSignals(False)

    def update_status_coordinates(self, text):
        self.statusBar.showMessage(text)

    # ---------------------------------------------------------
    # Zoom Sliders Event Handlers
    # ---------------------------------------------------------
    def on_zoom_slider_changed(self):
        if self.active_idx < 0:
            return
        val_min = self.sld_zoom_min.value() / 10.0
        val_max = self.sld_zoom_max.value() / 10.0
        
        if val_min < val_max:
            self.canvas_fitting.ax.set_xlim(val_min, val_max)
            self.canvas_fitting.draw()

    def on_overlap_slider_changed(self):
        if self.active_idx < 0:
            return
        val_min = self.sld_overlap_min.value() / 10.0
        val_max = self.sld_overlap_max.value() / 10.0
        
        if val_min < val_max:
            self.canvas_overlap.ax.set_xlim(val_min, val_max)
            self.canvas_overlap.draw()

    # ---------------------------------------------------------
    # Blaze Continuity Logic (Copy / Interpolate)
    # ---------------------------------------------------------
    def copy_blaze_selected(self):
        source_idx = self.blaze_selector.currentIndex()
        if source_idx == self.active_idx or source_idx < 0 or self.active_idx < 0:
            return
            
        source_order = self.orders[source_idx]
        dest_order = self.orders[self.active_idx]
        
        if not source_order.pts_x:
            QMessageBox.information(self, "No Points", "Source order contains no continuum points.")
            return
            
        dest_order.copy_blaze(source_order)
        self.canvas_fitting.draw_plot(reset_zoom=False)
        self.update_ui_state()
        self.update_overlap_plot()
        self.statusBar.showMessage(f"Copied continuum curve from Order {source_idx + 1}")

    def interpolate_blazes(self):
        idx1 = self.c1_selector.currentIndex()
        idx2 = self.c2_selector.currentIndex()
        
        if idx1 < 0 or idx2 < 0 or self.active_idx < 0:
            return
            
        order1 = self.orders[idx1]
        order2 = self.orders[idx2]
        dest_order = self.orders[self.active_idx]
        
        if not order1.is_fitted or not order2.is_fitted:
            QMessageBox.warning(self, "Interpolate Error", "Both bounding orders must be fitted first.")
            return
            
        dest_order.interpolate_blazes(order1, order2)
        self.canvas_fitting.draw_plot(reset_zoom=False)
        self.update_ui_state()
        self.update_overlap_plot()
        self.statusBar.showMessage(f"Interpolated blaze curve between Order {idx1 + 1} and Order {idx2 + 1}")

    # ---------------------------------------------------------
    # Reference Spectrum and Doppler Logic
    # ---------------------------------------------------------
    def load_ref_spectrum(self):
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Load Reference Spectrum File", "", "Data Files (*.dat *.txt *.*)"
        )
        if not filepath:
            return
            
        wavelength = []
        intensity = []
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        wavelength.append(float(parts[0]))
                        intensity.append(float(parts[1]))
                    except ValueError:
                        continue
                        
        if not wavelength:
            QMessageBox.critical(self, "Load Error", "Failed to parse reference spectrum data.")
            return
            
        self.ref_wavelength_orig = np.array(wavelength)
        self.ref_wavelength = np.array(wavelength)
        self.ref_intensity = np.array(intensity)
        
        self.doppler_field.setText("0.0")
        self.update_overlap_plot()
        self.statusBar.showMessage("Successfully loaded reference spectrum.")

    def doppler_shift(self):
        if self.ref_wavelength_orig is None:
            QMessageBox.warning(self, "No Reference", "Please load a reference spectrum first.")
            return
            
        try:
            shift = float(self.doppler_field.text())
        except ValueError:
            QMessageBox.warning(self, "Invalid Shift", "Please enter a valid numeric Doppler shift.")
            return
            
        # Shift in wavelength (Å)
        self.ref_wavelength = self.ref_wavelength_orig + shift
        self.update_overlap_plot()
        self.statusBar.showMessage(f"Shifted reference spectrum by {shift:.2f} Å.")

    def update_overlap_plot(self, reset_axes=False):
        draw_adjacent = self.chk_draw_adjacent.isChecked()
        from_y = float(self.from_selector.currentText())
        
        self.canvas_overlap.draw_overlap(
            self.orders, self.active_idx,
            draw_adjacent=draw_adjacent,
            ref_wavelength=self.ref_wavelength,
            ref_intensity=self.ref_intensity,
            from_y=from_y,
            reset_axes=reset_axes
        )
        
        # Realize sliders limits based on current view bounds
        xlim = self.canvas_overlap.ax.get_xlim()
        self.sld_overlap_min.blockSignals(True)
        self.sld_overlap_max.blockSignals(True)
        self.sld_overlap_min.setValue(int(xlim[0] * 10))
        self.sld_overlap_max.setValue(int(xlim[1] * 10))
        self.sld_overlap_min.blockSignals(False)
        self.sld_overlap_max.blockSignals(False)

    # ---------------------------------------------------------
    # Save/Load Continuum Blaze Points File
    # Matches original Java format perfectly!
    # ---------------------------------------------------------
    def save_blazes(self):
        if not self.orders:
            return
            
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save Continuum Blaze File", "", "Blaze Files (*.dat *.txt)"
        )
        if not filepath:
            return
            
        try:
            with open(filepath, 'w') as f:
                # Output format:
                # 1. Total number of orders
                f.write(f"{len(self.orders)}\n")
                # 2. Number of picked points for each order
                for order in self.orders:
                    f.write(f"{len(order.pts_x)}\n")
                # 3. coordinates (X Y) of all picked points
                for order in self.orders:
                    for x, y in zip(order.pts_x, order.pts_y):
                        f.write(f"{x} {y}\n")
                        
            self.statusBar.showMessage(f"Continuum blazes successfully saved to {os.path.basename(filepath)}")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Failed to save blazes: {str(e)}")

    def load_blazes(self):
        if not self.orders:
            QMessageBox.warning(self, "No Spectra", "Please load spectra files first.")
            return
            
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Load Continuum Blaze File", "", "Blaze Files (*.dat *.txt)"
        )
        if not filepath:
            return
            
        try:
            with open(filepath, 'r') as f:
                lines = [line.strip() for line in f if line.strip()]
                
            if not lines:
                return
                
            num_orders = int(lines[0])
            pts_counts = []
            curr_line = 1
            
            for _ in range(num_orders):
                pts_counts.append(int(lines[curr_line]))
                curr_line += 1
                
            for i in range(num_orders):
                count = pts_counts[i]
                if i < len(self.orders):
                    order = self.orders[i]
                    order.pts_x = []
                    order.pts_y = []
                    for _ in range(count):
                        parts = lines[curr_line].split()
                        order.pts_x.append(float(parts[0]))
                        order.pts_y.append(float(parts[1]))
                        curr_line += 1
                    order.fit()
                    
            self.update_ui_state()
            self.canvas_fitting.draw_plot(reset_zoom=False)
            self.update_overlap_plot()
            self.statusBar.showMessage(f"Continuum blazes successfully loaded from {os.path.basename(filepath)}")
        except Exception as e:
            QMessageBox.critical(self, "Load Error", f"Failed to parse blazes: {str(e)}")

    # ---------------------------------------------------------
    # Save Normalized Spectra Files
    # Appends -norm to filename in the output directory
    # ---------------------------------------------------------
    def save_norm_spectra(self):
        if not self.orders:
            return
            
        # Get output directory from user
        save_dir = QFileDialog.getExistingDirectory(self, "Select Save Directory")
        if not save_dir:
            return
            
        saved_count = 0
        for order in self.orders:
            if not order.is_fitted or order.norm_y is None:
                continue
                
            # Create output filename: basename-norm.ext
            base, ext = os.path.splitext(order.filename)
            norm_filename = f"{base}-norm{ext}"
            save_path = os.path.join(save_dir, norm_filename)
            
            try:
                with open(save_path, 'w') as f:
                    for w, ny in zip(order.wavelength, order.norm_y):
                        f.write(f"{w}\t{ny}\n")
                saved_count += 1
            except Exception as e:
                QMessageBox.warning(self, "Save Error", f"Failed to save {norm_filename}: {str(e)}")
                
        QMessageBox.information(
            self, "Save Complete", 
            f"Successfully normalized and exported {saved_count} of {len(self.orders)} orders."
        )
        self.statusBar.showMessage(f"Successfully saved {saved_count} normalized orders.")

    def open_spectra_merger(self):
        if not self.orders:
            QMessageBox.warning(self, "No Spectra", "Please load echelle spectra orders first.")
            return
            
        # Ensure at least some orders are fitted
        fitted_orders = [o for o in self.orders if o.is_fitted]
        if not fitted_orders:
            QMessageBox.warning(self, "No Fits", "You must fit at least one echelle order before merging.")
            return
            
        self.merger_window = SpectraMergerWindow(self.orders, self)
        self.merger_window.show()

# ---------------------------------------------------------
# SpectraMergerWindow Class
# Popup visual workspace for visually trimming and merging echelle orders
# ---------------------------------------------------------
class SpectraMergerWindow(QDialog):
    def __init__(self, orders, parent=None):
        super().__init__(parent)
        self.orders = orders
        self.parent_window = parent
        self.setWindowTitle("PEGASUS - Interactive Spectra Merger & Trimmer")
        self.resize(1200, 750)
        
        # Explicitly allow resizing, maximizing, and minimizing on all platforms
        self.setWindowFlags(self.windowFlags() | Qt.WindowMinMaxButtonsHint)
        
        # State variables
        self.active_order_idx = 0
        self.active_trim_drag = None  # 'left' or 'right'
        
        self.init_ui()
        self.draw_merger_plot()

    def init_ui(self):
        # Using standard Qt light theme
        
        main_layout = QHBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(12)
        
        # Left Panel (Sidebar Controls)
        sidebar = QWidget()
        sidebar.setMaximumWidth(320)
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setContentsMargins(0, 0, 0, 0)
        sidebar_layout.setSpacing(10)
        
        # Order Selector List
        list_group = QGroupBox("Echelle Orders")
        list_layout = QVBoxLayout(list_group)
        self.order_list_widget = QListWidget()
        self.order_list_widget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.order_list_widget.currentRowChanged.connect(self.select_order)
        list_layout.addWidget(self.order_list_widget)
        sidebar_layout.addWidget(list_group)
        
        # Trim Controls Group
        trim_group = QGroupBox("Wavelength Trim Bounds")
        trim_layout = QGridLayout(trim_group)
        
        self.txt_trim_min = QLineEdit()
        self.txt_trim_max = QLineEdit()
        self.txt_trim_min.editingFinished.connect(self.on_trim_text_changed)
        self.txt_trim_max.editingFinished.connect(self.on_trim_text_changed)
        
        trim_layout.addWidget(QLabel("Trim Min (Å):"), 0, 0)
        trim_layout.addWidget(self.txt_trim_min, 0, 1)
        trim_layout.addWidget(QLabel("Trim Max (Å):"), 1, 0)
        trim_layout.addWidget(self.txt_trim_max, 1, 1)
        
        self.btn_reset_trim = QPushButton("Reset Current Trim")
        self.btn_reset_trim.clicked.connect(self.reset_current_trim)
        trim_layout.addWidget(self.btn_reset_trim, 2, 0, 1, 2)
        sidebar_layout.addWidget(trim_group)
        
        # Action Buttons
        self.btn_recenter = QPushButton("🔍 Re-center View")
        self.btn_recenter.clicked.connect(self.recenter_view)
        sidebar_layout.addWidget(self.btn_recenter)
        
        sidebar_layout.addStretch()
        
        self.btn_merge = QPushButton("🚀 Merge & Save 1D Spectrum")
        self.btn_merge.setStyleSheet("""
            QPushButton {
                background-color: #27ae60;
                font-weight: bold;
                font-size: 14px;
                color: white;
                padding: 10px;
                border-color: #2ecc71;
            }
            QPushButton:hover {
                background-color: #2ecc71;
            }
        """)
        self.btn_merge.clicked.connect(self.merge_and_save)
        sidebar_layout.addWidget(self.btn_merge)
        
        main_layout.addWidget(sidebar)
        
        # Right Panel (Interactive Trimming Plot)
        plot_panel = QWidget()
        plot_layout = QVBoxLayout(plot_panel)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        
        # Setup Figure and Canvas
        self.fig = Figure(figsize=(7, 5), dpi=100)
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(111)
        self.fig.patch.set_facecolor('#ffffff')
        self.ax.set_facecolor('#ffffff')
        self.ax.tick_params(colors='#000000', grid_color='#e0e0e0')
        for spine in self.ax.spines.values():
            spine.set_color('#cccccc')
            
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        
        main_layout.addWidget(plot_panel)
        
        # Load orders into QListWidget
        for i, o in enumerate(self.orders):
            status = "Fitted" if o.is_fitted else "Unfitted"
            item = QListWidgetItem(f"Order {i+1:02d} ({status}) - {o.filename}")
            if not o.is_fitted:
                item.setFlags(item.flags() & ~Qt.ItemIsEnabled) # disable unfitted orders
            self.order_list_widget.addItem(item)
            
        # Select first fitted order
        for i, o in enumerate(self.orders):
            if o.is_fitted:
                self.order_list_widget.setCurrentRow(i)
                break
                
        # Connect Matplotlib Canvas events
        self.canvas.mpl_connect('button_press_event', self.on_button_press)
        self.canvas.mpl_connect('button_release_event', self.on_button_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion_notify)
        self.canvas.setCursor(Qt.CrossCursor)

    def select_order(self, row):
        if row < 0 or row >= len(self.orders):
            return
        self.active_order_idx = row
        order = self.orders[self.active_order_idx]
        
        # Populate text bounds
        self.txt_trim_min.setText(f"{order.trim_min:.3f}")
        self.txt_trim_max.setText(f"{order.trim_max:.3f}")
        
        self.draw_merger_plot(reset_zoom=True)

    def on_trim_text_changed(self):
        order = self.orders[self.active_order_idx]
        try:
            val_min = float(self.txt_trim_min.text())
            val_max = float(self.txt_trim_max.text())
            
            # Bound validation
            if val_min < val_max and order.wavelength[0] <= val_min <= order.wavelength[-1] and order.wavelength[0] <= val_max <= order.wavelength[-1]:
                order.trim_min = val_min
                order.trim_max = val_max
                self.draw_merger_plot(reset_zoom=False)
        except ValueError:
            pass

    def reset_current_trim(self):
        order = self.orders[self.active_order_idx]
        order.trim_min = order.wavelength[0]
        order.trim_max = order.wavelength[-1]
        self.txt_trim_min.setText(f"{order.trim_min:.3f}")
        self.txt_trim_max.setText(f"{order.trim_max:.3f}")
        self.draw_merger_plot(reset_zoom=False)

    def recenter_view(self):
        self.draw_merger_plot(reset_zoom=True)

    def draw_merger_plot(self, reset_zoom=False):
        xlim = self.ax.get_xlim() if not reset_zoom else None
        ylim = self.ax.get_ylim() if not reset_zoom else None
        
        self.ax.clear()
        self.ax.grid(True, color='#e0e0e0', linestyle='--')
        
        active_order = self.orders[self.active_order_idx]
        
        # Plot all loaded normalized orders in faint colors
        for i, o in enumerate(self.orders):
            if o.is_fitted and o.norm_y is not None:
                if i == self.active_order_idx:
                    # Highlight selected order in crimson/royal red
                    self.ax.plot(o.wavelength, o.norm_y, color='#e74c3c', linewidth=2.0, zorder=5, label=f"Selected: Order {i+1}")
                else:
                    # Faint grey for non-active orders
                    self.ax.plot(o.wavelength, o.norm_y, color='#bdc3c7', linewidth=1.0, zorder=2, alpha=0.6)
                    
        # Render left and right draggable dashed boundary lines (nice royal blue)
        self.line_min = self.ax.axvline(active_order.trim_min, color='#3498db', linestyle='--', linewidth=2.0, zorder=6, label='Trim Bounds')
        self.line_max = self.ax.axvline(active_order.trim_max, color='#3498db', linestyle='--', linewidth=2.0, zorder=6)
        
        # Render dynamic gray masking overlay (shade discarded areas with soft gray translucent masks)
        self.ax.axvspan(active_order.wavelength[0], active_order.trim_min, color='#bdc3c7', alpha=0.3, zorder=3)
        self.ax.axvspan(active_order.trim_max, active_order.wavelength[-1], color='#bdc3c7', alpha=0.3, zorder=3)
        
        # Labels and title
        self.ax.set_title(f"Visual Trimmer: Order {self.active_order_idx+1:02d} ({active_order.filename})", color='#000000', fontsize=12)
        self.ax.set_xlabel("Wavelength (Å)", color='#000000')
        self.ax.set_ylabel("Normalized Intensity", color='#000000')
        self.ax.legend(facecolor='#ffffff', edgecolor='#cccccc', labelcolor='#000000', loc='lower left')
        
        # Redraw
        if not reset_zoom and xlim is not None and ylim is not None:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        else:
            # Re-center viewport tightly on the active order
            dx = (active_order.wavelength[-1] - active_order.wavelength[0]) * 0.05
            self.ax.set_xlim(active_order.wavelength[0] - dx, active_order.wavelength[-1] + dx)
            self.ax.set_ylim(0.0, 1.25)
            
        self.canvas.draw()

    # ---------------------------------------------------------
    # Mouse Interaction Event Handlers
    # ---------------------------------------------------------
    def on_button_press(self, event):
        if event.inaxes != self.ax:
            return
            
        # Left click grabs vertical lines
        if event.button == 1:
            active_order = self.orders[self.active_order_idx]
            
            # Check horizontal display pixel distance
            left_px = self.ax.transData.transform([active_order.trim_min, 0])[0]
            right_px = self.ax.transData.transform([active_order.trim_max, 0])[0]
            click_px = event.x
            
            if abs(click_px - left_px) < 25:
                self.active_trim_drag = 'left'
                self.canvas.setCursor(Qt.SizeHorCursor)
            elif abs(click_px - right_px) < 25:
                self.active_trim_drag = 'right'
                self.canvas.setCursor(Qt.SizeHorCursor)

    def on_motion_notify(self, event):
        active_order = self.orders[self.active_order_idx]
        
        # Update dragging coordinates
        if self.active_trim_drag is not None and event.xdata is not None:
            if self.active_trim_drag == 'left':
                active_order.trim_min = max(active_order.wavelength[0], min(event.xdata, active_order.trim_max - 0.5))
                self.txt_trim_min.setText(f"{active_order.trim_min:.3f}")
            elif self.active_trim_drag == 'right':
                active_order.trim_max = min(active_order.wavelength[-1], max(event.xdata, active_order.trim_min + 0.5))
                self.txt_trim_max.setText(f"{active_order.trim_max:.3f}")
                
            self.draw_merger_plot(reset_zoom=False)
        else:
            # Hover cursor visual feedback (resize cursor when near line)
            if event.inaxes == self.ax:
                left_px = self.ax.transData.transform([active_order.trim_min, 0])[0]
                right_px = self.ax.transData.transform([active_order.trim_max, 0])[0]
                mouse_px = event.x
                
                if abs(mouse_px - left_px) < 25 or abs(mouse_px - right_px) < 25:
                    self.canvas.setCursor(Qt.SizeHorCursor)
                else:
                    self.canvas.setCursor(Qt.CrossCursor)

    def on_button_release(self, event):
        if self.active_trim_drag is not None:
            self.active_trim_drag = None
            self.canvas.setCursor(Qt.CrossCursor)
            self.draw_merger_plot(reset_zoom=False)

    # ---------------------------------------------------------
    # Scientific Median-Binning Merging Algorithm
    # ---------------------------------------------------------
    def merge_and_save(self):
        # 1. Collect trimmed science-ready data
        segments = []
        for i, order in enumerate(self.orders):
            if not order.is_fitted or order.norm_y is None:
                continue
                
            # Filter wavelengths strictly between trim bounds
            mask = (order.wavelength >= order.trim_min) & (order.wavelength <= order.trim_max)
            if np.any(mask):
                segments.append((order.wavelength[mask], order.norm_y[mask], i + 1))
                
        if not segments:
            QMessageBox.warning(self, "No Trimmed Data", "No fitted or trimmed echelle orders are available for merging.")
            return
            
        # 2. Re-grid over continuous spectral envelope
        all_wavelengths = np.concatenate([seg[0] for seg in segments])
        w_min = np.min(all_wavelengths)
        w_max = np.max(all_wavelengths)
        
        # Calculate optimal unified pixel bin spacing (delta_lambda)
        spacings = [np.diff(seg[0]).mean() for seg in segments if len(seg[0]) > 1]
        delta_w = np.mean(spacings) if spacings else 0.05
        
        # Generate target uniform grid
        num_grid_points = int((w_max - w_min) / delta_w) + 1
        merged_w = np.linspace(w_min, w_max, num_grid_points)
        merged_y = np.zeros_like(merged_w)
        
        # 3. Robust median combining of overlaps
        # Combine points using median to suppress cosmic rays and blaze drops
        half_dw = delta_w / 2.0
        for i, w in enumerate(merged_w):
            intensities_in_bin = []
            for seg_w, seg_y, _ in segments:
                # Find elements within bin spacing
                bin_mask = (seg_w >= w - half_dw) & (seg_w <= w + half_dw)
                if np.any(bin_mask):
                    intensities_in_bin.extend(seg_y[bin_mask])
                    
            if intensities_in_bin:
                merged_y[i] = np.median(intensities_in_bin)
            else:
                # Fallback interpolation for rare physical gaps
                merged_y[i] = np.interp(w, all_wavelengths, np.concatenate([seg[1] for seg in segments]))
                
        # 4. Save merged spectrum
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save Merged 1D Spectrum", "", "Data Files (*.dat *.txt)"
        )
        if not filepath:
            return
            
        try:
            with open(filepath, 'w') as f:
                for w, y in zip(merged_w, merged_y):
                    f.write(f"{w:.5f}\t{y:.6f}\n")
            QMessageBox.information(
                self, "Merge Successful", 
                f"Successfully merged {len(segments)} orders into a continuous 1D spectrum!\nSaved to: {os.path.basename(filepath)}"
            )
            self.accept()
        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Failed to save merged 1D spectrum: {str(e)}")

# ---------------------------------------------------------
# Main Application Launch
# ---------------------------------------------------------
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PegasusWindow()
    window.show()
    sys.exit(app.exec_())
