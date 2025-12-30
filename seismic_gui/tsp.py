import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from typing import List, Tuple, Union
import hashlib
import math
import os
import time
import glob
from scipy.ndimage import map_coordinates
from scipy.ndimage import label, generate_binary_structure
from scipy.ndimage import binary_dilation
from scipy.interpolate import interp1d
import copy 
import cv2
from scipy.spatial import KDTree
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
from IPython.display import display, clear_output
from collections import defaultdict
from copy import deepcopy
from math import inf, sqrt
import itertools
import xml.etree.ElementTree as ET
from shapely.geometry import MultiLineString, LineString
from shapely.affinity import scale
from scipy.spatial.distance import cdist, directed_hausdorff
from sklearn.metrics import r2_score
import pandas as pd
import gc  # good for dumping half made debugging plots



def find_base(file_name):
    if "." in file_name:
        base, _ = file_name.rsplit(".", 1)
    else:
        base = file_name  # no extension        
    
    return base



def integer_fitting(lines):
    if isinstance(lines[0], tuple):  # Handle single line
        lines = [lines]
        
    int_lines = [[(int(x),int(y)) for (x, y) in line] for line in lines]
    return int_lines  if len(int_lines) > 1 else int_lines[0]

def filter_lines_by_length(lines, l_threshold = 5):
    if isinstance(lines[0], tuple):  # Handle single line
        lines = [lines]
    # Filter out lines with length less than the threshold
    filtered_lines = [line for line in lines if len(line) >= l_threshold]
    return filtered_lines  if len(filtered_lines) > 1 else filtered_lines[0]


def filter_lines_by_length_reverse(lines, l_threshold = 5):
    if isinstance(lines[0], tuple):  # Handle single line
        lines = [lines]
    # Filter out lines with length less than the threshold
    filtered_lines = [line for line in lines if len(line) < l_threshold]
    return filtered_lines  if len(filtered_lines) > 1 else filtered_lines[0]



def filter_by_y_interval(lines, y_low, y_high):

    if isinstance(lines[0], tuple):  # Handle single line
        lines = [lines]    
    '''for snipping seismic profile between two y limits'''
    interval_lines = [[(x, y) for (x, y) in line if y >= y_low and y <= y_high] for line in lines]
    return  interval_lines  if len(interval_lines) > 1 else interval_lines[0]
    
    

# def rescale_lines_to_image(lines, dimensions, minmax = None, offset=(0, 0)):
    # """
    # Rescale a list of lines' coordinates to fit within the dimensions of a 2D image.

    # Parameters:
    # lines (list of list of tuples): List of lines, where each line is a list of (x, y) coordinates.
    # dimensions (tuple): Dimensions of the image as (width, height).

    # Returns:
    # list of list of tuples: Rescaled lines where coordinates are adjusted to fit within the image dimensions.
    # """
    
    
    # if not lines or not dimensions:
        # raise ValueError("Lines and dimensions must be provided and non-empty.")
        
    # if isinstance(lines[0], tuple):  # Handle single line
        # lines = [lines]        

    # width, height = dimensions
        # # Determine the minimum and maximum values of x and y in the lines
    # all_x = [coord[0] for line in lines for coord in line]
    # all_y = [coord[1] for line in lines for coord in line]
    # min_x, max_x = min(all_x), max(all_x)
    # min_y, max_y = min(all_y), max(all_y)

    # if minmax:
        # min_x_test, max_x_test, min_y_test, max_y_test = minmax
        # min_x = min(min_x, min_x_test)
        # max_x = max(max_x, max_x_test)
        # min_y = min(min_y, min_y_test)
        # max_y = max(max_y, max_y_test)

    # # Handle edge case where all x or y values are identical
    # x_range = max(max_x - min_x, 1e-9)
    # y_range = max(max_y - min_y, 1e-9)

    # # Scale and normalize each coordinate
    # rescaled_lines = [
        # [
            # (
                # ((x - min_x) / x_range * (width - 1)),
                # ((y - min_y) / y_range * (height - 1))
            # )
            # for x, y in line      # original values were int(),int() but removing that enhances resolution as matplotlib makes it int at the end
        # ]
        # for line in lines
    # ]

    # return  rescaled_lines  if len(rescaled_lines) > 1 else rescaled_lines[0]


def rescale_lines_to_image(lines, dimensions, minmax=None, offset=(0, 0)):
    """
    Rescales lines and applies an offset to preserve absolute image positioning.
    """
    if not lines or not dimensions:
        raise ValueError("Lines and dimensions must be provided and non-empty.")
    
    is_single_line = isinstance(lines[0], tuple)
    input_lines = [lines] if is_single_line else lines

    width, height = dimensions
    off_x, off_y = offset # This handles your x_start and y_start
    
    # 1. Find the bounds of the data
    all_x = [coord[0] for line in input_lines for coord in line]
    all_y = [coord[1] for line in input_lines for coord in line]
    
    min_x, max_x = min(all_x), max(all_x)
    min_y, max_y = min(all_y), max(all_y)

    if minmax:
        mx_min, mx_max, my_min, my_max = minmax
        min_x, max_x = min(min_x, mx_min), max(max_x, mx_max)
        min_y, max_y = min(min_y, my_min), max(max_y, my_max)

    x_range = max(max_x - min_x, 1e-9)
    y_range = max(max_y - min_y, 1e-9)

    # 2. Scale and then ADD the offset
    rescaled_lines = [
        [
            (
                ((x - min_x) / x_range * (width - 1)) + off_x,
                ((y - min_y) / y_range * (height - 1)) + off_y
            )
            for x, y in line
        ]
        for line in input_lines
    ]

    return rescaled_lines if not is_single_line else rescaled_lines[0]





def remove_duplicate_x_points(stripe_coordinates):
    """
    Replaces points with duplicate x-values in each line of stripe_coordinates,
    by computing the integer mean of y-values for each group of duplicate x-values.

    Parameters:
        stripe_coordinates (list of list of tuple): List of lines,
            where each line is a list of (x, y) tuples.

    Returns:
        list of list of tuple: The corrected stripe_coordinates with one point per unique x,
        using the average y for each x.
    """
    if isinstance(stripe_coordinates[0], tuple):  # Handle single line
        stripe_coordinates = [stripe_coordinates]      
    
    corrected_lines = []
    for line in stripe_coordinates:
        x_to_ys = {}  # Dictionary to track all y's for each x
        for x, y in line:
            if x not in x_to_ys:
                x_to_ys[x] = []
            x_to_ys[x].append(y)

        # Calculate the average y for each x (rounded to int)
        corrected_line = [(x, sum(ys) // len(ys)) for x, ys in x_to_ys.items()]
        corrected_line.sort(key=lambda point: point[0])  # sort by x value
        corrected_lines.append(corrected_line)

    return  corrected_lines  if len(corrected_lines) > 1 else corrected_lines[0]


def sort_by_x(stripe_coordinates):

    if isinstance(stripe_coordinates[0], tuple):  # Handle single line
        stripe_coordinates = [stripe_coordinates]      
    
    corrected_lines = []
    for line in stripe_coordinates:
        line.sort(key=lambda point: point[0])  # sort by x value
        corrected_lines.append(line)

    return  corrected_lines  if len(corrected_lines) > 1 else corrected_lines[0]
    
    

def correct_lines_trend(stripe_coordinates):
    """
    Adjusts the order of lines in stripe_coordinates to ensure x-axis values are 
    in an increasing order for each line.

    Parameters:
        stripe_coordinates (list of list of tuple): List of lines, 
            where each line is a list of (x, y) tuples.

    Returns:
        list of list of tuple: The corrected stripe_coordinates with x-values in increasing order.
    """
    corrected_lines = []
    for line in stripe_coordinates:
        # Check if the x-values are in a decreasing trend
        if len(line) > 1 and all(line[i][0] > line[i+1][0] for i in range(len(line) - 1)):
            corrected_lines.append(line[::-1])  # Reverse the line
        else:
            corrected_lines.append(line)  # Keep the line as is
    return corrected_lines


def split_lines_by_x_trend(stripe_coordinates):
    """
    Splits lines based on changes in the x-axis incrementing trend, 
    ensuring that split segments overlap by one point.

    Parameters:
        stripe_coordinates (list of list of tuple): List of lines, where each line is a list of tuple coordinates (x, y).

    Returns:
        list of list of tuple: Modified list of lines after splitting.
    """
    modified_lines = []

    for line in stripe_coordinates:
        if len(line) < 2:
            # Skip lines with less than 2 points
            modified_lines.append(line)
            continue

        # Initialize the current sub-line and determine the initial trend
        sub_line = [line[0]]
        increasing = None  # Keeps track of the trend (None for initialization)

        for i in range(1, len(line)):
            prev_x = line[i - 1][0]
            curr_x = line[i][0]

            # Determine the current trend
            if curr_x > prev_x:
                current_trend = True  # Increasing
            elif curr_x < prev_x:
                current_trend = False  # Decreasing
            else:
                current_trend = increasing  # No change, continue the same trend

            # Check if the trend has reversed
            if increasing is not None and current_trend != increasing:
                # Add the current sub-line to modified_lines
                modified_lines.append(sub_line)
                # Start a new sub-line, ensuring it overlaps with the last point of the previous segment
                sub_line = [sub_line[-1]]  # Start with the last point of the current sub-line

            # Update the trend and add the current point to the sub-line
            increasing = current_trend
            sub_line.append(line[i])

        # Append the last sub-line
        if sub_line:
            modified_lines.append(sub_line)

    return modified_lines


# Find the maximum y-value to flip vertically
def flip_layers(lines):
    if isinstance(lines[0], tuple):
        lines = [lines]

    max_y = max(y for line in lines for _, y in line)

    # Now flip y-values vertically
    lines_ = [
        [(x, max_y - y) for x, y in line]
        for line in lines
    ]
    
    return  lines_ if len(lines_) > 1 else lines_[0]






def cv2_layers(digofile):
    # Assuming digofile_positive is your ndarray of 0s and 1s
    # Step 1: Ensure the data is binary
    binary_image = (digofile * 255).astype(np.uint8)  # Scale to 0-255

    # Step 2: Find contours
    # Use cv2.RETR_EXTERNAL to get only the outer boundaries of each stripe
    contours, _ = cv2.findContours(binary_image, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE) # seems more detailed
#     contours, _ = cv2.findContours(binary_image, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
#     contours, _ = cv2.findContours(binary_image, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Step 3: Extract stripe coordinates
    stripe_coordinates_ = []
    for contour in contours:
        # Each contour is a numpy array of shape (n, 1, 2)
        coords = contour.squeeze(axis=1).tolist()  # Remove redundant dimensions and convert to list
        stripe_coordinates_.append(coords)

    stripe_coordinates = flip_layers(stripe_coordinates_)  # flipping vertically for better plotting
    stripe_coordinates = sort_lines_by_x(stripe_coordinates)  # new; tricky little think I lost!
    stripe_coordinates = split_lines_by_x_trend(stripe_coordinates)
    stripe_coordinates = correct_lines_trend(stripe_coordinates)  # reverse any reversely incrementing layer
    stripe_coordinates = remove_duplicate_x_points(stripe_coordinates)
    stripe_coordinates = bresenham_line(stripe_coordinates)
    # below _reversed was excluded because occasionally manual layers are long
    min_value = min(digofile.shape[0], digofile.shape[1])                   # this filtering is wrong, because whatever is identified in manual interpretation should be preserved
    stripe_coordinates = filter_lines_by_length_reverse(stripe_coordinates, l_threshold=min_value)  # filters unnecessary spikes due to wrong contour detection of masks # because it seems too good to be true
#     Now, stripe_coordinates contains the list of coordinates for each stripe

    return stripe_coordinates
    


def euclidean_distance(point1, point2):
    x1, y1 = point1
    x2, y2 = point2
    return ((x2 - x1)**2 + (y2 - y1)**2) ** 0.5

    

def get_positive_negative_masks(arr):
    """Return binary masks for positive and negative values in a 2D array.

    Args:
        arr (np.ndarray): A 2D NumPy array.

    Returns:
        mask_p (np.ndarray): Binary mask where arr > 0.
        mask_n (np.ndarray): Binary mask where arr < 0.
    """
    # Ensure it's a 2D array
    if arr.ndim != 2:
        raise ValueError("Input array must be 2D.")
    
    mask_p = (arr > 0).astype(np.uint8)
    mask_n = (arr < 0).astype(np.uint8)

    return mask_p, mask_n


def separate_positive_negative_coords(arr):
    """Return coordinates of positive and negative values in a 2D array.

    Args:
        arr (np.ndarray): A 2D NumPy array.

    Returns:
        coordinates_p (list of tuple): Coordinates where arr > 0.
        coordinates_n (list of tuple): Coordinates where arr < 0.
    """
    # Ensure it's a 2D array
    if arr.ndim != 2:
        raise ValueError("Input array must be 2D.")
    
    # Get coordinates
    coordinates_p = list(zip(*np.where(arr > 0)))
    coordinates_n = list(zip(*np.where(arr < 0)))

    return coordinates_p, coordinates_n


def interval_scale(old_y, new_max = 951, old_max = 511):
    """Return seismic dimension y intervals based on default 512 pixels plot 
    old_y = snipping y based on old (512) figure plot 
    new_max = seismic global maximum when it starts from 0
    """

    new_y = old_y * (new_max / old_max)
    return int(new_y)
    
    

def magnet(lines, cluster_len_threshold=10, distance_threshold=1):
    """
    Merge lines whose points are within 'distance_threshold' pixels of each other.
    Stop merging if merged line exceeds 'cluster_len_threshold' length.
    """
    lines = [line.copy() for line in lines]

    offsets = [
        (dx, dy)
        for dx in range(-distance_threshold, distance_threshold + 1)
        for dy in range(-distance_threshold, distance_threshold + 1)
    ]

    def build_point_map(lines):
        point_map = {}
        for idx, line in enumerate(lines):
            if not line or len(line) >= cluster_len_threshold:
                continue
            for point in line:
                point_map[point] = idx
        return point_map

    changed = True
    while changed:
        changed = False
        point_map = build_point_map(lines)

        for i in range(len(lines)):
            if not lines[i] or len(lines[i]) >= cluster_len_threshold:
                continue

            merged_indices = set()

            for point in lines[i]:
                for dx, dy in offsets:
                    neighbor = (point[0] + dx, point[1] + dy)
                    j = point_map.get(neighbor)

                    if (j is not None and j != i and
                        j not in merged_indices and
                        lines[j] and len(lines[j]) < cluster_len_threshold and
                        len(lines[i]) < cluster_len_threshold):

                        lines[i].extend(lines[j])
                        lines[j] = []
                        merged_indices.add(j)
                        changed = True

                        # Enforce cluster_len_threshold immediately
                        if len(lines[i]) >= cluster_len_threshold:
                            break
                if len(lines[i]) >= cluster_len_threshold:
                    break

    return [line for line in lines if line]

    
    
    
def generate_cutting_grid(lines, square_size=10):
    """
    Generate a grid of bounding boxes (squares) that partition the bounding
    box of the dataset based on a given square size.

    Parameters:
        lines (list of list of tuple): List of polylines, each a list of (x, y) tuples.
        square_size (int): Width and height of each square.

    Returns:
        grid (np.ndarray): 2D array where each element is a tuple:
                           ((x_min, y_min), (x_max, y_max))
                           representing the frame of that grid cell.
    """
    # Flatten all points
    all_points = [pt for line in lines for pt in line if pt is not None and len(pt) == 2]
    if not all_points:
        raise ValueError("No valid points provided.")

    # Convert to numpy array
    points_array = np.array(all_points)

    # Determine global bounding box
    x_min, y_min = np.min(points_array, axis=0)
    x_max, y_max = np.max(points_array, axis=0)

    # Align to the square grid
    x_min_aligned = (x_min // square_size) * square_size
    y_min_aligned = (y_min // square_size) * square_size
    x_max_aligned = ((x_max + square_size - 1) // square_size) * square_size
    y_max_aligned = ((y_max + square_size - 1) // square_size) * square_size

    # Determine number of rows and columns
    n_cols = int((x_max_aligned - x_min_aligned) // square_size)
    n_rows = int((y_max_aligned - y_min_aligned) // square_size)

    # Initialize grid array
    grid = np.empty((n_rows, n_cols), dtype=object)

    for row in range(n_rows):
        for col in range(n_cols):
            cell_x_min = x_min_aligned + col * square_size #- 1e-3
            cell_y_min = y_min_aligned + row * square_size #- 1e-3
            cell_x_max = cell_x_min + square_size #+ 1e-3
            cell_y_max = cell_y_min + square_size #+ 1e-3
            grid[row, col] = ((cell_x_min, cell_y_min), (cell_x_max, cell_y_max))

    return grid



def allocate_segments_to_grid(lines, square_size=10):
    """
    Allocates continuous sub-lines from input lines into square grid cells
    and detects integrity disruptions (junctions) between cells.

    Parameters:
        lines (list of list of tuple): Input lines, each a list of (x, y) points.
        square_size (int): The size of each square grid cell.

    Returns:
        squares (np.ndarray): 2D grid where each cell contains a list of sub-lines.
        junctions (list of tuple): List of (point1, point2) tuples indicating 
                                   line continuation cut across grid cells.
    """
    lines = [line.copy() for line in lines]
    grid = generate_cutting_grid(lines, square_size)
    squares = np.empty_like(grid, dtype=object)
    for idx in np.ndindex(grid.shape):
        squares[idx] = []

    junctions = []

    for line in lines:
        for row in range(grid.shape[0]):
            for col in range(grid.shape[1]):
                (x_min, y_min), (x_max, y_max) = grid[row, col]

                current_subline = []
                last_point = None
                last_inside = False

                for pt in line:
                    x, y = pt
                    inside = x_min <= x <= x_max and y_min <= y <= y_max

                    if inside:
                        current_subline.append(pt)
                    else:
                        if last_inside and last_point is not None:
                            # Integrity cut from last_point (inside) to pt (outside)
                            junctions.append((last_point, pt))
                        if len(current_subline) >= 1:
                            squares[row, col].append(current_subline)
                        current_subline = []

                    last_inside = inside
                    last_point = pt

                # End of line case: if still inside
                if len(current_subline) > 1:
                    squares[row, col].append(current_subline)

    return squares#, junctions






def flatten_squares(squares):

    all_lines = []

    for row in range(squares.shape[0]):
        for col in range(squares.shape[1]):
            cell = squares[row, col]
            if cell:  # not empty
                all_lines.extend(cell)

    return all_lines


def find_junctions(square_lines, lines):
    """
    For each line in square_lines, determine if its start or end point is isolated.
    If so, find the adjacent point in the corresponding original line from `lines`
    and return the junctions as (point1, point2) pairs.

    Parameters:
        square_lines (list of list of tuple): Cut or segmented lines.
        lines (list of list of tuple): Original full lines before cutting.

    Returns:
        junctions (list of list of tuple): Each element is a list [point1, point2],
                                           where point1 is isolated and point2 is
                                           its adjacent point in the original line.
    """
    def point_to_line_map(lines):
        """Maps each point to its corresponding line index in lines."""
        point_to_line_idx = {}
        for idx, line in enumerate(lines):
            if not line:
                continue
            for pt in line:
                point_to_line_idx[pt] = idx
        return point_to_line_idx

    point_to_line_idx_lines = point_to_line_map(lines)

    # Collect all start and end points from original lines
    original_starts = set(line[0] for line in lines if line)
    original_ends = set(line[-1] for line in lines if line)

    junctions = []

    for sq_line in square_lines:
        if not sq_line:
            continue

        start = sq_line[0]
        end = sq_line[-1]

        if start not in original_starts:
            # It's an isolated start point
            point1 = start
            idx1 = point_to_line_idx_lines.get(point1)
            if idx1 is not None:
                for i, pt in enumerate(lines[idx1]):
                    if pt == point1 and i > 0:
                        point2 = lines[idx1][i - 1]
                        junctions.append([point2, point1])
                        break

        if end not in original_ends:
            # It's an isolated end point
            point1 = end
            idx1 = point_to_line_idx_lines.get(point1)
            if idx1 is not None:
                for i, pt in enumerate(lines[idx1]):
                    if pt == point1 and i < len(lines[idx1]) - 1:
                        point2 = lines[idx1][i + 1]
                        junctions.append([point1, point2])
                        break
    
    #remove duplicate junctions due to the fact that one inside end to start in square i can be start to end in square j
    unique_pairs = set()
    junctions_unique = []

    for pair in junctions:
        # Convert to a frozenset of tuples to make it order-independent and hashable
        pair_set = frozenset(pair)
        if pair_set not in unique_pairs:
            unique_pairs.add(pair_set)
            junctions_unique.append(pair)  # Keep the original structure: list of two tuples

    return junctions_unique


def merge_lines_by_junctions(square_lines, junctions):
    """
    Merge line segments in square_lines based on junction information.
    Continues merging until no further reduction in number of lines occurs.

    Parameters:
        square_lines (list of list of tuple): Fragmented or cut lines.
        junctions (list of list of tuple): Junctions as returned from `find_junctions`,
                                           each a [point1, point2] list.

    Returns:
        merged_lines (list of list of tuple): Updated square_lines with merged segments.
    """
    def point_to_line_map(lines):
        """Maps each point to its corresponding line index."""
        point_to_line_idx = {}
        for idx, line in enumerate(lines):
            if not line:
                continue
            for pt in line:
                point_to_line_idx[pt] = idx
        return point_to_line_idx

    merged_lines = square_lines.copy()
    prev_len = -1

    while prev_len != len(merged_lines):
        prev_len = len(merged_lines)
        point_to_line_idx = point_to_line_map(merged_lines)
        to_remove = set()

        for point1, point2 in junctions:
            idx1 = point_to_line_idx.get(point1)
            idx2 = point_to_line_idx.get(point2)

            if idx1 is None or idx2 is None or idx1 == idx2:
                continue

            # Merge idx2 into idx1                
            merged_lines[idx1] = merged_lines[idx1] + merged_lines[idx2]
            merged_lines[idx2] = []  # Mark idx2 as empty
#             to_remove.add(idx2)

            # Update the mapping for merged points
            for pt in merged_lines[idx1]:
                point_to_line_idx[pt] = idx1
            
            # Remove all points that were mapped to idx2
            point_to_line_idx = {pt: idx for pt, idx in point_to_line_idx.items() if idx != idx2}


        # Remove empty lines after each full junction pass
        merged_lines = [line for i,line in enumerate(merged_lines) if line]

    return merged_lines




def merge_lines_by_reference(reference_lines, lines, threshold=0.9):
    """
    Maps each line in `lines` to a reference line in `reference_lines` if it's a subset
    (based on point overlap with threshold), and merges all such lines by x-order.

    Parameters:
        reference_lines (list of list of (x, y)): Base lines to relate to.
        lines (list of list of (x, y)): Input lines to map and merge.
        threshold (float): Minimum overlap ratio to consider a line a subset.

    Returns:
        merged_and_unmatched_lines (list of list of (x, y)): Merged lines from matches + untouched unmatched lines.
    """

    # Convert each reference line to a set of points for fast comparison
    ref_sets = [set(map(tuple, r_line)) for r_line in reference_lines]

    # Map from ref_line index to list of matched lines
    ref_to_lines = defaultdict(list)

    # Track which lines are matched
    matched_flags = [False] * len(lines)

    for i, line in enumerate(lines):
        line_set = set(map(tuple, line))
        best_match = -1
        best_score = 0

        for idx, ref_set in enumerate(ref_sets):
            intersection = line_set & ref_set
            if not line:
                continue
            overlap_ratio = len(intersection) / len(line)
            if overlap_ratio >= threshold and overlap_ratio > best_score:
                best_score = overlap_ratio
                best_match = idx

        if best_match != -1:
            ref_to_lines[best_match].append(line)
            matched_flags[i] = True

    # Merge matched lines by x-order
    merged_lines = []
    for matched_group in ref_to_lines.values():
        all_pts = [pt for line in matched_group for pt in line]
        unique_pts = sorted(set(map(tuple, all_pts)), key=lambda p: p[0])
        merged_lines.append(unique_pts)

    # Add unmatched lines directly
    unmatched_lines = [lines[i] for i, matched in enumerate(matched_flags) if not matched]

    # Combine and return all
    return merged_lines  #+ unmatched_lines


def tiled_and_tsp_merge(lines, window_size=50, overlap=1, angle_thresh = 50, threshold_distance = 50):
    """
    Tiles the space and applies modified_spatial3 per tile, then merges segments 
    back into original line structure.

    Parameters:
        lines: list of list of (x, y)
        modified_spatial3: function(lines_subset) → lines_subset_modified
        window_size: tiling size
        overlap: extra margin to include neighbors

    Returns:
        mod_lines: modified version of input lines (same length, merged)
    """

    all_points = [pt for line in lines for pt in line]
    if not all_points:
        return []

    # Compute bounds + margin frame
    x_coords = [pt[0] for pt in all_points]
    y_coords = [pt[1] for pt in all_points]
    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)

    frame_min_x = min_x - overlap
    frame_max_x = max_x + overlap
    frame_min_y = min_y - overlap
    frame_max_y = max_y + overlap

    # Storage to collect modified pieces by original line index
    line_segments = defaultdict(list)

    # Tiling over space
    x = frame_min_x
    while x < frame_max_x:
        y = frame_min_y
        while y < frame_max_y:
            x_start, x_end = x, x + window_size
            y_start, y_end = y, y + window_size

            # Extract segments inside this tile
            local_lines = []
            local_line_ids = []

            for i, line in enumerate(lines):
                segment = [
                    pt for pt in line
                    if (x_start - overlap <= pt[0] <= x_end + overlap) and
                       (y_start - overlap <= pt[1] <= y_end + overlap)
                ]
                if segment:
                    local_lines.append(segment)
                    local_line_ids.append(i)

            if local_lines:
                # Apply spatial transformation
                modified_segments = tsp_initiation(local_lines, angle_thresh , threshold_distance, False)
                local_lines = conjoin_lines(local_lines)
#                 modified_segments = local_lines  # test

                # Reassign segments back to their original line id
                for line_id, segment in zip(local_line_ids, modified_segments):
                    line_segments[line_id].extend(segment)

            y += window_size - overlap
        x += window_size - overlap

    # Reconstruct lines: sort merged segments by x
    mod_lines = []
    for i in range(len(lines)):
        merged = sorted(set(line_segments[i]), key=lambda pt: pt[0])
        mod_lines.append(merged)

    return mod_lines



#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#   |t|s|p|_|i|n|i|t|i|a|t|i|o|n|
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
def point_to_line_map(lines):
    # Flatten all points and map each to its source line index
    point_to_line_idx = {} # map new

    for idx, line in enumerate(lines):
        if not line:
            continue
        for pt in line:
            point_to_line_idx[pt] = idx
            
    return point_to_line_idx


def find_first_last_peak_or_trough(line, sigma=3):

    # Convert to numpy array
    line = np.asarray(line)
    y = line[:, 1]  # extract y-values

    # Smooth the signal
    y_smooth = gaussian_filter1d(y, sigma=sigma)

    # Find peaks and troughs
    peaks, _ = find_peaks(y_smooth)
    troughs, _ = find_peaks(-y_smooth)

    # Combine and sort indices
    extrema = np.sort(np.concatenate((peaks, troughs)))

    if len(extrema) == 0:
        return (None, None)
    elif len(extrema) == 1:
        return (int(extrema[0]), int(extrema[0]))
    else:
        return (int(extrema[0]), int(extrema[-1]))



def smoothed_angle(line, head=True, k=1):
    if len(line) < k + 1:
        return 0
    pts = line[:k+1] if head else line[-(k+1):]
    dx = pts[-1][0] - pts[0][0]
    dy = pts[-1][1] - pts[0][1]
    return np.arctan2(dy, dx)

def angle_between(p1, p2):   # always between 0 to 90 degree
    dx, dy = abs(p2[0] - p1[0]), abs(p2[1] - p1[1])  # take absolute differences
    if dx == 0 and dy == 0:
        return 0.0  # identical points
    return math.degrees(math.atan2(dy, dx))



def tsp_maps(lines):
    # Flatten all points and map each to its source line index
    point_to_line_idx = {} # map new
    line_idx_to_mean_dy = {} # map new  #(first_k, last_k)
    line_idx_to_k12 = {}

    for idx, line in enumerate(lines):
        if not line or len(line) < 2:
            continue
        mean_dy = sum(abs(line[i+1][1] - line[i][1]) for i in range(len(line)-1)) / max(len(line)-1, 1)  # dy is about mean of two consecutive stripe points
        line_idx_to_mean_dy[idx] = mean_dy
        k12 = find_first_last_peak_or_trough(line)
        line_idx_to_k12[idx] = k12
        for pt in line:
            pt_tuple = tuple(pt)
            point_to_line_idx[pt_tuple] = idx
            
    return point_to_line_idx, line_idx_to_mean_dy, line_idx_to_k12


def tsp_terminations(lines):
    terminations_start = [line[0] for i, line in enumerate(lines)]
    terminations_start_set = set(terminations_start)
    terminations_end = [line[-1] for i, line in enumerate(lines)]    
    terminations_end_set = set(terminations_end)
    terminations = terminations_end + terminations_start
    terminations = sorted(terminations, key=lambda t: t[0])
    
    return terminations_start, terminations_start_set, terminations_end, terminations_end_set, terminations


def calculate_x_delta(points):
    """
    Calculate the difference between max and min x-coordinates.
    
    Parameters:
    points: list of (x, y) tuples
    
    Returns:
    float: difference between max_x and min_x
    """
    x_coords = [x for x, y in points]
    return max(x_coords) - min(x_coords)


def tsp_initiation(lines, angle_thresh = 50, threshold_distance = 50, verbose = False):   # np.pi / 18 = 10 degree # revisited
    
    turtle_walk = [
    '🐢    ',
    ' 🐢   ',
    '  🐢  ',
    '   🐢 ',
    '    🐢',
    '   🐢 ',
    '  🐢  ',
    ' 🐢   ',
    ]
    if not lines:
        return []
    if isinstance(lines[0], tuple):  # Handle single line
        lines = [lines]
        
    #angle_thresh = math.radians(angle_thresh)  # no need for new angle_between function
    lines = [line.copy() for line in lines]
    coords_set = {coord for line in lines for coord in line}
    
    point_to_line_idx = point_to_line_map(lines)
    terminations_start, terminations_start_set, terminations_end, terminations_end_set, terminations = tsp_terminations(lines)
    
    total_points = len(terminations)
    processed_stripes = lines.copy()
    processed_points = set()
    tree = KDTree(terminations)

    max_outer_iterations = total_points  # Or another safe upper bound
    inner_step_count = 0
    for end_pt in terminations_end:
        inner_step_count += 1
        #clear_output(wait=True)  if not verbose else None # Clears entire cell output
        #print("end_pt",end_pt,"in iter", inner_step_count, "of",len(terminations_end), "ends")
        if end_pt not in terminations_end_set: # due to the fact terminations_end_set becomes updated when successful
            print(f"bypassed by updated ends set", current_point) if verbose else None
            continue

        current_line = [end_pt]
        success = False  # Flag to check if we connected to anything
        temp_processed = set(current_line)  # Track points only if connection is valid
        current_point = end_pt


        distances, indices = tree.query(current_point, k = total_points)

        found_valid = False
        kd_point_processed = set()
        for idx in np.atleast_1d(indices):
            kd_point = terminations[idx]
            
            
            if distances[idx] > threshold_distance:
                break
            
            
            if kd_point == current_point:
#                 if (72, 449) == current_point and (74, 446) == kd_point:
#                     print(f"bypassed by self_match")
                print(f"bypassed by self_match", current_point, kd_point) if verbose else None
                continue  

            # kd_point not being a start point
            if kd_point not in terminations_start_set:
#                 if (72, 449) == current_point and (74, 446) == kd_point:
#                     print(f"bypassed by start-to-end bypass")                
                print(f"bypassed by start-to-end bypass", current_point, kd_point) if verbose else None
                continue

            # kd_point and current_point belong to the same line index
            line_idx_kd = point_to_line_idx.get(kd_point, float('inf'))
            line_idx_current = point_to_line_idx.get(current_point, float('inf'))
            if line_idx_kd == line_idx_current or line_idx_kd == float('inf') or line_idx_current == float('inf'):
#                 if (72, 449) == current_point and (74, 446) == kd_point:
#                     print(f"bypassed by same line index")                
                print(f"bypassed by same line index", current_point, kd_point) if verbose else None
                continue

            # Enforce increasing x-direction
            curr_x = current_point[0]
            next_x = kd_point[0]
#            if (next_x - curr_x) <= 0:  
            if (next_x - curr_x) <= 0:  # new; giving tolerance to the length of kd line => + calculate_x_delta(processed_stripes[line_idx_kd]) => removed because it produces error
                temp_processed.add(kd_point)
#                 if (72, 449) == current_point and (74, 446) == kd_point:
#                     print(f"bypassed by x-reversals")                
                print(f"bypassed by x-reversals", current_point, kd_point) if verbose else None
                continue

#             x_kd, y_kd = kd_point  # no need
#             if x_kd == x_min or y_kd == y_min or x_kd == x_max or y_kd == y_max:
            
            
            angle_current = angle_between(current_point, kd_point)  # angles sound obsoleted due to errorneous lower angle of intersecting lines
            if abs(angle_current) > angle_thresh:
                temp_processed.add(kd_point)
#                 if (72, 449) == current_point and (74, 446) == kd_point:
#                     print(f"bypassed by angle")                
                print(f"bypassed by angle", current_point, kd_point) if verbose else None
                continue

            # Check Bresenham intersection
            intersection_break = False
            connection_line = [current_point, kd_point]  
            another_set = bresenham_line(connection_line)
            intersection_flag = False
            for pt in another_set:
                pt_ = (pt[0], pt[1]-1) # this is for bresenham intersection ladder effect
                if (pt in coords_set or pt_ in coords_set) and pt not in connection_line:  # this is intersection condition! noticed processed_stripes are all in bresenham line condition
                    print(f"bypassed by intersection", current_point, kd_point) if verbose else None
                    intersection_flag = True
                    break  # this break is so important, if you continue here, it only gets to other intersection point
            if intersection_flag:
                continue
        
        
            ######## 
            if kd_point in kd_point_processed: # if condition 1 or 2 has been flagged a kd_point for a piggy backed current_line 
                continue
        ###################### preprocess => no merging happens in preprocessing, 2 lines come in, 2 lines come out
            # tsp lateral condition 1: the lower angles are favored in lateral paths 
            point_to_line_idx = point_to_line_map(processed_stripes)
            if kd_point in processed_points:
                line_idx_kd = point_to_line_idx.get(kd_point)
                line_idx_current = point_to_line_idx.get(current_point)
                old_line_idx = line_idx_kd
                old_kd_line = processed_stripes[old_line_idx]
                for p,pt in enumerate(old_kd_line):
                    if pt == kd_point:
                        kd_point_idx_old = p
                        old_current_point = old_kd_line[p - 1]
                        break

#                     angle_current = angle_between(current_point, kd_point)  # angles sound obsoleted due to errorneous lower angle of intersecting lines
#                     angle_current_old = angle_between(old_current_point, kd_point)
                dist_current = euclidean_distance(current_point, kd_point)
                dist_current_old = euclidean_distance(old_current_point, kd_point)


                if  dist_current_old <= dist_current:  # lower angles are favored in lateral paths
                    continue # shooting kd_point to garbage bin for current_line !
                elif dist_current < dist_current_old:
                    # 1) erase wrong path from sets and records              
                    old_current_line = processed_stripes[old_line_idx]
                    processed_points.difference_update(old_current_line)          
                    # 2) mario jump from current_point to end of previous tsp! by replacing point up to kd_point
                    processed_stripes[old_line_idx] = current_line + old_current_line[kd_point_idx_old:]
                    # 3) add new current line to sets
                    processed_points.update(processed_stripes[old_line_idx])
                    point_to_line_idx = point_to_line_map(processed_stripes)  # update delegated points index to adopting line
                    # 4) piggy back old_current_line on current_line for saving computation for what was wrong
                    if kd_point_idx_old > 0:
                        current_line = old_current_line[:kd_point_idx_old]
                        kd_point_processed.add(kd_point)
                        found_valid = True  # This make sure while keeps going!
                        break
                    elif kd_point_idx_old == 0:
                        current_line = processed_stripes[old_line_idx]
                        kd_point_processed.add(kd_point)
                        found_valid = False # stop growing and find a new current_line
                        success = False # because we appended tsp lateral, we don't need to do it again
                        break        
        ##########################
        
        
        
        

            # Valid point found                     
            temp_processed.add(kd_point)
#             print("success") if verbose else None
            current_line = bresenham_line(connection_line)
            processed_points.update(connection_line)
            processed_stripes[line_idx_current].extend(processed_stripes[line_idx_kd])            
            processed_stripes[line_idx_current] = bresenham_line(processed_stripes[line_idx_current])
            processed_stripes[line_idx_kd] = []  # test

            
#             intersections = fast_five(processed_stripes)
#             if intersections:
#                 return print(f"Error#000: processed_stripes have unexpected intersections:", intersections)
            
            # Updating mappings because of conjoin function
            point_to_line_idx, line_idx_to_mean_dy, line_idx_to_k12 = tsp_maps(processed_stripes)
            
            # Updating termination parameters
            terminations_end_set.discard(current_point)
            terminations_start_set.discard(kd_point)
            terminations.remove(current_point)
            terminations.remove(kd_point)
            total_points = len(terminations)
            tree = KDTree(terminations)
            break
            
#             clear_output(wait=True)  if not verbose else None # Clears entire cell output
#             print(f"{len(processed_stripes)} from {total_points} {turtle_walk[inner_step_count % len(turtle_walk)]}")  if not verbose else None
         

#     print(f"\nNumber of processed_stripes: {len(processed_stripes)}")
    return [line for line in processed_stripes]# if line]


def tsp(lines, angle_thresh = 45, threshold_distance = 10, verbose = False):
    """Repeatedly apply tsp_initiation() until the length of the result stabilizes."""
    current_lines = lines
    previous_len = len(current_lines)
    
    while_iter = 0
    while True:
        while_iter += 1
        mod_lines = tsp_initiation(current_lines, angle_thresh, threshold_distance, verbose)  # Apply transformation
        current_len = len(mod_lines)
        
        # Stop if length hasn't changed
        if current_len == previous_len:
            break
        
        # Update for next iteration
        current_lines = mod_lines
        previous_len = current_len
#         print("while_iter",while_iter)
    
    return current_lines



def tsp_squares(squares, angle_threshold = 45):

    all_lines = []

    for row in range(squares.shape[0]):
        for col in range(squares.shape[1]):
            cell = squares[row, col]
            cell = tsp(cell, False,angle_threshold)
            if cell:  # not empty
                all_lines.extend(cell)

    return all_lines
    
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#   
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+    


#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#   |b|r|e|s|e|n|h|a|m|_|i|n|t|e|r|s|e|c|t|i|o|n|
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
def five_condition_intersection_finder(lines, bresenham = True):
    """
    Given a list of lines, where each line is a list of (x, y) tuple coordinates,
    this function returns all points that satisfy either:
      1) The point (x, y) lies on any other line besides its own.
      OR
      2) All below are true for the point (x, y) on current line:
         a) (x+1, y-1) lies on the same (current) line.
         b) (x, y-1) lies on any other line besides its own.
         c) (x+1, y) lies on any other line besides its own.
      OR
      3) All below are true for the point (x, y) on current line:
         a) (x-1, y-1) lies on the same (current) line.
         b) (x, y-1) lies on any other line besides its own.
         c) (x-1, y) lies on any other line besides its own.

    Invalid or malformed points (e.g., not a 2-tuple of ints) are skipped.
    """
    point_to_lines = defaultdict(set)
    lines = [line.copy() for line in lines if line]
    if bresenham:
        lines = bresenham_line(lines)
        #lines = [line for line in lines if line]
    cleaned_lines = []
    
    if lines is None:
        return []    

    for idx, raw_line in enumerate(lines):
        this_line_pts = set()
        if not isinstance(raw_line, list):
            cleaned_lines.append(this_line_pts)
            continue
        for pt in raw_line:
            try:
                x, y = pt
            except Exception:
                continue
            if not (isinstance(x, int) and isinstance(y, int)):
                continue
            this_line_pts.add((x, y))
            point_to_lines[(x, y)].add(idx)
        cleaned_lines.append(this_line_pts)

    intersection_points = set()

    for i, current_set in enumerate(cleaned_lines):
        if not current_set:
            continue
        for (x, y) in current_set:
            # Condition 1: (x, y) lies on another line
            if point_to_lines[(x, y)] - {i}:
                intersection_points.add((x, y))
                continue  # no need to evaluate the rest

            # Right-side condition group
            right_diag = (x + 1, y - 1) in current_set
            below = point_to_lines.get((x, y - 1), set()) - {i}
            right = point_to_lines.get((x + 1, y), set()) - {i}
            cond_group_2 = right_diag and bool(below) and bool(right)

            # Left-side condition group
            left_diag = (x - 1, y - 1) in current_set
            left = point_to_lines.get((x - 1, y), set()) - {i}
            cond_group_3 = left_diag and bool(below) and bool(left)

            if cond_group_2 or cond_group_3:
                intersection_points.add((x, y))

    return list(intersection_points)


def fast_five(lines):
    """
    Optimized version of five_condition_intersection_finder:
    Returns the first point satisfying the conditions and stops immediately.
    """
    point_to_lines = defaultdict(set)
    lines = [line.copy() for line in lines]
    lines = bresenham_line(lines)
    cleaned_lines = []

    for idx, raw_line in enumerate(lines):
        this_line_pts = set()
        if not isinstance(raw_line, list):
            cleaned_lines.append(this_line_pts)
            continue
        for pt in raw_line:
            try:
                x, y = pt
            except Exception:
                continue
            if not (isinstance(x, int) and isinstance(y, int)):
                continue
            this_line_pts.add((x, y))
            point_to_lines[(x, y)].add(idx)
        cleaned_lines.append(this_line_pts)

    for i, current_set in enumerate(cleaned_lines):
        if not current_set:
            continue
        for (x, y) in current_set:
            # Condition 1
            if point_to_lines[(x, y)] - {i}:
                return (x, y)

            # Check below once for reuse
            below = point_to_lines.get((x, y - 1), set()) - {i}

            # Condition 2
            if (x + 1, y - 1) in current_set and below and point_to_lines.get((x + 1, y), set()) - {i}:
                return (x, y)

            # Condition 3
            if (x - 1, y - 1) in current_set and below and point_to_lines.get((x - 1, y), set()) - {i}:
                return (x, y)

    return None  # No intersection found




def detect_frozen_intersections(intersection_history, current_intersections, frozen_threshold=3):
    """Detect if we're stuck with the same intersections recurring"""
    if len(intersection_history) < frozen_threshold:
        return False
    
    # Check if current intersections match any of last N historical sets
    current_set = frozenset(current_intersections)
    for historical_set in intersection_history[-frozen_threshold:]:
        if historical_set == current_set:
            return True
    return False
    

def remove_intersections_segmentation_hard(lines, intersection_points, min_segment_length=2):
    """Remove all segments containing intersections, keep only clean parts above min_segment_length"""
    intersection_set = set(intersection_points)
    cleaned_lines = []
    
    for line in lines:
        if not line:
            continue
        
        segments = []
        current_segment = []
        
        for pt in line:
            if pt in intersection_set:
                if len(current_segment) >= min_segment_length:
                    segments.append(current_segment)
                current_segment = []
            else:
                current_segment.append(pt)
        
        if len(current_segment) >= min_segment_length:
            segments.append(current_segment)
        
        cleaned_lines.extend(segments)  # Add all valid segments
    
    return cleaned_lines


    

def remove_intersections_segmentation_new(lines):
    """Aggressively remove segments containing intersections, keeping only largest clean segments"""
    flag_point = -999
    max_iterations = 100
    frozen_threshold = 3
    intersection_history = []
    
    lines = [line.copy() for line in lines]
    lines = bresenham_line(lines)
    
    while_iter = 0
    
    while while_iter < max_iterations:
        current_intersections = five_condition_intersection_finder(lines)
#         current_intersections, longest_lines_indices = five_condition_intersection_enhanced(lines)
        if not current_intersections:
            break
        
        ### 
        # Check for frozen state
        frozen = detect_frozen_intersections(
            intersection_history, 
            current_intersections,
            frozen_threshold
        )
        
        if frozen:
            print(f"Frozen intersections detected - applying aggressive segmentation, iter:", while_iter)
            lines = remove_intersections_segmentation_hard(lines, current_intersections)
            lines = bresenham_line(lines)
            intersection_history = []  # Reset history after treatment
            continue        

        # Track intersection history
        intersection_history.append(frozenset(current_intersections))
        if len(intersection_history) > frozen_threshold:
            intersection_history.pop(0)            
            
            
        ###    
            
        intersection_set = set(current_intersections)
        cleaned_lines = []

        for idx, line in enumerate(lines):
            if not line:
                continue
                
#             if idx in longest_lines_indices:
#                 cleaned_lines.append(line) # shorcut for the longest lines in intersections
#                 continue

            # Find all segments not containing intersection points
            segments = []
            current_segment = []

            for pt in line:
                if pt in intersection_set:
                    if current_segment and len(current_segment) >= 2:  # Finalize current segment
                        segments.append(current_segment)
                        current_segment = []
                else:
                    current_segment.append(pt)

            if current_segment and len(current_segment) >= 2:  # Add last segment if exists
                segments.append(current_segment)

            if segments:
                cleaned_lines += segments
                
        lines = bresenham_line(cleaned_lines)
        while_iter += 1
        
        
    final_intersections = five_condition_intersection_finder(lines)
    if final_intersections:
        print("Final cleanup: removing lingering intersections")
        lines = remove_intersections_segmentation_hard(lines, final_intersections)
        lines = bresenham_line(lines)        
    

    return lines

#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def generate_diagonal_offset_layers(n_layers=3):
    """
    Generate offsets for rasterized layers around the center point (0, 0),
    excluding the center. Each layer includes points at Manhattan distance = layer index.
    """
    offsets = [(0,0)]
    for d in range(1, n_layers + 1):
        for dx in range(-d, d + 1):
            dy = d - abs(dx)
            offsets.append((dx, dy))
            if dy != 0:
                offsets.append((dx, -dy))
    return offsets

def generate_square_offset_layers(n_layers=3):
    """
    Generate square-shaped offset layers around the center (0, 0),
    excluding the center. Each layer is a square ring at Chebyshev distance `d`.
    """
    offsets = [(0,0)]
    for d in range(1, n_layers + 1):
        for dx in range(-d, d + 1):
            for dy in range(-d, d + 1):
                if abs(dx) == d or abs(dy) == d:  # edge of the square
                    offsets.append((dx, dy))
    return offsets



def tsp_point(lines, angle_threshold = 45):
    turtle_walk = [
    '🐢    ',
    ' 🐢   ',
    '  🐢  ',
    '   🐢 ',
    '    🐢',
    '   🐢 ',
    '  🐢  ',
    ' 🐢   ',
    ]
    if not lines:
        return []
    
    angle_threshold = math.radians(angle_threshold)
              
    if isinstance(lines[0], tuple):  # Handle single line
        lines = [lines]
        
    lines = deepcopy(lines)
    # Flatten all points and map each to its source line index
    flat_points = [] # temporary parameter
    
    
    for idx, line in enumerate(lines):
        if not line:
            continue
        for pt in line:
            pt_tuple = tuple(pt)
            flat_points.append(pt_tuple)
    
    flat_points = list(set(flat_points))  # so the rest of gpt 5, looks exactly the same except for delta_y bypass
    flat_points.sort(key=lambda p: (p[0], p[1]))
    total_points = len(flat_points)
    processed_lines = []
    processed_points = set()
    attempted_starts = set()
    tree = KDTree(flat_points)
    node_length = 20


    max_outer_iter_count = total_points  # Or another safe upper bound
    outer_iter_count = 0
    while len(processed_points) < total_points and outer_iter_count < max_outer_iter_count:
        outer_iter_count += 1
        
        # Find a new valid start point
        start_candidates = [p for p in flat_points if p not in processed_points and p not in attempted_starts]
        if not start_candidates:
            break

        start_point = min(start_candidates, key=lambda p: (p[0], p[1]))
        current_line = [start_point]

        success = False  # Flag to check if we connected to anything
        attempted_starts.add(start_point)

        x_range = max(p[0] for p in flat_points) - min(p[0] for p in flat_points)
        max_inner_iter_count = int(x_range)    
        inner_iter_count = 0
        reverse = False  # start tsp with normal (record) direction
        record_current_line = False
        kd_point_processed = set()
        while inner_iter_count < max_inner_iter_count:  # current_line growth
            inner_iter_count += 1
            
            
            
            if len(flat_points)  == len(processed_points):
                return processed_lines
                
            
#             if len(current_line) % node_length == 0:  # reverse condition (tsp second condition)
#                 reverse = not reverse
                
            if reverse and not record_current_line:
                record_current_line = current_line.copy()
                record_current_line_set = set(record_current_line)
                current_line = [current_line[-1]]  # restart current_line
            
            if record_current_line and not reverse:
                current_line = record_current_line
                record_current_line = False
            
            
            current_point = current_line[-1]

            distances, indices = tree.query(current_point, k = total_points)

            found_valid = False
            for idx in np.atleast_1d(indices):
                kd_point = flat_points[idx]
                
                if kd_point == current_point:
                    continue  # Skip self
                

                # Enforce increasing x-direction
                if not reverse:
                    if len(current_line) >= 2:
                        curr_x = current_line[-1][0]
                        next_x = kd_point[0]
                        if (next_x - curr_x) < 0:
#                             if (72, 449) == current_point and (74, 446) == kd_point:
#                                 print(f"bypassed by x-reversals")
                            continue 
                elif reverse:
                    if len(current_line) >= 2:
                        curr_x = current_line[-1][0]
                        next_x = kd_point[0]
                        if (next_x - curr_x) > 0:
    #                         print(f"bypassed by x-reversals")
                            continue 
        
                angle_current = angle_between(current_point, kd_point)  # angles sound obsoleted due to errorneous lower angle of intersecting lines
                if abs(angle_current) > angle_threshold:
#                     if (72, 449) == current_point and (74, 446) == kd_point:
#                         print(f"bypassed by angle")
                    continue 

         
                ######## 
                if kd_point in kd_point_processed: # if condition 1 or 2 has been flagged a kd_point for a piggy backed current_line 
                    continue
                    
                #######
                if reverse: 
                    if kd_point not in record_current_line_set: # crash
                        for p, pt in enumerate(record_current_line):
                            if pt == current_point:
                                false_record_kd_point_idx = p
                                shouldhave_kd_point = record_current_line[p - 1] # notice we traverse normal not reverse in record_current_line
                                break
                                
                        angle_record = angle_between(shouldhave_kd_point, current_point) # direction of angles are normal so not to cause trouble
                        angle_reverse = angle_between(kd_point, current_point)
                        
                        if  angle_record <= angle_reverse:  # lower angles are favored in lateral paths
                            continue # give reverse a chance to find normal tsp
                        elif angle_reverse < angle_record:  # no record of reverse current_line but rerun of kd_point for normal current_line
                            if false_record_kd_point_idx > 0:
                                current_line = record_current_line[:false_record_kd_point_idx]
                                kd_point_processed.add(current_point)  # reverse crashed current_point is kd_point that should be skipped
                                current_point = shouldhave_kd_point
                                reverse = False
                                found_valid = True  # This make sure while keeps going!
                                break
                            elif false_record_kd_point_idx == 0:
                                current_line = [start_point]
                                kd_point_processed.add(current_point)
                                reverse = False
                                found_valid = False # stop growing and find a new current_line
                                success = False # because we appended tsp lateral, we don't need to do it again
                                break        
                        
                    
                
                # tsp lateral condition 1: the lower angles are favored in lateral paths 
                if kd_point in processed_points:
                    line_idx_kd = point_to_line_idx.get(kd_point, float('inf'))
                    line_idx_current = point_to_line_idx.get(current_point, float('inf'))
                    old_line_idx = point_to_line_idx.get(kd_point, None)
                    old_kd_line = processed_lines[old_line_idx]
                    for p,pt in enumerate(old_kd_line):
                        if pt == kd_point:
                            kd_point_idx_old = p
                            old_current_point = old_kd_line[p - 1]
                            break
                            
#                     angle_current = angle_between(current_point, kd_point)  # angles sound obsoleted due to errorneous lower angle of intersecting lines
#                     angle_current_old = angle_between(old_current_point, kd_point)
                    dist_current = euclidean_distance(current_point, kd_point)
                    dist_current_old = euclidean_distance(old_current_point, kd_point)

                    
                    if  dist_current_old <= dist_current:  # lower angles are favored in lateral paths
                        continue # shooting kd_point to garbage bin for current_line !
                    elif dist_current < dist_current_old:
                        # 1) erase wrong path from sets and records              
                        old_current_line = processed_lines[old_line_idx]
                        processed_points.difference_update(old_current_line)          
                        # 2) mario jump from current_point to end of previous tsp! by replacing point up to kd_point
                        processed_lines[old_line_idx] = current_line + old_current_line[kd_point_idx_old:]
                        # 3) add new current line to sets
                        processed_points.update(processed_lines[old_line_idx])
                        point_to_line_idx = point_to_line_map(processed_lines)  # update delegated points index to adopting line
                        # 4) piggy back old_current_line on current_line for saving computation for what was wrong
                        if kd_point_idx_old > 0:
                            current_line = old_current_line[:kd_point_idx_old]
                            kd_point_processed.add(kd_point)
                            found_valid = True  # This make sure while keeps going!
                            break
                        elif kd_point_idx_old == 0:
                            current_line = processed_lines[old_line_idx]
                            kd_point_processed.add(kd_point)
                            found_valid = False # stop growing and find a new current_line
                            success = False # because we appended tsp lateral, we don't need to do it again
                            break
                           
                ########
    
                #Check Bresenham intersection   ==> very important, intersections should be after reprocess
                connection_line = [current_point, kd_point]           
                
                test_lines = processed_lines + [connection_line]         
                if fast_five(test_lines):  # new intersections!
#                     if (72, 449) == current_point and (74, 446) == kd_point:
#                         print(f"bypassed by intersection", current_point, kd_point)
                    continue      
    
    
                    
                # Valid point found  # kd for loop continuation
                current_line.append(kd_point)  
                
                
                #   +-+-+-+-+-+-+-+-+-+
                #   |a|n|i|m|a|t|i|o|n|
                #   +-+-+-+-+-+-+-+-+-+
#                 if inner_iter_count % 10:
#                     test_lines = processed_lines + [current_line]
#                     plot_layers(test_lines, pixels = 256, saving_directory="frames", marked_layers = test_lines[-1], close_option = True)



                found_valid = True  # This make sure while keeps going!
                success = True
                break  # kd_point for loop break

            if not found_valid and not reverse: # not breaking in reverse
                break  # break while loop as tsp is finished for start_point
            elif not found_valid and reverse:
                reverse = not reverse

        if success:  # success in finishing while loop and current_line growth
            processed_lines.append(current_line)
            processed_points.update(current_line)
            point_to_line_idx = point_to_line_map(processed_lines)  # new
            
            #   +-+-+-+-+-+-+-+-+-+
            #   |a|n|i|m|a|t|i|o|n|
            #   +-+-+-+-+-+-+-+-+-+
#             if outer_iter_count % 10:
#                 plot_layers(processed_lines, pixels = 256, saving_directory="frames", marked_layers = processed_lines[-1], close_option = True)
            
    
    
#             clear_output(wait=True)  # Clears entire cell output
#             print(f"{len(processed_lines)} from {total_points} {turtle_walk[outer_iter_count % len(turtle_walk)]}")


    print(f"\nNumber of processed flat_points: {len(processed_lines)}")
    return processed_lines



def tsp_point_recursive(lines, angle_threshold = 45):
    """Repeatedly apply tsp_point() until there is no intersections."""
            
    current_lines = lines
    previous_len = len(current_lines)
    
    while_iter = 0
    while True:
        while_iter += 1
        mod_lines = tsp_point(current_lines, angle_threshold)  # Apply transformation
        current_len = len(mod_lines)
        
        # Stop if length hasn't changed
        if current_len == previous_len:
            break
        
        # Stop if there are no intersections
        if not fast_five(mod_lines):
            return mod_lines
        
        # Update for next iteration
        current_lines = mod_lines
        previous_len = current_len
#         print("while_iter",while_iter)
    
    return current_lines            
            
            
            
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#   |m|o|d|i|f|i|e|d|_|s|p|a|t|i|a|l|3|
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+            
            
def find_nearest_point(end_pt, intersections):
    if not intersections:
        return None

    # Filter out the end_pt itself
    candidates = [pt for pt in intersections if pt != end_pt]
    if not candidates:
        return None  # No other points to compare

    nearest = min(candidates, key=lambda pt: (pt[0] - end_pt[0])**2 + (pt[1] - end_pt[1])**2)
    return nearest


def point_to_line_map_dict(lines):
    point_to_line_idx = defaultdict(set)
    for idx, line in enumerate(lines):
        if not line:
            continue
        for pt in line:
            point_to_line_idx[pt].add(idx)
    return dict(point_to_line_idx)



def remove_consecutive_duplicate_points(lines):
    if isinstance(lines[0], tuple):  # Handle single line
        lines = [lines]
    
    cleaned_lines = []
    for line in lines:
        if not line:  # Skip empty lines
            cleaned_lines.append(line)
            continue
            
        # Keep first point, then check consecutive duplicates
        unique_points = [line[0]]
        for i in range(1, len(line)):
            if line[i] != line[i-1]:  # Only add if different from previous
                unique_points.append(line[i])
        cleaned_lines.append(unique_points)
    
    return cleaned_lines if len(cleaned_lines) > 1 else cleaned_lines[0]


def bresenham_line(lines):
    if not lines:
        return []
    # If input is a single line (list of tuples/lists), wrap it in a list
    if isinstance(lines[0], tuple):
        lines = [lines]

    def bresenham_segment(p1, p2):
        x1, y1 = map(int, tuple(p1))  # ensure p1 is tuple
        x2, y2 = map(int, tuple(p2))
        points = []
        dx = abs(x2 - x1)
        dy = abs(y2 - y1)
        sx = 1 if x1 < x2 else -1
        sy = 1 if y1 < y2 else -1
        err = dx - dy

        while True:
            points.append((x1, y1))
            if x1 == x2 and y1 == y2:
                break
            e2 = err * 2
            if e2 > -dy:
                err -= dy
                x1 += sx
            if e2 < dx:
                err += dx
                y1 += sy

        return points

    refined_lines = []
    for line in lines:
        # Normalize all points to tuples
        line = [tuple(p) for p in line]
        refined_line = []
        for i in range(len(line) - 1):
            refined_line.extend(bresenham_segment(line[i], line[i + 1]))
        refined_lines.append(refined_line)

    refined_lines = remove_consecutive_duplicate_points(refined_lines)
    
    if not refined_lines:
        return None

    if isinstance(refined_lines[0], tuple):  # this is needed because of consecutive function above
        refined_lines = [refined_lines]    
    return refined_lines  if len(refined_lines) > 1 else refined_lines[0]


def find_all_intersections(lines_):
    
    lines = lines_.copy()
    lines = bresenham_line(lines)
    
    if isinstance(lines[0], tuple):
        lines = [lines]
    
#     print(lines)
    point_counts = {}
    for line in lines:
        seen_in_line = set()
        for point in line:
            # This avoids counting duplicates within a single line multiple times
#             print(point)
            if point not in seen_in_line:
                point_counts[point] = point_counts.get(point, 0) + 1
                seen_in_line.add(point)

    intersections = [point for point, count in point_counts.items() if count > 1]
    return intersections



def modified_spatial3_ends(lines, global_frames = False, angle_threshold = 50, grid = 50):
    
    turtle_walk = [
    '🐢    ',
    ' 🐢   ',
    '  🐢  ',
    '   🐢 ',
    '    🐢',
    '   🐢 ',
    '  🐢  ',
    ' 🐢   ',
    ]
    if not lines:
        return []
    if isinstance(lines[0], tuple):  # Handle single line
        lines = [lines]
#     lines = deepcopy(lines)    
    
    # Flatten all points and keep track of which line they come from
    all_points = []
    line_indices = []
    for i, line in enumerate(lines):
        for pt in line:
            all_points.append(pt)
            line_indices.append(i)

    angle_threshold = math.radians(angle_threshold)
            
    #
#     lines = bresenham_line(lines)
    
    # Calculate x_min, x_max, y_min, y_max from the given lines
    x_min = min(x for x, _ in all_points)
    x_max = max(x for x, _ in all_points)
    y_min = min(y for _, y in all_points)
    y_max = max(y for _, y in all_points)
    
    minmax_set = set([x_min, x_max, y_min, y_max])
    # Making bresenham frame lines:
    f1 = [(x_min, y_min), (x_min, y_max)]  # left frame
    f2 = [(x_min, y_max), (x_max, y_max)]
    f3 = [(x_max, y_max), (x_max, y_min)]  # right frame
    f4 = [(x_max, y_min), (x_min, y_min)] 
    
    frame_lines = []
    
#     if global_frames:
    
#         f1_g, f2_g, f3_g, f4_g = global_frames
#         F = 0
#         if any(f1) in f1_g:
#             F += 1
#             frame_lines.append(f1)
#             lines.append(bresenham_line(f1))                
#         if any(f2) in f2_g:
#             F += 1            
#             frame_lines.append(f2)
#             lines.append(bresenham_line(f2))   
#         if any(f3) in f3_g:
#             F += 1            
#             frame_lines.append(f3)
#             lines.append(bresenham_line(f3))         
#         if any(f4) in f4_g:
#             F += 1            
#             frame_lines.append(f4)
#             lines.append(bresenham_line(f4))
#         #
#     elif not global_frames:
    lines.append(bresenham_line(f1))                
    lines.append(bresenham_line(f2))    
    lines.append(bresenham_line(f3))         
    lines.append(bresenham_line(f4))        
        
    
            
    # map
#     point_to_line_idx = point_to_line_map_dict(lines)
    # Build KDTree from all points
    tree = KDTree(all_points)
    mod_lines = [line.copy() for line in lines]  # deep copy to avoid mutation
    lines_to_remove = set()
    
    for i, line in enumerate(lines):
        if not line or line[-1][0] in minmax_set:
            continue  # skip empty lines
        end_pt = line[-1]

       
        start_pt = line[0]
        x_end = end_pt[0]
        x_start = start_pt[0]

#         current_line_indices =  point_to_line_idx[end_pt]
        
        
        # Query multiple nearest neighbors in case we need to skip some [end points]
#         dists, idxs = tree.query(end_pt, k=min(grid,len(all_points)))
        dists, idxs = tree.query(start_pt, k= len(all_points))        
#         mask = dists < grid  # Replace 50 with your desired distance threshold
#         idxs = idxs[mask]        

        for turtle_count, idx in enumerate(idxs):
        
            target_pt = all_points[idx]
#             target_line_index = line_indices[idx]  # right!?
#             target_line_indices =  point_to_line_idx[target_pt]
#             if len(target_line_indices) > 1: # shared point
#                 continue
#             else:
#                 target_line_index = list(target_line_indices)[0]
                   
            if target_pt[0] <= x_end:
                #print(f"bypassed by x direction")
                continue  # skip if x is not greater
                
            angle_current = angle_between(end_pt, target_pt)  # angles sound obsoleted due to errorneous lower angle of intersecting lines
            if abs(angle_current) > angle_threshold:
                #print(f"bypassed by angle")
                continue    
               
                
            if turtle_count == (len(idxs)-1):  # avoiding stuck target_pt
                #print(f"turtle break")
                break
                
#             connection_line = [end_pt, target_pt]
#             connection_line = bresenham_line(connection_line)
#             # Build mod_lines excluding current line
#             mod_lines_excl = mod_lines[:i] + mod_lines[i+1:]

#             # Check for intersection
#             mod_lines_excl.append(connection_line)
#             intersections = five_condition_intersection_finder(mod_lines_excl)
#             if intersections and target_pt in intersections: # the reason this works, is because there is no way in bresenham lines, two intersections happen for the nearest point
#                 nearest = find_nearest_point(end_pt, intersections)
            
#             if end_pt == (738, 676):
#                 return print(f"line 0 wants to append {target_pt} to {end_pt}")
            mod_lines[i].append(target_pt)
            mod_lines[i] = bresenham_line(mod_lines[i])
        
        ### while adding connection_line to all_points is right, it slows the method with lowest modification
#             connection_line = bresenham_line([end_pt, target_pt])
#             internal_points = connection_line[1:-1]
#             if internal_points: all_points.extend(internal_points) # order doesn't matter 
#             tree = KDTree(all_points)        
        ###
            #clear_output(wait=True)  # Clears entire cell output
            #print(f"ends: {i} from {len(lines)} {turtle_walk[turtle_count % len(turtle_walk)]}")                
            break  # bypassing intersection. why? I believe the nearest point is a shared point 
                     # and cannot diagonally intercepts lines, because there is no way that it 
                     # passes through the line and a shared point on the other line wouldn't be detected as target_pt    
                                
#                 break  # connection added, move to next line


##################
#     print("removable indices:", lines_to_remove)
    # Remove all merged lines from mod_lines
    mod_lines = [line for idx, line in enumerate(mod_lines)]
    return mod_lines

def modified_spatial3_starts(lines, global_frames = False, angle_threshold = 50, grid = 50):
    
    turtle_walk = [
    '🐢    ',
    ' 🐢   ',
    '  🐢  ',
    '   🐢 ',
    '    🐢',
    '   🐢 ',
    '  🐢  ',
    ' 🐢   ',
    ]
    if not lines:
        return []
    if isinstance(lines[0], tuple):  # Handle single line
        lines = [lines]
#     lines = deepcopy(lines)    
    
    # Flatten all points and keep track of which line they come from
    all_points = []
    line_indices = []
    for i, line in enumerate(lines):
        for pt in line:
            all_points.append(pt)
            line_indices.append(i)

    angle_threshold = math.radians(angle_threshold)
            
    #
#     lines = bresenham_line(lines)
    
    # Calculate x_min, x_max, y_min, y_max from the given lines
    x_min = min(x for x, _ in all_points)
    x_max = max(x for x, _ in all_points)
    y_min = min(y for _, y in all_points)
    y_max = max(y for _, y in all_points)
    
    minmax_set = set([x_min, x_max, y_min, y_max])    
    # Making bresenham frame lines:
    f1 = [(x_min, y_min), (x_min, y_max)]  # left frame
    f2 = [(x_min, y_max), (x_max, y_max)]
    f3 = [(x_max, y_max), (x_max, y_min)]  # right frame
    f4 = [(x_max, y_min), (x_min, y_min)] 
    
#     frame_lines = []
    
#     if global_frames:
    
#         f1_g, f2_g, f3_g, f4_g = global_frames
#         F = 0
#         if any(f1) in f1_g:
#             F += 1
#             frame_lines.append(f1)
#             lines.append(bresenham_line(f1))                
#         if any(f2) in f2_g:
#             F += 1            
#             frame_lines.append(f2)
#             lines.append(bresenham_line(f2))    
#         if any(f3) in f3_g:
#             F += 1            
#             frame_lines.append(f3)
#             lines.append(bresenham_line(f3))         
#         if any(f4) in f4_g:
#             F += 1            
#             frame_lines.append(f4)
#             lines.append(bresenham_line(f4))
#         #
#     elif not global_frames:
    lines.append(bresenham_line(f1))                
    lines.append(bresenham_line(f2))    
    lines.append(bresenham_line(f3))         
    lines.append(bresenham_line(f4))        
        
    
            
    # map
#     point_to_line_idx = point_to_line_map_dict(lines)
    # Build KDTree from all points
    tree = KDTree(all_points)
    mod_lines = [line.copy() for line in lines]  # deep copy to avoid mutation
    lines_to_remove = set()
    
    
    for i, line in enumerate(lines):
        if not line or line[0][0] in minmax_set:
            continue  # skip empty lines
        end_pt = line[-1]
        start_pt = line[0]
        x_end = end_pt[0]
        x_start = start_pt[0]    
        
#         current_line_indices =  point_to_line_idx[start_pt]
        
        # Query multiple nearest neighbors in case we need to skip some [start points]                
#         dists, idxs = tree.query(start_pt, k= min(grid,len(all_points)))
        dists, idxs = tree.query(start_pt, k= len(all_points))        
#         mask = dists < grid  # Replace 50 with your desired distance threshold
#         idxs = idxs[mask]        
        
        for turtle_count, idx in enumerate(idxs):
            target_pt = all_points[idx]
#             target_line_index = line_indices[idx]  # wrong? below is right?

#             target_line_indices =  point_to_line_idx[target_pt]  # it is alright to have a shared point at target_pt
#             if len(target_line_indices) > 1: # shared point
#                 continue
#             else:
#                 target_line_index = list(target_line_indices)[0]
                
            if target_pt[0] >= x_start:
                #print("bypassed by x direction")
                continue  # skip if x is not lower

            angle_current = angle_between(target_pt, start_pt)  # angles sound obsoleted due to errorneous lower angle of intersecting lines
            if abs(angle_current) > angle_threshold:
                #print("bypassed by angle")                
                continue    
    
            if turtle_count == (len(idxs)-1):  # avoiding stuck target_pt:
                #print("turtle break")
                break
                
#             connection_line = [start_pt, target_pt]
#             connection_line = bresenham_line(connection_line)
#             # Build mod_lines excluding current line
#             mod_lines_excl = mod_lines[:i] + mod_lines[i+1:]

#             # Check for intersection
#             mod_lines_excl.append(connection_line)
#             intersections = five_condition_intersection_finder(mod_lines_excl)
#             if intersections and target_pt in intersections: # the reason this works, is because there is no way in bresenham lines, two intersections happen for the nearest point
#                 nearest = find_nearest_point(start_pt, intersections)

#             if turtle_count < (len(idxs)-1):  # avoiding stuck target_pt
#                 mod_lines[i].extend(connection_line)

            mod_lines[i].insert(0, target_pt)
            mod_lines[i] = bresenham_line(mod_lines[i])
        ### while adding connection_line to all_points is right, it slows the method with lowest modification        
#             connection_line = bresenham_line([target_pt, start_pt])
#             internal_points = connection_line[1:-1]
#             if internal_points: all_points.extend(internal_points) # order doesn't matter 
#             tree = KDTree(all_points)
        ###
        
#             clear_output(wait=True)  # Clears entire cell output
#             print(f"starts: {i} from {len(lines)} {turtle_walk[turtle_count % len(turtle_walk)]}")                
            break  # bypassing intersection. why? I believe the nearest point is a shared point 
                     # and cannot diagonally intercepts lines, because there is no way that it 
                     # passes through the line and a shared point on the other line wouldn't be detected as target_pt    
                                
#                 break  # connection added, move to next line


##################
#     print("removable indices:", lines_to_remove)
    # Remove all merged lines from mod_lines
    mod_lines = [line for idx, line in enumerate(mod_lines)]
    return mod_lines


def conjoin_lines(lines):
    """
    Given a list of lines (each line is a list of tuple coordinates), this function repeatedly finds and merges
    pairs of lines where the end of one line matches the start of another. Merging always preserves the order:
        merged_line = line_with_shared_end + line_with_shared_start
    After merging, the longer of the two original lines is replaced by the merged line, and the shorter line
    is removed. The process repeats until no further merges are possible.

    Parameters
    ----------
    lines : list of list of tuple
        A list where each element represents a line, itself a list of (x, y) coordinate tuples.

    Returns
    -------
    list of list of tuple
        A list of all new lines that were formed by merging (in the order they were created).
        The original `lines` list is modified in-place: merged lines replace one of the originals,
        and the other original is removed.

    Notes
    -----
    - If a line is empty or not a list of tuples, it is ignored in the merging logic.
    - Direction matters: we only merge if the end point of one line exactly equals the start point of another.
    - After each merge, mappings of start/end points are rebuilt so that multiple chained merges can occur.
    """
    
    lines = [line.copy() for line in lines]
    conjoined_lines = []

    # Helper to build mappings from point -> list of line-indices
    def build_maps(current_lines):
        start_map = defaultdict(list)  # point -> [indices where it's a start]
        end_map = defaultdict(list)    # point -> [indices where it's an end]
        for idx, ln in enumerate(current_lines):
            # skip invalid or empty lines
            if not isinstance(ln, list) or len(ln) == 0:
                continue
            start_pt = ln[0]
            end_pt = ln[-1]
            # only consider tuple-like points
            if isinstance(start_pt, tuple):
                start_map[start_pt].append(idx)
            if isinstance(end_pt, tuple):
                end_map[end_pt].append(idx)
        return start_map, end_map

    # Repeat until no merges happen in a pass
    while True:
        start_map, end_map = build_maps(lines)
        merged_this_round = False

        # Find any matching end -> start pair
        for shared_pt, end_indices in end_map.items():
            if shared_pt not in start_map:
                continue
            start_indices = start_map[shared_pt]

            # Try all combinations of (i, j) where i ends at shared_pt and j starts at shared_pt
            for i in end_indices:
                for j in start_indices:
                    if i == j:
                        # skip the trivial case where a line's end == its own start
                        continue

                    # We have line_i ending at shared_pt, and line_j starting at shared_pt
                    line_i = lines[i]
                    line_j = lines[j]
                    # Form the merged sequence: end of i meets start of j
                    merged = list(line_i) + list(line_j)

                    # Decide which index to replace based on length
                    # The longer original line will be replaced by `merged`; the shorter will be removed
                    if len(line_j) > len(line_i):
                        replace_idx = j
                        remove_idx = i
                    else:
                        replace_idx = i
                        remove_idx = j

                    # Update lines list: replace & remove
                    lines[replace_idx] = merged
                    # Always remove the line with the smaller index last to preserve correct indices
                    # But we want to pop in a way that does not disturb replace_idx.
                    if remove_idx < replace_idx:
                        lines.pop(remove_idx)
                        replace_idx -= 1  # shifting because an earlier item was removed
                    else:
                        lines.pop(remove_idx)

                    # Record the newly formed line
                    conjoined_lines.append(merged)

                    # Mark that we did a merge
                    merged_this_round = True
                    #clear_output(wait=True)  # Clears entire cell output
                    #print("connected start to end:", i, "to", j)

                    # Break out to rebuild mappings and search again from scratch
                    break
                if merged_this_round:
                    break
            if merged_this_round:
                break

        # If no merges in this iteration, we're done
        if not merged_this_round:
            break

    # Combine and remove empty/duplicate lines
    seen = set()
    result = []
    for line in conjoined_lines + lines:
        if not line:  # skip empty lines
            continue
        key = tuple(line)  # convert to tuple for hashable comparison
        if key not in seen:
            seen.add(key)
            result.append(line)
            
    return result




def modified_spatial3(lines, angle_threshold = 50, grid = 50):
    all_points = []
    for i, line in enumerate(lines):
        for pt in line:
            all_points.append(pt)    
    
    # Calculate x_min, x_max, y_min, y_max from the given lines
    x_min = min(x for x, _ in all_points)
    x_max = max(x for x, _ in all_points)
    y_min = min(y for _, y in all_points)
    y_max = max(y for _, y in all_points)
    
    
    # Making bresenham frame lines:
    f1 = [(x_min, y_min), (x_min, y_max)]  # left frame
    f2 = [(x_min, y_max), (x_max, y_max)]
    f3 = [(x_max, y_max), (x_max, y_min)]  # right frame
    f4 = [(x_max, y_min), (x_min, y_min)] 
    
    global_frames = [f1, f2, f3, f4]
    
    lines.append(bresenham_line(f2))    
    lines.append(bresenham_line(f4))
    #
    lines.append(bresenham_line(f3)) 
    lines.append(bresenham_line(f1))
    
    spatial3_result_ends = modified_spatial3_ends(lines, global_frames, angle_threshold, grid)
    
    spatial3_result_starts = modified_spatial3_starts(spatial3_result_ends[:-4], global_frames, angle_threshold, grid)    
    return spatial3_result_starts[:-4]

    
def tsp_squares_03(squares, angle_threshold = 45, grid = 50):

    all_lines = []

    for row in range(squares.shape[0]):
        for col in range(squares.shape[1]):
            cell = squares[row, col]
            cell = modified_spatial3(cell, angle_threshold, grid)
            if cell:  # not empty
                all_lines.extend(cell)

    return all_lines



def tiled_and_modified_spatial3_merge(lines, window_size=100, overlap=20):
    """
    Tiles the space and applies modified_spatial3 per tile, then merges segments 
    back into original line structure.

    Parameters:
        lines: list of list of (x, y)
        modified_spatial3: function(lines_subset) → lines_subset_modified
        window_size: tiling size
        overlap: extra margin to include neighbors

    Returns:
        mod_lines: modified version of input lines (same length, merged)
    """

    all_points = [pt for line in lines for pt in line]
    if not all_points:
        return []

    # Compute bounds + margin frame
    x_coords = [pt[0] for pt in all_points]
    y_coords = [pt[1] for pt in all_points]
    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)

    frame_min_x = min_x - overlap
    frame_max_x = max_x + overlap
    frame_min_y = min_y - overlap
    frame_max_y = max_y + overlap

    # Making bresenham frame lines:
    f1 = [(frame_min_x, frame_min_y), (frame_min_x, frame_max_y)]  # left frame
    f2 = [(frame_min_x, frame_max_y), (frame_max_x, frame_max_y)]
    f3 = [(frame_max_x, frame_max_y), (frame_max_x, frame_min_y)]  # right frame
    f4 = [(frame_max_x, frame_min_y), (frame_min_x, frame_min_y)] 
    global_frames = [f1, f2, f3, f4]
    # Storage to collect modified pieces by original line index
    line_segments = defaultdict(list)

    # Tiling over space
    x = frame_min_x
    while x < frame_max_x:
        y = frame_min_y
        while y < frame_max_y:
            x_start, x_end = x, x + window_size
            y_start, y_end = y, y + window_size

            # Extract segments inside this tile
            local_lines = []
            local_line_ids = []

            for i, line in enumerate(lines):
                segment = [
                    pt for pt in line
                    if (x_start - overlap <= pt[0] <= x_end + overlap) and
                       (y_start - overlap <= pt[1] <= y_end + overlap)
                ]
                if segment:
                    local_lines.append(segment)
                    local_line_ids.append(i)

            if local_lines:
                # Apply spatial transformation
#                 modified_segments = local_lines
                #
                modified_segments = modified_spatial3_ends(local_lines, global_frames, 50, 100)
                modified_segments = modified_spatial3_starts(modified_segments[:-4], global_frames, 50, 100)[:-4]
#                 modified_segments = filter_dead_lines(modified_segments, threshold=3)
                
                # Reassign segments back to their original line id
                for line_id, segment in zip(local_line_ids, modified_segments):
                    line_segments[line_id].extend(segment)

            y += window_size - overlap
        x += window_size - overlap

    # Reconstruct lines: sort merged segments by x
    mod_lines = []
    for i in range(len(lines)):
        merged = sorted(set(line_segments[i]), key=lambda pt: pt[0])
        mod_lines.append(merged)

    return mod_lines

#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#   
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


def remove_frame_lines(lines):
    all_points = []
    for i, line in enumerate(lines):
        for pt in line:
            all_points.append(pt)    
    
    # Calculate x_min, x_max, y_min, y_max from the given lines
    x_min = min(x for x, _ in all_points)
    x_max = max(x for x, _ in all_points)
    y_min = min(y for _, y in all_points)
    y_max = max(y for _, y in all_points)

    # Define the frame lines
    f1 = [(x_min, y_min), (x_min, y_max)]
    f2 = [(x_min, y_min), (x_max, y_min)]
    f3 = [(x_min, y_max), (x_max, y_max)]
    f4 = [(x_max, y_min), (x_max, y_max)]

    frame_lines = [f1, f2, f3, f4]
    frame_lines = bresenham_line(frame_lines)

    # Filter out any lines that match one of the frame lines
    filtered_lines = [line for line in lines if line not in frame_lines]

    return filtered_lines



def filter_lines_by_angle(lines, threshold_angle=70):
    """
    this helps to remove frame lines and more
    """
    filtered_lines = []
    for line in lines:
        if not line or len(line) < 2:
            continue  # skip empty or too short lines

        (x1, y1), (x2, y2) = line[0], line[-1]
        dx = x2 - x1
        dy = y2 - y1

        angle_rad = math.atan2(dy, dx)
        angle_deg = abs(math.degrees(angle_rad))

        # Ensure angle is in [0, 90] range since slope direction doesn't matter
        angle_from_horizontal = min(angle_deg, 180 - angle_deg)

        if angle_from_horizontal <= threshold_angle:
            filtered_lines.append(line)

    return filtered_lines


def filter_lines_by_bresenham_growth(lines, threshold=3):
    """
    another effective function to filter frame lines and mistakes
    """
    filtered_lines = []
    for line in lines:
        if len(line)>1:
            if len(bresenham_line(line)) <= len(line) + threshold:
                filtered_lines.append(line)
    return filtered_lines



def filter_dead_lines(lines, threshold=5):
    """
    Filters out lines that do not show any slope change within the first `threshold` segments.
    """
    def angle_between_2(p1, p2):
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        return math.atan2(dy, dx)  # in radians

    filtered_lines = []

    for line in lines:
        if len(line) < threshold + 1:
            continue  # skip lines that are too short to evaluate

        is_alive = False
        prev_angle = angle_between_2(line[0], line[1])

        for i in range(1, min(threshold, len(line) - 1)):
            angle = angle_between_2(line[i], line[i + 1])
            if abs(angle - prev_angle) > 1e-5:  # detect any change
                is_alive = True
                break
            prev_angle = angle

        if is_alive:
            filtered_lines.append(line)

    return filtered_lines



def flip_lines_horizontally(lines):

    if not lines:
        return []

    # Find the global maximum x
    all_x = [pt[0] for line in lines for pt in line if pt is not None and len(pt) == 2]
    if not all_x:
        return []

    max_x = max(all_x)

    # Flip each x across the max_x line
    flipped_lines = [
        [(max_x - x, y) for (x, y) in line]
        for line in lines
    ]

    return flipped_lines



def exclude_border_lines(lines):
    """
    Excludes lines that:
    1) Contain at least two points on the global x/y min/max values.
    2) Have a slope equal to a frame line: vertical or horizontal.
    """
    # Step 1: Get global min/max values
    all_points = [pt for line in lines for pt in line if line]
    if not all_points:
        return []

    x_values = [x for x, _ in all_points]
    y_values = [y for _, y in all_points]

    x_min, x_max = min(x_values), max(x_values)
    y_min, y_max = min(y_values), max(y_values)

    minmax_set = {x_min, x_max, y_min, y_max}
    filtered_lines = []

    for line in lines:
        if len(line) < 2:
            filtered_lines.append(line)
            continue

        # --- Condition 1: count how many points fall on minmax x or y ---
        count = 0
        for x, y in line:
            if x in (x_min, x_max) or y in (y_min, y_max):
                count += 1
                if count >= 2:
                    break

        # --- Condition 2: check slope of line (start to end) ---
        (x1, y1), (x2, y2) = line[0], line[-1]
        dx = x2 - x1
        dy = y2 - y1

        is_horizontal = dy == 0
        is_vertical = dx == 0

        # Exclude line only if both conditions met
        if count >= 2 and (is_horizontal or is_vertical):
            continue  # exclude this line

        filtered_lines.append(line)

    return filtered_lines




#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#   |c|o|n|n|e|c|t|_|l|i|n|e|s|_|l|a|t|e|r|a|l|l|y|
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def is_shared_point(lines, point):
    count = 0
    for line in lines:
        if point in line:
            count += 1
            if count >= 2:
                return True
    return False


def find_non_shared_or_boundary_lines(lines):
    """
    Finds indices of lines that either:
    1. Have no shared start/end point, and
    2. Have no start/end point that lies on the global x/y boundary.
    """
    if not lines:
        return []

    # Compute global x/y min/max
    all_points = [pt for line in lines for pt in line if line]
    x_min = min(x for x, _ in all_points)
    x_max = max(x for x, _ in all_points)
    y_min = min(y for _, y in all_points)
    y_max = max(y for _, y in all_points)

    def is_on_boundary(pt):
        x, y = pt
        return x in (x_min, x_max) or y in (y_min, y_max)

    # Build map from point -> set of line indices it's used in
    point_to_lines = defaultdict(set)
    for idx, line in enumerate(lines):
        for pt in line:
            point_to_lines[pt].add(idx)

    result_indices = []
    for idx, line in enumerate(lines):
        if not line:
            continue
        start, end = line[0], line[-1]

        start_shared = len(point_to_lines[start]) >= 2
        end_shared = len(point_to_lines[end]) >= 2

        start_bad = not start_shared and not is_on_boundary(start)
        end_bad = not end_shared and not is_on_boundary(end)

        if start_bad or end_bad:
            result_indices.append(idx)

    return result_indices



def connect_lines_laterally_ends(lines):
    """
    Iteratively extend each line by traversing shared endpoints until no line grows further.
    Repeats the extension process for all lines until a full pass yields zero growth.

    Args:
        lines (list of list of tuple):
            Initial list of lines, where each line is a list of (x, y) tuples in strictly increasing x order.
            Each line’s endpoints intersect another line, unless they lie on a global min/max of x or y.

    Returns:
        list of list of tuple:
            Final list of extended lines after no further growth is possible.
    """
    skip_lines = set() 
    # Make a working copy of lines
    current_lines = [list(ln) for ln in lines]

    while True:
        # 1) Build point→[line_indices] mapping for this iteration
        point_to_lines = {}
        for idx, ln in enumerate(current_lines):
            for pt in ln:
                point_to_lines.setdefault(pt, []).append(idx)

        # 2) Compute global maxima in x and y for all points in current_lines
        all_x = [pt[0] for ln in current_lines for pt in ln]
        all_y = [pt[1] for ln in current_lines for pt in ln]
        global_max_x = max(all_x)
        global_min_x = min(all_x)
        global_max_y = max(all_y)
        global_min_y = min(all_y)

        new_lines = []
        any_growth = False

        # 3) Try to extend each line
        for idx, ln in enumerate(current_lines):
            extended = list(ln)  # copy of this line
            while True:
                endpt = extended[-1]

                # Stop if endpoint is at a global max
                if (
                    endpt[0] == global_max_x 
                    or endpt[0] == global_min_x 
                    or endpt[1] == global_max_y 
                    or endpt[1] == global_min_y
                ):
                    break

                mapped_idxs = point_to_lines.get(endpt, [])
                if len(mapped_idxs) == 1:
                    # Only itself shares this point
                    skip_lines.add(idx)
                    print(f"line {idx} shows no shared point on its end")
                    break

                if len(mapped_idxs) == 2:
                    # Exactly two lines share this endpoint: pick the “other” one
                    other_idx = mapped_idxs[0] if mapped_idxs[1] == idx else mapped_idxs[1]
                    other_ln = current_lines[other_idx]

                    # If endpoint is also the other line's endpoint, stop
                    if endpt == other_ln[0] or endpt == other_ln[-1]:
                        break

                    # Otherwise, traverse from shared point to other line's end
                    pos_in_other = other_ln.index(endpt)
                    segment = other_ln[pos_in_other + 1 :]
                    extended.extend(segment)
                    continue

                # More than two lines share this endpoint
                candidates = []
                for j in mapped_idxs:
                    ln_j = current_lines[j]
                    if (
                        endpt not in (ln_j[0], ln_j[-1])
                        and endpt[0] != global_max_x
                        and endpt[0] != global_min_x
                        and endpt[1] != global_max_y
                        and endpt[1] != global_min_y
                    ):
                        candidates.append(j)

                if not candidates:
                    # No valid “through-lines” left: stop
                    break

                # Pick the longest candidate line
                other_idx = max(candidates, key=lambda j: len(current_lines[j]))
                other_ln = current_lines[other_idx]
                pos_in_other = other_ln.index(endpt)
                segment = other_ln[pos_in_other + 1 :]
                extended.extend(segment)
                # Loop again with new endpoint

            new_lines.append(extended)
            if len(extended) > len(ln):
                any_growth = True

        # 4) If no line grew in this pass, we are done
        if not any_growth:
            break

        # Otherwise, use these newly extended lines for the next iteration
        current_lines = new_lines

    return   [line for i,line in enumerate(current_lines) if i not in skip_lines]   # current_lines


def flip_lines_by_x_axis(lines):
    flipped = []
    for ln in lines:
        # Reverse the original line (so that after negating x, the x-values are increasing)
        rev = reversed(ln)
        # Negate x in each point
        flipped_line = [(-x, y) for x, y in rev]
        flipped.append(flipped_line)
    return flipped

def connect_lines_laterally(R):
    R = connect_lines_laterally_ends(R)
    R = flip_lines_by_x_axis(R)
    R = connect_lines_laterally_ends(R)
    R = flip_lines_by_x_axis(R)
    return R

# => maximum downlaps function


def filter_frame_lines(lines):
    """
    Given a list of lines (each line is a list of (x, y) tuples), compute the four
    bounding-box corner points responsible for frame lines, then return a new list of lines
    excluding any line that contains at least two of those corner points.
    """
    if not lines:
        return [], []

    # Flatten and separate coordinates
    xs, ys = [], []
    for line in lines:
        for pt in line:
            try:
                x, y = pt
            except (TypeError, ValueError):
                raise ValueError(f"Invalid point {pt!r}; expected a tuple (x, y).")
            xs.append(x)
            ys.append(y)

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    bounding_points = [
        (min_x, min_y),
        (min_x, max_y),
        (max_x, min_y),
        (max_x, max_y),
    ]

    # Build filtered list: keep only lines with < 2 bounding points
    filtered_lines = []
    for line in lines:
        # count how many corner points lie on this line
        count_corners = sum(1 for pt in line if pt in bounding_points)
        if count_corners < 2:
            filtered_lines.append(line)

    return filtered_lines



def maximum_downlaps(lines, M):
    '''estimates list of lines as the most of downlapped lines in the profile.
    this technique is used to extract sequence geometries from seismic layers
    M = number of desired maximums'''
    
    filtered_lines = []
    downlaps = np.zeros((len(lines)))
    for i, line1 in enumerate(lines):  # downlapped line
        if line1 == []:
            continue
        for j, line2 in enumerate(lines): # candidate line to downlap its start or end on line2
            if line2 == []:
                continue
            if i != j:
#                 for point2 in line2:
                start_shared = False
                end_shared = False
                start_point, end_point = line2[0], line2[-1]
#                 Check if start point or end point is shared with another line
                if start_point in line1 and (start_point != line1[0] or start_point != line1[-1]): # notice a maximum downlapped line on one side (up or below) has maximum start downlaps and on the other side end downlaps
                    start_shared = True
                if end_point in line1 and (end_point != line1[0] or end_point != line1[-1]):
                    end_shared = True
#                     if point2 in line1:
                    # Break early if both start and end are shared  (very nice chatgpt!)
                if start_shared or end_shared:  # or: one downlap, and: two downlaps
                    downlaps[i] += 1  # i is for line1 (calculates for all lines, even zeros show no if) 
        # after creating number of downlaps ... >        
        downlaps[i] =  len(line1)  #* downlaps[i]   # both length (notice because lines are rasterized, euclidean length is very close to bresenham length (simple count of pixels)) and num of down/on laps mark a distinct feature
    
    I = np.flip(np.argsort(downlaps)[-M:]) # max downlaps indices
    print("score for length * number of downlaps", np.flip(np.sort(downlaps)[-M:]))  # flip is for making order by maximum on first
    maxdownlaps = []
    i = 0  # predefine i in case the loop doesn't run
    for i in I:
        if downlaps[i] != 0:  # zero marks len(lines) = 1 because num of down/on laps for such line is always zero
            maxdownlaps.append(lines[i])
        else:
            break
    
    maxdownlaps = filter_frame_lines(maxdownlaps)
    return maxdownlaps, I[:i]  # notice how when breaks i shows last i [TG!!]!


def extract_feature(lines, boundaries):
    """
    Extracts feature lines that are between the top_line and bottom_line derived from the boundaries.

    Args:
        lines (list of list of tuple): List of (x, y) coordinate lines.
        boundaries (list of list of tuple): Subset of lines used to determine top and bottom bounding lines.

    Returns:
        tuple:
            bounding_lines (tuple of list of tuple): (top_line, bottom_line)
            filtered_lines (list of list of tuple): Lines between top and bottom lines.
    """
    if not lines or not boundaries:
        return (None, None), []

    # Corrected definition: top_line = max y, bottom_line = min y
    bottom_line = min(boundaries, key=lambda line: line[0][1])
    top_line = max(boundaries, key=lambda line: line[0][1])

    top_y = top_line[0][1]
    bottom_y = bottom_line[0][1]

    filtered = []
    for line in lines:
        line_y = line[0][1]
        if top_line == bottom_line:
            if line_y < bottom_y:
                filtered.append(line)
        else:
            if bottom_y < line_y < top_y:
                filtered.append(line)

    return (top_line, bottom_line), filtered


#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#   
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#   |m|e|r|g|e|_|s|i|m|i|l|a|r|_|l|i|n|e|s|
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def merge_similar_lines_keep_longest(lines, threshold=0.9):
    """
    Group similar lines based on threshold similarity, 
    and return only the longest line from each group.

    Parameters
    ----------
    lines : list[list[tuple]]
        Input lines, each a list of (x, y) points.
    threshold : float, default=0.9
        Merge lines if overlap ratio >= threshold.

    Returns
    -------
    list[list[tuple]]
        One longest line from each similarity group.
    """
    if not lines:
        return []

    # Convert each line into a set for fast overlap checks
    line_sets = [set(line) for line in lines]
    used = [False] * len(lines)
    groups = []

    for i in range(len(lines)):
        if used[i]:
            continue

        group = [lines[i]]
        used[i] = True
        changed = True

        while changed:
            changed = False
            for j in range(len(lines)):
                if used[j] or i == j:
                    continue

                intersection = line_sets[i] & line_sets[j]
                ratio = len(intersection) / min(len(line_sets[i]), len(line_sets[j]))

                if ratio >= threshold:
                    group.append(lines[j])
                    used[j] = True
                    changed = True

        groups.append(group)

    # From each group, keep only the longest line
    result = []
    for group in groups:
        longest = max(group, key=len)
        result.append(longest)

    return result


def select_unique_lines_by_exact_start_end(lines):
    def to_tuple(pt):
        return tuple(pt)

    seen_pairs = {}
    selected_lines = []

    for idx, line in enumerate(lines):
        if not line:
            continue
        start = to_tuple(line[0])
        end = to_tuple(line[-1])
        key = (start, end)

        # Only add the first occurrence of a start-end pair
        if key not in seen_pairs:
            seen_pairs[key] = idx
            selected_lines.append(line)

    return selected_lines


def order_interpreted_lines(lines):
    """
    Process lines so that each line is adjusted such that no point in a selected base line
    is above any other line at the same x position.
    """

    # Step 1: Map each point to its line index
    point_to_line = {}
    for i, line in enumerate(lines):
        for pt in line:
            point_to_line[pt] = i

    # Step 2: Group points with similar x values
    x_groups = defaultdict(list)
    for i, line in enumerate(lines):
        for x, y in line:
            x_groups[x].append((i, x, y))  # Store line index for easy trace

    # Keep track of processed base lines
    processed = set()

    # Step 3: Define base line detection function
    def detect_baseline(lines, processed, x_groups):
        """
        Return the index of the line that has the most number of lowest y-values
        across shared x values, excluding already processed lines.
        """
        score = [0] * len(lines)
        for x_val, entries in x_groups.items():
            # Group entries by x and find the lowest y at that x
            min_y = min(y for _, _, y in entries)
            for idx, _, y in entries:
                if y == min_y and idx not in processed:
                    score[idx] += 1
        # Select the line with the maximum score
        best_idx = np.argmax(score)
        return best_idx if best_idx not in processed else None

    # Step 4: Define correction function
    def correct_baseline(lines, base_idx, x_groups):
        """
        Modify the baseline in lines such that its y-values at shared x are
        not above other lines.
        """
        corrected = []
        base_line = lines[base_idx]
        for x, y in base_line:
            # For the current x, find the min y of all lines at that x
            min_y = y
            for i, x_i, y_i in x_groups.get(x, []):
                if i != base_idx:
                    min_y = min(min_y, y_i)
            corrected.append((x, min_y))
        lines[base_idx] = corrected

    # Step 5–7: Loop until all lines are corrected
    total_lines = len(lines)
    while len(processed) < total_lines:
        base_idx = detect_baseline(lines, processed, x_groups)
        if base_idx is None:
            break  # No more unprocessed lines
        correct_baseline(lines, base_idx, x_groups)
        processed.add(base_idx)

    # Step 8: Return modified lines
    return lines



def merge_similar_excluding_exact_start_end(lines, threshold=0.9): # this function uses start/end match function and similarity merging function 
    """
    Combines:
    1. Protected lines based on exact (start, end) match
    2. Remaining lines merged by point-wise similarity

    Parameters:
        lines (list of list of tuple): Input lines
        threshold (float): Point similarity threshold

    Returns:
        list of list of tuple: Final result = protected + merged non-protected
    """

    def to_tuple(pt):
        return tuple(pt)

    # STEP 1: Find protected lines using exact start-end match
    protected_lines = select_unique_lines_by_exact_start_end(lines)
    protected_set = set()
    for line in protected_lines:
        if not line:
            continue
        key = (to_tuple(line[0]), to_tuple(line[-1]))
        protected_set.add(key)

    # STEP 2: Separate protected and unprotected lines
    unprotected_lines = []
    for line in lines:
        if not line:
            continue
        key = (to_tuple(line[0]), to_tuple(line[-1]))
        if key not in protected_set:
            unprotected_lines.append(line)

    # STEP 3: Run similarity merge on unprotected lines
    def merge_similar_lines(lines, threshold=0.9):
        lines = [list(line) for line in lines]
        line_sets = [set(map(to_tuple, line)) for line in lines]
        used = [False] * len(lines)
        merged_lines = []

        for i in range(len(lines)):
            if used[i]:
                continue
            current_line = set(line_sets[i])
            used[i] = True
            changed = True
            while changed:
                changed = False
                for j in range(len(lines)):
                    if used[j] or i == j:
                        continue
                    intersection = current_line & line_sets[j]
                    ratio = len(intersection) / min(len(current_line), len(line_sets[j]))
                    if ratio >= threshold:
                        current_line.update(line_sets[j])
                        used[j] = True
                        changed = True
            merged_lines.append(list(current_line))
        return merged_lines

    merged_unprotected = merge_similar_lines(unprotected_lines, threshold=threshold)

    # STEP 4: Combine results
    return protected_lines + merged_unprotected


#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#   
#   +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#   +-+-+-+-+-+-+-+-+-+-+-+-+
#   |s|m|o|o|t|h|_|l|i|n|e|s|
#   +-+-+-+-+-+-+-+-+-+-+-+-+

# x-group sort level # Thank God!! !! ~ finally, without global first point sorting, now x-group sorting made smooth function to be fine!

def smooth_lines(lines, x_tol=0.0, sigma=2.0):
    """
    Smooth lines without disturbing vertical order based on local x-groupings.

    Args:
        lines: List of lines, where each line is a list of (x, y) tuples.
        x_tol: Tolerance for grouping similar x-values.
        sigma: Smoothing intensity.

    Returns:
        Smoothed list of lines with the same structure (original input order).
    """
    if not lines:
        return []

    # --- no global sorting, keep original ---
    sorted_indices = list(range(len(lines)))
    sorted_lines = lines

    # --- flatten sorted lines and group by x (within x_tol) ---
    all_points = []
    for line_idx, line in enumerate(sorted_lines):
        for pt_idx, (x, y) in enumerate(line):
            all_points.append((x, y, line_idx, pt_idx))

    all_points.sort(key=lambda p: p[0])
    x_groups = []
    current = []
    last_x = None
    for x, y, l_idx, p_idx in all_points:
        if last_x is None or abs(x - last_x) <= x_tol:
            current.append((x, y, l_idx, p_idx))
        else:
            x_groups.append(current)
            current = [(x, y, l_idx, p_idx)]
        last_x = x
    if current:
        x_groups.append(current)

    # --- compute deltas per (line_idx, pt_idx), with local vertical sorting ---
    delta_map = defaultdict(float)
    for group in x_groups:
        # sort points vertically (by y) within this x-group
        group_sorted = sorted(group, key=lambda p: p[1])
        y_vals = [y for _, y, _, _ in group_sorted]
        smoothed_y = gaussian_filter1d(y_vals, sigma=sigma, mode='nearest')
        for (x, y, li, pj), smy in zip(group_sorted, smoothed_y):
            delta_map[(li, pj)] = smy - y

    # --- apply deltas to each line ---
    smoothed_sorted = []
    for line_idx, line in enumerate(sorted_lines):
        new_line = []
        for pt_idx, (x, y) in enumerate(line):
            d = delta_map.get((line_idx, pt_idx), 0.0)
            new_line.append((x, y + d))
        smoothed_sorted.append(new_line)

    # --- restore original order (trivial here, but kept for clarity) ---
    n = len(lines)
    smoothed = [None] * n
    for sorted_pos, orig_idx in enumerate(sorted_indices):
        smoothed[orig_idx] = smoothed_sorted[sorted_pos]

    return smoothed





def compute_crossing_scores(lines, intersections):
    """
    Compute crossing scores for lines based on shared intersection points.

    Args:
        lines (list of list of (x, y)): The set of polylines.
        intersections (list of (x, y)): The points where lines intersect.

    Returns:
        sorted_indices (list of int): Indices of lines with non-zero crossing score, sorted by score descending.
        sorted_lines (list of list of (x, y)): The lines corresponding to sorted indices.
        sorted_scores (list of int): Crossing scores corresponding to the sorted_indices.
    """
    
    def find_point_index(line, point):
        try:
            return line.index(point)
        except ValueError:
            return -1

    point_to_lines = defaultdict(list)
    for idx, line in enumerate(lines):
        line_set = set(line)
        for pt in intersections:
            if pt in line_set:
                point_to_lines[pt].append(idx)

    crossing_scores = defaultdict(int)

    for pt, shared_line_indices in point_to_lines.items():
        if len(shared_line_indices) < 2:
            continue

        for i_idx, j_idx in itertools.combinations(shared_line_indices, 2):
            line_i = lines[i_idx]
            line_j = lines[j_idx]

            i_pos = find_point_index(line_i, pt)
            j_pos = find_point_index(line_j, pt)

            if i_pos == -1 or j_pos == -1:
                continue

            def get_y_neighbors(line, pos):
                neg_y = line[pos - 1][1] if pos > 0 else None
                pos_y = line[pos + 1][1] if pos + 1 < len(line) else None
                return neg_y, pos_y

            i_neg_y, i_pos_y = get_y_neighbors(line_i, i_pos)
            j_neg_y, j_pos_y = get_y_neighbors(line_j, j_pos)

            if None in [i_neg_y, i_pos_y, j_neg_y, j_pos_y]:
                continue

            def relative_order(y1, y2):
                if y1 < y2:
                    return "above"
                elif y1 > y2:
                    return "below"
                else:
                    return "equal"

            neg_order = relative_order(i_neg_y, j_neg_y)
            pos_order = relative_order(i_pos_y, j_pos_y)

            if "equal" in [neg_order, pos_order]:
                continue

            if neg_order != pos_order:
                crossing_scores[i_idx] += 1
                crossing_scores[j_idx] += 1

    nonzero_scores = [(idx, score) for idx, score in crossing_scores.items() if score > 0]
    nonzero_scores.sort(key=lambda x: -x[1])

    sorted_indices = [idx for idx, _ in nonzero_scores]
    sorted_lines = [lines[idx] for idx in sorted_indices]
    sorted_scores = [score for _, score in nonzero_scores]

    return sorted_indices, sorted_lines, sorted_scores


def remove_crossing_lines_until_threshold(lines, intersections, N=10, score_threshold=10):
    """
    Iteratively remove top N lines with highest crossing scores until the max score ≤ threshold.

    Args:
        lines (list of list of (x, y)): List of polyline segments.
        intersections (list of (x, y)): List of intersection points.
        N (int): Number of lines to remove per iteration.
        score_threshold (int): Max allowed crossing score to stop iteration.

    Returns:
        remaining_lines (list of list of (x, y)): Lines left after removals.
        removed_indices (list of int): Indices from original lines that were removed.
        final_scores (list of int): Crossing scores of the remaining lines.
    """

    # Track original indices
    indexed_lines = list(enumerate(copy.deepcopy(lines)))
    removed_indices = []

    while True:
        current_lines = [line for _, line in indexed_lines]

        # Get crossing scores for current set
        sorted_indices, _, sorted_scores = compute_crossing_scores(current_lines, intersections)

        if not sorted_scores:
            break  # No crossings left
        if sorted_scores[0] <= score_threshold:
            break  # Crossing level acceptable

        # Determine which actual original indices to remove
        to_remove = set(sorted_indices[:N])

        # Record their original indices and remove them
        updated_lines = []
        for i, (original_index, line) in enumerate(indexed_lines):
            if i in to_remove:
                removed_indices.append(original_index)
            else:
                updated_lines.append((original_index, line))

        indexed_lines = updated_lines

        if not indexed_lines:
            break  # All lines removed

    # Final output
    remaining_lines = [line for _, line in indexed_lines]
    final_indices, _, final_scores = compute_crossing_scores(remaining_lines, intersections)

    return remaining_lines, removed_indices, final_scores


#   +-+-+-+-+-+-+-+-+-+-+-+-+
#   
#   +-+-+-+-+-+-+-+-+-+-+-+-+            


def filter_lines_on_global_bounds(lines, T=5):
    if not lines:
        return []

    # Flatten all points to find global min/max for x and y
    all_points = [pt for line in lines for pt in line]
    all_x = [pt[0] for pt in all_points]
    all_y = [pt[1] for pt in all_points]

    min_x, max_x = min(all_x), max(all_x)
    min_y, max_y = min(all_y), max(all_y)

    def is_on_boundary(pt):
        x, y = pt
        return (
            abs(x - min_x) <= T or
            abs(x - max_x) <= T or
            abs(y - min_y) <= T or
            abs(y - max_y) <= T
        )

    # Keep lines if BOTH start AND end points are on (tolerant) boundary
    filtered = [line for line in lines if is_on_boundary(line[0]) and is_on_boundary(line[-1])]

    return filtered


def filter_lines_between_maxdownlaps(lines, maxdownlaps, threshold=0.9):
    """
    Filters out lines that fall above the top line or below the bottom line
    in more than `threshold` proportion of their points.

    Args:
        lines (list of list of (x, y)): Input lines to filter
        maxdownlaps (list of list of (x, y)): Lines used to define top and bottom bounds
        threshold (float): Proportion of points outside bounds needed to exclude a line

    Returns:
        list: Filtered lines
    """
    if not maxdownlaps:
        return []

    # Identify top and bottom lines by start y-value
    start_ys = [(line[0][1], idx) for idx, line in enumerate(maxdownlaps) if line]
    top_idx = max(start_ys)[1]
    bottom_idx = min(start_ys)[1]

    top_line = maxdownlaps[top_idx]
    bottom_line = maxdownlaps[bottom_idx]

    # Create x→y maps for top and bottom
    top_map = {x: y for x, y in top_line}
    bottom_map = {x: y for x, y in bottom_line}

    # Find common x-values to compare against
    valid_x = set(top_map.keys()) & set(bottom_map.keys())

    def is_out_of_bounds(x, y):
        if x not in valid_x:
            return False  # Cannot determine, so treat as safe
        return y > top_map[x] or y < bottom_map[x]

    # Filter lines
    filtered = []
    for line in lines:
        if not line:
            continue
        out_count = sum(1 for x, y in line if is_out_of_bounds(x, y))
        if out_count / len(line) <= threshold:
            filtered.append(line)

    return filtered



def filter_lower_sharedpoints(lines):
    """
    For each line, remove duplicate points from consecutive groups with the same x-value.
    Special rule: If any point in the group appears as an endpoint (start or end) in any line,
    then preserve all such endpoint instances (in order) in that line.
    Otherwise, choose the point in that group that maps to the most lines.
    
    Args:
        lines (list of list of tuple): Each line is a list of (x, y) coordinates.
        
    Returns:
        list of list of tuple: New lines with redundant points removed.
    """
    # Build a set of global endpoints (first and last point of every line)
    global_endpoints = set()
    for line in lines:
        if line:
            global_endpoints.add(line[0])
            global_endpoints.add(line[-1])
    
    # Build mapping: point -> set of line indices where it appears.
    point_to_lines = defaultdict(set)
    for i, line in enumerate(lines):
        for pt in line:
            point_to_lines[pt].add(i)
    
    new_lines = []
    for line in lines:
        if not line:
            new_lines.append([])
            continue
        
        filtered_line = []
        i = 0
        n = len(line)
        # Process the line in consecutive groups of points sharing the same x.
        while i < n:
            run = [line[i]]
            j = i + 1
            while j < n and line[j][0] == line[i][0]:
                run.append(line[j])
                j += 1
            
            # If any point in the run is a global endpoint, preserve all unique occurrences in order.
            ge_points = []
            for pt in run:
                if pt in global_endpoints and pt not in ge_points:
                    ge_points.append(pt)
            if ge_points:
                # Append all global endpoints in this group.
                filtered_line.extend(ge_points)
            else:
                # No global endpoints in this run:
                # Choose the point with the highest mapping count.
                candidate = max(run, key=lambda pt: len(point_to_lines[pt]))
                filtered_line.append(candidate)
            
            i = j
        new_lines.append(filtered_line)
    
    return new_lines





def update_point_to_lines(point_to_lines, line, line_index):
    """Update the mapping of points to their corresponding lines."""
    for point in set(line):
        point_to_lines[point].add(line_index)

        
def parasitic_lines(lines_):
    lines = copy.deepcopy(lines_)
    point_to_lines = defaultdict(list)

    # 🛠️ Step 1: Ensure all lines are sorted by increasing x-values
    lines = sort_lines_by_x(lines)
    
    # Map each point to the lines they belong to
    point_to_lines = defaultdict(set)
    for i, line in enumerate(lines):
        update_point_to_lines(point_to_lines, line, i)
    
    x_max = max(point[0] for line in lines for point in line)  # Maximum x value
    x_min = min(point[0] for line in lines for point in line)  # Maximum x value
    
    for i in range(len(lines)):
        new_line = lines[i]  # Start with a copy of the original line
        first = min(new_line, key=lambda point: point[0])
        last =  max(new_line, key=lambda point: point[0])
        extended = True
        prev_length = -1  # Track previous length to avoid redundant loops
        prev_first = -1
        prev_last = -1
        iteration_count = 0  # Track how many times the while loop runs

        while (last[0] != x_max or first[0] != x_min) and iteration_count < 100:
            iteration_count += 1
            new_line = lines[i]  # Start with a copy of the original line
            first = min(new_line, key=lambda point: point[0])
            last =  max(new_line, key=lambda point: point[0])
            prev_length = len(new_line)
            extended = False

            

            if prev_first == first and prev_last == last:
                print(f"⚠️  No change in line {i} after {iteration_count} iterations! Exiting loop.")
                break 
                
            print(f"🔄 Iteration {iteration_count} for line {i} - First: {first}, Last: {last}")
                
            prev_first = first
            
            ##### for start            
            # Check if first point exists in other lines
            nn = []
            all_new_points = []
            indices = []
            for idx in point_to_lines[first]:
                if idx != i and len(lines[idx]) > 2:
                    other_line = lines[idx]
                    first_x = first[0]
                    new_points = [point for point in other_line if point[0] < first_x] # definitely there will be no repeats because x is different
                    try:
                        first_new_points =  min(new_points, key=lambda point: point[0])
                        nn.append( first_new_points[0] ) # x-values                    
                        indices.append(idx)
                        all_new_points.append(new_points)
                    except:
#                         print("first, line:", i)
                        pass
            
            if len(nn) > 0:
                I = np.argmin(nn)
                idx = indices[I]
                other_line = lines[idx]
                new_points = all_new_points[I]
                new_line = new_points + new_line
                update_point_to_lines(point_to_lines, new_points, i)  # 🔄 Update mapping # adding other line index to the point

                if len(new_line) > prev_length:
                    extended = True
                    print(f"📌 Extending line {i} backward using line {idx}")                    
                    lines[i] = new_line  # 🔥 Ensure memory is updated         

            ####  for last 
            prev_length = len(new_line)
            new_line = lines[i]
            last =  max(new_line, key=lambda point: point[0])
            prev_last = last
            # Check if last point exists in other lines
            nn = []
            all_new_points = []
            indices = []
            for idx in point_to_lines[last]:
                if idx != i and len(lines[idx]) > 2:
                    other_line = lines[idx]
                    last_x = last[0]
                    new_points = [point for point in other_line if point[0] > last_x]
                    try:
                        last_new_points =  max(new_points, key=lambda point: point[0])
                        nn.append( last_new_points[0] ) # x-values
                        indices.append(idx)
                        all_new_points.append(new_points)
                    except:
#                         print("last, line:", i)
                        pass     

            if len(nn) > 0:
                I = np.argmax(nn)
                idx = indices[I]
                other_line = lines[idx]
                new_points = all_new_points[I]
                new_line = new_line + new_points
                update_point_to_lines(point_to_lines, new_points, i)  # 🔄 Update mapping # adding other line index to the point

                if len(new_line) > prev_length:
                    extended = True
                    print(f"📌 Extending line {i} forward using line {idx}")                    
                    lines[i] = new_line  # 🔥 Ensure memory is updated

        
    
    return lines


def sort_lines_by_x(lines):  # used in spatial3  while loop
    """Ensure each line is sorted by increasing x-values."""
    return [sorted(line, key=lambda p: p[0]) for line in lines]
    
    
    
    
def correct_bad_lines(all_lines, bad_line_indices):
    """
    For each bad line, find the best pair of consecutive good lines (by index in all_lines)
    that cover the most x-points where the bad line y-value falls between them.
    Then correct the bad line where it violates the sandwich.

    Parameters:
        all_lines (list of list of (x, y)): All lines.
        bad_line_indices (list of int): Indices in all_lines considered as bad.

    Returns:
        corrected_all_lines: same structure as all_lines, with corrected bad lines.
        good_line_pairs: list of [lower_idx, upper_idx] used for correcting each bad line.
    """

    corrected_all_lines = copy.deepcopy(all_lines)

    bad_line_set = set(bad_line_indices)
    good_line_indices = [i for i in range(len(all_lines)) if i not in bad_line_set]

    # Group good line points by x for fast lookup
    x_to_y_line = defaultdict(list)
    for idx in good_line_indices:
        for (x, y) in all_lines[idx]:
            x_to_y_line[x].append((y, idx))

    # Sort y-values for each x
    for x in x_to_y_line:
        x_to_y_line[x].sort()

    result_pairs = []

    for bad_idx in bad_line_indices:
        bad_line = all_lines[bad_idx]
        x_to_y = dict(bad_line)

        pair_counter = defaultdict(int)

        for x, y in x_to_y.items():
            if x not in x_to_y_line:
                continue

            y_line_sorted = x_to_y_line[x]
            for i in range(len(y_line_sorted) - 1):
                y1, idx1 = y_line_sorted[i]
                y2, idx2 = y_line_sorted[i + 1]

                lower_y, upper_y = sorted([y1, y2])
                if lower_y <= y <= upper_y:
                    pair = tuple(sorted([idx1, idx2]))
                    pair_counter[pair] += 1

        if pair_counter:
            best_pair = max(pair_counter.items(), key=lambda x: x[1])[0]
            lower_idx, upper_idx = best_pair
            result_pairs.append([lower_idx, upper_idx])

            # Get the two good lines
            lower_line = dict(all_lines[lower_idx])
            upper_line = dict(all_lines[upper_idx])

            # Now correct bad line in x-values outside the sandwich
            corrected_line = []
            for x, y in bad_line:
                if x not in lower_line or x not in upper_line:
                    corrected_line.append((x, y))
                    continue

                y_low = lower_line[x]
                y_up = upper_line[x]
                low, up = sorted([y_low, y_up])

                if low <= y <= up:
                    corrected_line.append((x, y))  # No correction needed
                else:
                    corrected_y = (y_low + y_up) / 2
                    corrected_line.append((x, corrected_y))

            corrected_all_lines[bad_idx] = corrected_line
        else:
            result_pairs.append([-1, -1])  # No enclosing pair found

    return corrected_all_lines, result_pairs
    


def generate_code(points):
    """ inputs can be list of int or float or list of tuple (x,y) or ?!"""
    # Convert list of tuples to a byte string
    points_bytes = str(points).encode('utf-8')
    # Create a SHA-256 hash
    hash_object = hashlib.sha256(points_bytes)
    # Return a shortened version of the hex digest (e.g., first 8 characters)
    return hash_object.hexdigest()[:8]
    

# Thank God!! perfected axes
def plot_layers(
    layers,
    pixels=512,
    save_option=False,
    marked_layers=False,
    axis_off=True,
    close_option=False,
    fixed_name=None,
    inline=-999,
    initial_inline=0,
    crossline=0,
    fs=[14, 12, 11],
    line_color='red',
    x_start=0,
    y_start=0,
    twt_spacing=4,
    seismic_profile_shape=(950, 460),
    return_fig=False,
    svg_option=False
):

    # ---------------------------
    # Normalize inputs
    # ---------------------------
    if isinstance(layers[0], tuple):
        layers = [layers]

    if marked_layers:
        if isinstance(marked_layers[0], tuple):
            marked_layers = [marked_layers]

    dpi = 100

    if isinstance(pixels, int):
        width_px, height_px = pixels, pixels
    else:
        width_px, height_px = pixels

    fig, ax = plt.subplots(figsize=(width_px / dpi, height_px / dpi), dpi=dpi)

    # ============================================================
    # AXIS OFF  (image-like coordinates)
    # ============================================================
    if axis_off:
        ax.axis('off')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

        def process_layer(layer):
            pts = list(set(layer))
            return sorted(pts, key=lambda p: p[0])

        all_x, all_y = [], []

        # regular layers
        for layer in layers:
            cl = process_layer(layer)
            if cl:
                x = [p[0] for p in cl]
                y = [p[1] for p in cl]
                all_x.extend(x)
                all_y.extend(y)
                ax.plot(x, y, c='black', linewidth=1)

        # marked layers
        if marked_layers:
            for layer in marked_layers:
                cl = process_layer(layer)
                if cl:
                    x = [p[0] for p in cl]
                    y = [p[1] for p in cl]
                    all_x.extend(x)
                    all_y.extend(y)
                    ax.plot(x, y, c=line_color, linewidth=2)

        if all_x and all_y:
            ax.set_xlim(min(all_x), max(all_x))
            ax.set_ylim(min(all_y), max(all_y))

    # ============================================================
    # AXIS ON  (physical seismic coordinates)
    # ============================================================
    else:
        # seismic geometry
        n_traces, n_samples = seismic_profile_shape
        twt_height = n_samples * twt_spacing
        x_end = x_start + n_traces
        y_end = y_start + twt_height

        # debug
        print(*(f"{name}: {val}" for name, val in [("x_start", x_start), ("y_start", y_start), ("twt_spacing", twt_spacing), ("n_traces", n_traces), ("n_samples", n_samples), ("twt_height", n_samples * twt_spacing), ("x_end", x_start + n_traces), ("y_end", y_start + n_samples * twt_spacing)]))

        # --------------------------------------------------------
        # SCALE BOTH layers AND marked_layers IDENTICALLY
        # --------------------------------------------------------
        scaled_layers = rescale_lines_to_image(
            layers,
            dimensions=(x_end - x_start, y_end - y_start),
            offset=(x_start, y_start)
        )

        scaled_marked = None
        if marked_layers:
            scaled_marked = rescale_lines_to_image(
                marked_layers,
                dimensions=(x_end - x_start, y_end - y_start),
                offset=(x_start, y_start)
            )

        # ensure list-of-lists
        if isinstance(scaled_layers[0], tuple):
            scaled_layers = [scaled_layers]
        if scaled_marked and isinstance(scaled_marked[0], tuple):
            scaled_marked = [scaled_marked]

        # --------------------------------------------------------
        # Apply SAME coordinate transforms to both
        # --------------------------------------------------------
        
        #### transform and flipping was unncessary when y_start and y_end is in plotting and offset added in rescaling
        # def transform(lines):
            # out = []
            # max_y = max(y for ln in lines for _, y in ln)

            # for ln in lines:
                # tmp = [(x, max_y - y) for x, y in ln]   # removed  + x_start and  + y_start
                # out.append(tmp)

            # out = flip_layers(out)
            # return out 

        sample_point = scaled_layers[0][0]
        #x,samp, y_samp = sample_point
        print("sample_point:", sample_point)

        #scaled_layers = transform(scaled_layers)
        
        
        #if scaled_marked:
            #scaled_marked = transform(scaled_marked)
            
        sample_point = scaled_layers[0][0]
        #x,samp, y_samp = sample_point
        print("sample_point after flipping:", sample_point)            
            
        def process_layer(layer):
            pts = list(set(layer))
            return sorted(pts, key=lambda p: p[0])            

        all_x, all_y = [], []

        # plot regular layers
        for layer in scaled_layers:
            cl = process_layer(layer)
            if cl:
                x = [p[0] for p in cl]
                y = [p[1] for p in cl]
            all_x.extend(x)
            all_y.extend(y)
            ax.plot(x, y, c='black', linewidth=1)

        # plot marked layers (NOW CORRECTLY SCALED)
        if scaled_marked:
            for layer in scaled_marked:
                cl = process_layer(layer)
                if cl:
                    x = [p[0] for p in cl]
                    y = [p[1] for p in cl]
                all_x.extend(x)
                all_y.extend(y)
                ax.plot(x, y, c=line_color, linewidth=2)

        # labels & title
        ax.set_xlabel("CDP number", fontsize=fs[1])
        ax.set_ylabel("TWT [ms]", fontsize=fs[1])

        inline_plot = inline + initial_inline
        if crossline == 0:
            ax.set_title(f"{fixed_name} inline = {inline_plot}", fontsize=fs[0])
        else:
            ax.set_title(f"{fixed_name} crossline = {inline_plot}", fontsize=fs[0])

        # y ticks
        num_ticks = int(width_px / 85)
        twt_ticks = np.linspace(y_start, y_end, num_ticks)
        ax.set_yticks(twt_ticks)
        ax.set_yticklabels(list(reversed([f"{int(v)}" for v in twt_ticks])))

        if all_x and all_y:
            ax.set_xlim(min(all_x), max(all_x))
            ax.set_ylim(min(all_y), max(all_y))

    # ============================================================
    # SAVE / RETURN
    # ============================================================
    if save_option:
        if inline:
            tag = "inline" if crossline == 0 else "crossline"
            inline_plot = inline + initial_inline
            filename = f"layers_plot_{tag}_{inline_plot}.png"
        else:
            filename = f"layers_plot_{generate_code(all_x)}.png"

        fig.savefig(filename, bbox_inches='tight', pad_inches=0, dpi=dpi)
        base = find_base(filename)
        if svg_option:
            fig.savefig(base + ".svg", bbox_inches='tight', pad_inches=0, dpi=dpi)

    if return_fig:
        return fig

    if not close_option:
        plt.show()

    plt.close(fig)





def plot_layers_GUI(layers, marked_layers=False, pixels=512, line_color='red'):
    """Return a matplotlib Figure and Axes with layers and optional marked layers plotted."""
    # Normalize input
    if isinstance(layers[0], tuple):
        layers = [layers]
    if marked_layers and isinstance(marked_layers[0], tuple):
        marked_layers = [marked_layers]

    dpi = 100
    fig_size = pixels / dpi
    fig = Figure(figsize=(fig_size, fig_size), dpi=dpi)
    ax = fig.add_subplot(111)

    all_x_coords = []
    all_y_coords = []

    def process_layer(layer):
        unique_pts = list(set(layer))
        return sorted(unique_pts, key=lambda pt: pt[0])

    # Plot regular layers in black
    for layer in layers:
        cleaned = process_layer(layer)
        if cleaned:
            x = [pt[0] for pt in cleaned]
            y = [pt[1] for pt in cleaned]
            all_x_coords.extend(x)
            all_y_coords.extend(y)
            ax.plot(x, y, c='black', linewidth=1)

    # Plot marked layers in red
    if marked_layers:
        for layer in marked_layers:
            cleaned = process_layer(layer)
            if cleaned:
                x = [pt[0] for pt in cleaned]
                y = [pt[1] for pt in cleaned]
                all_x_coords.extend(x)
                all_y_coords.extend(y)
                ax.plot(x, y, c=line_color, linewidth=2)

    # Set axis limits
    if all_x_coords and all_y_coords:
        ax.set_xlim(min(all_x_coords), max(all_x_coords))
        ax.set_ylim(min(all_y_coords), max(all_y_coords))

    ax.axis('off')
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    return fig, ax






def generate_code_seismic(values):
    """Generate a hash code from a list of float or int values."""
    flat = ','.join(map(str, np.round(values, 2)))  # round to reduce hash noise
    return hashlib.md5(flat.encode()).hexdigest()[:8]



def plot_seismic(
    seismic_profile,
    pixels=(512, 512),  # can now be int or tuple (width_pixels, height_pixels)
    marked_points=False,  # list of (x, y) points in CDP/TWT space
    saving_directory=False,
    save_option=False,
    axis_off=False,
    close_option=False,
    seismic_relief=99
):

    # Validate input
    if not isinstance(seismic_profile, np.ndarray) or seismic_profile.ndim != 2:
        raise ValueError("seismic_profile must be a 2D numpy array")

    # Determine figure size
    dpi = 100
    if isinstance(pixels, int):
        figure_size_in_inches = (pixels / dpi, pixels / dpi)
    elif isinstance(pixels, (tuple, list)) and len(pixels) == 2:
        figure_size_in_inches = (pixels[0] / dpi, pixels[1] / dpi)
    else:
        raise ValueError("pixels must be an int or a tuple/list of length 2")

    fig = plt.figure(figsize=figure_size_in_inches, dpi=dpi)
    ax = fig.add_subplot(1, 1, 1)

    # Compute value range for color scale
    seismic_relief = np.clip(seismic_relief, 51, 100)
    vm = np.percentile(seismic_profile, seismic_relief)

    # Plot the seismic profile (transpose for CDP x TWT)
    im = ax.imshow(
        seismic_profile.T,
        cmap="RdBu",
        vmin=-vm,
        vmax=vm,
        aspect='auto'
    )

    # Overlay marked points, if any
    if marked_points:
        if isinstance(marked_points[0], (int, float)):
            marked_points = [marked_points]
        for x, y in marked_points:
            ax.plot(x, y, 'ro', markersize=4)

    # Axis configuration
    if axis_off:
        ax.axis('off')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    else:
        ax.set_xlabel('CDP number')
        ax.set_ylabel('TWT [ms]')
        ax.set_title('Seismic Profile')

    # Save if requested
    if save_option:
        all_values = list(seismic_profile.flatten())
        if marked_points:
            all_values += [coord for pt in marked_points for coord in pt]
        code = generate_code_seismic(all_values)
        plt.savefig(f"seismic_plot_{code}.png", bbox_inches='tight', pad_inches=0, dpi=dpi)
        #plt.savefig(f"seismic_plot_{code}.svg", format="svg", bbox_inches='tight', pad_inches=0, dpi=dpi)
        #plt.savefig(f"{saving_directory}/seismic_plot_{code}.eps", format="eps", bbox_inches='tight', pad_inches=0, dpi=dpi)            
    
    if saving_directory:
        all_values = list(seismic_profile.flatten())
        if marked_points:
            all_values += [coord for pt in marked_points for coord in pt]
        code = generate_code_seismic(all_values)
        plt.savefig(f"{saving_directory}/seismic_plot_{code}.png", bbox_inches='tight', pad_inches=0, dpi=dpi)
        #plt.savefig(f"{saving_directory}/seismic_plot_{code}.svg", format="svg", bbox_inches='tight', pad_inches=0, dpi=dpi)
        #plt.savefig(f"{saving_directory}/seismic_plot_{code}.eps", format="eps", bbox_inches='tight', pad_inches=0, dpi=dpi)        


    # Show unless suppressed
    if not close_option:
        plt.show()

    plt.close(fig)
    
    
    

def remove_subset_lines(lines):
    """
    Removes lines that are subsets of other equal or longer lines.

    Args:
        lines (list of list of tuple): Each line is a list of (x, y) coordinates,
                                       strictly increasing in positive x direction.

    Returns:
        list of list of tuple: Filtered lines with subset lines removed.
    """
    # Sort lines by length in descending order to check longer lines first
    sorted_lines = sorted(lines, key=lambda l: -len(l))
    filtered_lines = []
    
    for line in sorted_lines:
        # Check if this line is a subset of any already accepted line
        is_subset = False
        for accepted_line in filtered_lines:
            # Since both lines are strictly increasing in x, we can check subset by coordinates
            # We need to find if all points of 'line' exist in 'accepted_line'
            i = 0  # index for line
            j = 0  # index for accepted_line
            matched = 0
            
            while i < len(line) and j < len(accepted_line):
                if line[i] == accepted_line[j]:
                    matched += 1
                    i += 1
                    j += 1
                elif line[i][0] > accepted_line[j][0]:  # compare x-coordinates
                    j += 1
                else:
                    break  # line[i] not found in accepted_line
            
            if matched == len(line):
                is_subset = True
                break
        
        if not is_subset:
            filtered_lines.append(line)
    
    return filtered_lines


def plot_layers_scatter(
    layers,
    pixels=512,
    saving_directory=False,
    marked_layers=False,
    axis_off=True,
    close_option=False,
    line_color='red'
):
    # Normalize input to list of layers
    if isinstance(layers[0], tuple):
        layers = [layers]

    if marked_layers and isinstance(marked_layers[0], tuple):
        marked_layers = [marked_layers]

    # Plotting configuration
    dpi = 100
    figure_size_in_inches = pixels / dpi
    fig = plt.figure(figsize=(figure_size_in_inches, figure_size_in_inches), dpi=dpi)

    all_x_coords = []
    all_y_coords = []

    # Helper to clean and sort points by x
    def process_layer(layer):
        unique_pts = list(set(layer))
        return sorted(unique_pts, key=lambda pt: pt[0])

    # Plot regular layers in black (CHANGED TO SCATTER)
    for layer in layers:
        cleaned = process_layer(layer)
        if cleaned:
            x = [pt[0] for pt in cleaned]
            y = [pt[1] for pt in cleaned]
            all_x_coords.extend(x)
            all_y_coords.extend(y)
            plt.scatter(x, y, c='black', s=1)  # Changed this line only

    # Plot marked layers in red, on top (CHANGED TO SCATTER)
    if marked_layers:
        for layer in marked_layers:
            cleaned = process_layer(layer)
            if cleaned:
                x = [pt[0] for pt in cleaned]
                y = [pt[1] for pt in cleaned]
                all_x_coords.extend(x)
                all_y_coords.extend(y)
                plt.scatter(x, y, c=line_color, s=4)  # Changed this line only

    # Set limits based on all points
    if all_x_coords and all_y_coords:
        plt.xlim(min(all_x_coords), max(all_x_coords))
        plt.ylim(min(all_y_coords), max(all_y_coords))

    # Optional axis off
    if axis_off:
        plt.axis('off')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

    # Save option
    if saving_directory:
        code = generate_code(all_x_coords)
        plt.savefig(f"{saving_directory}/layers_plot_{code}.png", bbox_inches='tight', pad_inches=0, dpi=dpi)
        #plt.savefig(f"{saving_directory}/layers_plot_{code}.svg", format="svg", bbox_inches='tight', pad_inches=0, dpi=dpi)        
        #plt.savefig(f"{saving_directory}/layers_plot_{code}.eps", format="eps", bbox_inches='tight', pad_inches=0, dpi=dpi)                

    
    # Only show if close_option is False
    if not close_option:
        plt.show()

    plt.close(fig)





############################################## applying subsidence ########################################
def calculate_trend_line(lines):
    # Flatten the list of lines into a single list of points
    points = [coord for line in lines for coord in line]
    
    # Separate into X and Y coordinates
    x_coords = [p[0] for p in points]
    y_coords = [p[1] for p in points]
    
    # Handle edge case: single point or no points
    if len(points) < 2:
        raise ValueError("Insufficient points to calculate a trend line.")
    
    # Perform linear regression
    x = np.array(x_coords)
    y = np.array(y_coords)
    A = np.vstack([x, np.ones(len(x))]).T  # Create design matrix
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]  # Solve y = mx + c
    
    return m, c, x_coords, y_coords  # Return slope and intercept, and coordinates

def flip_positive_trend(lines, m):
    
    if m <= 0:
        # Do nothing if the slope is not positive
        return lines
    
    # Flatten all x-coordinates to find the range
    all_x = [x for line in lines for x, y in line]
    x_min, x_max = min(all_x), max(all_x)
    
    # Flip x-coordinates and shift them to start from 0
    flipped_lines = []
    for line in lines:
        flipped_line = [(-(x - x_max), y) for x, y in line]  # Flip x and shift (shifting is not allowed if we were to keep x values for real data)
        flipped_lines.append(flipped_line)
    
    return flipped_lines



def normalize_ssm(ssm):
    """
    Normalizes the y-values of SSM to be between 0 and 1.

    Args:
        ssm (list of tuples): List of (x, y) tuples.

    Returns:
        list of tuples: SSM with normalized y-values.
    """
    # Extract all y-values from the SSM
    xs, ys = zip(*ssm)
    y_min = min(ys)
    y_max = max(ys)

    # Prevent division by zero if all y-values are the same
    if y_max == y_min:
        return [(x, 0.0) for x in xs]

    # Normalize each y-value
    normalized_ssm = [(x, (y - y_min) / (y_max - y_min)) for x, y in ssm]
    return normalized_ssm



def calculate_top_and_bottom(lines):
    # Flatten the dataset into a list of (x, y) points
    all_points = [point for line in lines for point in line]

    # Group points by their x-values
    grouped_points = defaultdict(list)
    for x, y in all_points:
        grouped_points[x].append(y)

    # Calculate top and bottom lines
    top = [(x, max(y_values)) for x, y_values in sorted(grouped_points.items())]
    bottom = [(x, min(y_values)) for x, y_values in sorted(grouped_points.items())]

    return top, bottom

#merge ends to start
def merge_ends_to_starts(lines):

    merged_lines = []
    visited = [False] * len(lines)

    for i, line in enumerate(lines):
        if visited[i]:
            continue
        
        # Start a new merged line
        merged_line = line[:]
        visited[i] = True
        
        # Try to merge with other lines
        while True:
            merged = False
            for j, other_line in enumerate(lines):
                if visited[j]:
                    continue
                
                # Check for last point to first point connection
                if merged_line[-1] == other_line[0]:
                    merged_line.extend(other_line[1:])
                    visited[j] = True
                    merged = True
                    break
            
            if not merged:  # No more lines to merge
                break
        
        # Add the merged line to the result
        merged_lines.append(merged_line)

    return merged_lines




def compress_lines(lines, compression_percentage = 0.5):
    """
    Compresses a list of lines along the y-axis by a given percentage and returns the result as a list of lines.

    Args:
        lines (list of list of tuples): The input lines, where each line is a list of tuple coordinates.
        compression_percentage (float): The percentage to compress along the y-axis (e.g., 0.9 for 90%).

    Returns:
        list of list of tuples: The compressed lines in the same structure as the input.
    """
    # Create a Shapely MultiLineString from the list of lines
    multiline = MultiLineString(lines)
    
    # Calculate the compression factor along the y-axis
    compression_factor = compression_percentage
    
    # Apply scaling to compress along the y-axis
    compressed_multiline = scale(multiline, xfact=1, yfact=compression_factor, origin='center')
    
    # Convert the compressed MultiLineString back to a list of lines
    compressed_lines = [list(line.coords) for line in compressed_multiline.geoms]
    
    return compressed_lines



def heterogeneous_compress(lines, ssm):
    """
    Compresses a list of lines along the y-axis using heterogeneous compression ratios
    defined by a list of (x, ratio) tuples (SSM). Points not sharing x-values with SSM
    coordinates remain uncompressed.

    Args:
        lines (list of list of tuples): Input lines, where each line is a list of (x, y) tuples.
        ssm (list of tuples): List of (x, ratio) tuples defining compression ratios for x-values.

    Returns:
        list of list of tuples: The compressed lines in the same structure as the input.
    """
    # Convert SSM to a dictionary for fast lookups
    ssm_dict = dict(ssm)
    
    def compress_point(x, y):
        """Compress a point based on the compression ratio for its x-value."""
        if x in ssm_dict:
            ratio = ssm_dict[x]
            return (x, y * ratio)  # Apply compression only to y
        return (x, y)  # Keep the point unchanged if x is not in SSM

    compressed_lines = []
    for line in lines:
        # Compress each point in the line
        compressed_line = [compress_point(x, y) for x, y in line]
        compressed_lines.append(compressed_line)

    return compressed_lines


def increase_line_range(SSM, num_points=512, kind='linear', plot=False, y_offset=49):
    x, y = zip(*SSM)  # Unzip input data into x and y

    # Generate equally spaced x-values for interpolation
    x_new = np.linspace(min(x), max(x), num_points)
    
    # Create interpolation function
    interp_func = interp1d(x, y, kind=kind, fill_value='extrapolate')
    
    # Interpolate the y-values and apply the vertical offset
    y_new = interp_func(x_new) + y_offset  # Apply offset directly after interpolation
    
    # Combine x_new and y_new into a list of tuples
    interpolated_points = list(zip(x_new, y_new))
    
    # Plotting (if required)
    if plot:
        plt.figure(figsize=(8, 6))
        plt.plot(x, y, 'o', label='Original Points', markersize=5)
        plt.plot(x_new, y_new, '-', label=f'Interpolated ({kind})', linewidth=2)
        plt.legend()
        plt.title('Interpolation of Points with Offset')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)
        plt.show()

    return interpolated_points



def flip_and_compress(lines, ssm):
    """
    Compresses a list of lines along the y-axis using heterogeneous compression ratios
    defined by a list of (x, ratio) tuples (SSM), with the compression applied from bottom
    to top (flipping the original y-values).

    Args:
        lines (list of list of tuples): Input lines, where each line is a list of (x, y) tuples.
        ssm (list of tuples): List of (x, ratio) tuples defining compression ratios for x-values.

    Returns:
        list of list of tuples: The compressed lines in the same structure as the input.
    """
    # Convert SSM to a dictionary for fast lookups
    ssm_dict = dict(ssm)
    
    def compress_point(x, y):
        """Compress a point based on the compression ratio for its x-value."""
        if x in ssm_dict:
            ratio = ssm_dict[x]
            return (x, -y * ratio)  # Apply compression to the negated y (flip the compression direction)
        return (x, -y)  # Keep the point unchanged if x is not in SSM

    # Flip all lines (negate y-values) to apply compression from bottom to top
    flipped_lines = [[(x, -y) for x, y in line] for line in lines]
    
    # Apply the compression to the flipped lines
    compressed_lines = []
    for line in flipped_lines:
        # Compress each point in the line after flipping
        compressed_line = [compress_point(x, y) for x, y in line]
        compressed_lines.append(compressed_line)

    # Flip the lines back to their original y-axis orientation
    final_lines = [[(x, -y) for x, y in line] for line in compressed_lines]

    return final_lines


def remove_duplicate_shared_points(lines):  # used by spatial3

    point_map = defaultdict(list)

    # Map each point to the indices of lines it appears in
    for i, line in enumerate(lines):
        for point in line:
            point_map[point].append(i)

    # Identify shared points
    shared_points = {p for p, l in point_map.items() if len(l) > 1}

    # Retain only one instance of each shared point
    used_shared_points = set()
    new_lines = [[] for _ in lines]

    for i, line in enumerate(lines):
        for point in line:
            if point in shared_points:
                if point not in used_shared_points:
                    new_lines[i].append(point)
                    used_shared_points.add(point)
            else:
                new_lines[i].append(point)

    return new_lines

def remove_none_and_empty_lines(lines):
    # Filter out any line that is either None or an empty list
    return [line for line in lines if line is not None and line != []]



def lines_to_TX(lines, T, X, integer = True):
    """
    Transforms a list of Bresenham lines into a (T, X) array by averaging, scaling, and interpolating.
    """
    # Step 1: Compute x_groups
    x_groups = {}
    for line in lines:
        for x, y in line:
            if x not in x_groups:
                x_groups[x] = []
            x_groups[x].append(y)
    
    # Average y values for each x
    x_groups = {x: np.mean(ys) for x, ys in x_groups.items()}
    
    # Sort x_groups by x values
    sorted_x = sorted(x_groups.keys())
    sorted_y = [x_groups[x] for x in sorted_x]
    
    # Step 2: Scale x-axis to X points
    if len(sorted_x) > X:
        # Downsample
        x_new = np.linspace(min(sorted_x), max(sorted_x), X)
        f = interp1d(sorted_x, sorted_y, kind='linear', assume_sorted=True)
        y_new = f(x_new)
    elif len(sorted_x) < X:
        # Upsample
        x_new = np.linspace(min(sorted_x), max(sorted_x), X)
        f = interp1d(sorted_x, sorted_y, kind='linear', assume_sorted=True)
        y_new = f(x_new)
    else:
        # No scaling needed
        x_new = sorted_x
        y_new = sorted_y
    
    # Step 3: Scale number of lines to T lines
    # Create a 2D array of y values for each line
    y_lines = []
    for line in lines:
        line_x = [x for x, y in line]
        line_y = [y for x, y in line]
        
        # Remove duplicate x values to avoid interpolation errors
        unique_x, unique_indices = np.unique(line_x, return_index=True)
        unique_y = np.array(line_y)[unique_indices]
        
        # Interpolate only if there are at least 2 unique x values
        if len(unique_x) >= 2:
            f = interp1d(unique_x, unique_y, kind='linear', fill_value="extrapolate")
            y_lines.append(f(x_new))
        else:
            # If there's only one unique x value, repeat the y value for all x_new
            y_lines.append(np.full_like(x_new, unique_y[0]))
    
    # Downsample or upsample the lines to T lines
    if len(y_lines) > T:
        # Downsample
        indices = np.linspace(0, len(y_lines) - 1, T).astype(int)
        y_lines = [y_lines[i] for i in indices]
    elif len(y_lines) < T:
        # Upsample
        x_old = np.arange(len(y_lines))
        x_new = np.linspace(0, len(y_lines) - 1, T)
        y_lines = [np.interp(x_new, x_old, [line[i] for line in y_lines]) for i in range(X)]
        y_lines = np.array(y_lines).T.tolist()
    
    # Convert to a (T, X) array
    output_array = np.array(y_lines)
    
    if integer == True:
        # Making all values to become integers
        output_array = output_array.astype(int)
    
    return output_array


def TX_to_lines(output_array):
    """
    Convert a (T, X) numpy array into a list of lines,
    where each line is a list of (x, y) tuples.
    """
    T, X = output_array.shape
    lines = []
    for t in range(T):
        line = [(x, output_array[t, x]) for x in range(X)]
        lines.append(line)
    return lines



def remove_steep_sections(lines, threshold=10):
    """
    Remove lines that contain at least one steep segment.
    """
    lines_to_be_removed = []

    for idx, line in enumerate(lines):
        for i in range(len(line) - 1):
            (x1, y1), (x2, y2) = line[i], line[i + 1]
            dx = abs(x2 - x1)
            dy = abs(y2 - y1)
            if dx == 1 and dy >= threshold:
                lines_to_be_removed.append(idx)
                break  # No need to check further points in this line

    # Keep only lines not in lines_to_be_removed
#     mod_lines = [line for idx, line in enumerate(lines) if idx not in lines_to_be_removed]
    return lines_to_be_removed

def correct_steep_frames(lines, threshold=10):
    """
    Modify steep sections in lines based on dx/dy and boundary proximity.
    """
    # Flatten all x values to get global min/max
    all_x = [x for line in lines for x, y in line]
    global_x_min = min(all_x)
    global_x_max = max(all_x)

    mod_lines = []

    for line in lines:
        new_line = line.copy()
        for i in range(len(line) - 1):
            x1, y1 = new_line[i]
            x2, y2 = new_line[i + 1]
            dx = abs(x2 - x1)
            dy = abs(y2 - y1)
            X_min = x2 - global_x_min
            X_max = global_x_max - x2
            if dx == 1 and dy >= threshold and (X_min <= threshold or X_max <= threshold):
                # Replace y2 with y1
                new_line[i + 1] = (x2, y1)
        mod_lines.append(new_line)

    return mod_lines


def remove_steep_frames(lines, threshold=3):
    """
    Modify steep sections in lines based on dx/dy and boundary proximity.
    """
    # Flatten all x values to get global min/max
    all_x = [x for line in lines for x, y in line]
    global_x_min = min(all_x)
    global_x_max = max(all_x)

    mod_lines = []
    lines_to_remove = []

    for idx, line in enumerate(lines):
        new_line = line.copy()
        for i in range(len(line) - 1):
            x1, y1 = new_line[i]
            x2, y2 = new_line[i + 1]
            dx = abs(x2 - x1)
            dy = abs(y2 - y1)
            X_min = x2 - global_x_min
            X_max = global_x_max - x2
            if dx == 1 and dy >= threshold and (X_min <= threshold or X_max <= threshold):
                lines_to_remove.append(idx)

    return [line for idx,line in enumerate(lines) if idx not in lines_to_remove]    #mod_lines


def find_line_tops(lines, threshold=10):
    # Step 1: Flatten all points with their line index
    all_points = []
    for idx, line in enumerate(lines):
        for x, y in line:
            all_points.append((x, y, idx))

    # Step 2: Group points by similar x values
    x_groups = defaultdict(list)
    for x, y, idx in all_points:
        found_group = False
        for gx in x_groups:
            if abs(gx - x) <= threshold:
                x_groups[gx].append((x, y, idx))
                found_group = True
                break
        if not found_group:
            x_groups[x].append((x, y, idx))

    # Step 3: Find which line owns min_y and max_y in each x_group
    top_counter = defaultdict(int)
    bottom_counter = defaultdict(int)

    for group_points in x_groups.values():
        # Find top (max y) and bottom (min y) in this group
        max_y_point = max(group_points, key=lambda p: p[1])  # p[1] is y
        min_y_point = min(group_points, key=lambda p: p[1])

        top_counter[max_y_point[2]] += 1  # Line index of max y
        bottom_counter[min_y_point[2]] += 1  # Line index of min y

    # Step 4: Find global top and bottom line
    global_top_line = max(top_counter.items(), key=lambda x: x[1])[0]
    global_bottom_line = max(bottom_counter.items(), key=lambda x: x[1])[0]

    return global_top_line, global_bottom_line





###################################### Thank God!! spstR function replacement for spatial3 and R functions ###############################
def N_slope(lines, N=5):
    """Return dictionaries of first-N and last-N slopes for each line index.
    Distinguishes horizontal/vertical cases using ε and ±∞.
    """
    N_first_slope, N_last_slope = {}, {}
    eps = 1e-6  # small epsilon to encode horizontal direction

    mother_N = N
    for i, line in enumerate(lines):
        if len(line) < 2:
            N_first_slope[i] = "short"
            N_last_slope[i] = "short"
            continue
            
        if len(line) < N:
            N = len(line)
        else:
            N = mother_N

        # --- first slope ---
        temp_N = N
        x1, y1 = line[0]
        idx2 = min(N-1, len(line)-1)
        x2, y2 = line[idx2]
        while (x2 == x1 or y2 == y1) and temp_N <= len(line):  
            temp_N += 1
            x1, y1 = line[0]
            idx2 = min(temp_N-1, len(line)-1)
            x2, y2 = line[idx2]

        if x2 == x1:  # vertical
            N_first_slope[i] = math.inf if y2 > y1 else -math.inf
        elif y2 == y1:  # horizontal
            N_first_slope[i] = eps if x2 > x1 else -eps
        else:
            N_first_slope[i] = (y2 - y1) / (x2 - x1)

        # --- last slope ---
        temp_N = N
        idx1 = -min(N, len(line))
        x1, y1 = line[idx1]
        x2, y2 = line[-1]
        while (x2 == x1 or y2 == y1) and temp_N <= len(line):  
            temp_N += 1                
            idx1 = -min(temp_N, len(line))
            x1, y1 = line[idx1]
            x2, y2 = line[-1]

        if x2 == x1:  # vertical
            N_last_slope[i] = math.inf if y2 > y1 else -math.inf
        elif y2 == y1:  # horizontal
            N_last_slope[i] = eps if x2 > x1 else -eps
        else:
            N_last_slope[i] = (y2 - y1) / (x2 - x1)



    print("N_slopes <=> finished")
    return N_first_slope, N_last_slope



def bresenham(x1, y1, x2, y2):
    """Generate integer coordinates on a line from (x1,y1) to (x2,y2)."""
    points = []
    dx, dy = abs(x2 - x1), abs(y2 - y1)
    sx, sy = (1 if x1 < x2 else -1), (1 if y1 < y2 else -1)
    err = dx - dy
    while True:
        points.append((x1, y1))
        if x1 == x2 and y1 == y2:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x1 += sx
        if e2 < dx:
            err += dx
            y1 += sy
    return points



def endpoint(x0, y0, slope, x_min, x_max, y_min, y_max, end):
    """
    Return the endpoint (integer) of a line starting at (x0, y0) with given slope
    that intersects the bounding box sides. Only checks eligible edges depending 
    on whether this is the first (end=0) or last (end=-1) point of the line.
    
    Eligible edges:
      - end == 0 (first): left, top, bottom
      - end == -1 (last): right, top, bottom
    Picks the intersection with the shortest Euclidean distance.
    """
    candidates = []
    #     print("slope, end",slope, end)
    if slope == "short":
        print("this should not be printed in spstR function")
        return (x0, y0)  # fallback

    # --- Handle infinite slopes (vertical lines) ---
    if np.isinf(slope):
        # For vertical lines, x doesn't change; extend straight up or down
        x_target = x0
        if slope > 0:  # upward vertical
            y_target = y_max if end == -1 else y_min
        else:  # downward vertical
            y_target = y_min if end == -1 else y_max
        # Clip to bounds
        y_target = min(max(y_target, y_min), y_max)
        x_target = min(max(x_target, x_min), x_max)
        return (int(round(x_target)), int(round(y_target)))

    # --- Regular finite slope logic ---
    if slope > 0:
        if end == 0:  # first, positive slope -> bottom or left
            eligible = ["bottom", "left"]
        else:  # last, positive slope -> top or right
            eligible = ["top", "right"]
    else:  # slope < 0
        if end == 0:  # first, negative slope -> top or left
            eligible = ["top", "left"]
        else:  # last, negative slope -> bottom or right
            eligible = ["bottom", "right"]

    # left edge
    if "left" in eligible:
        x_target = x_min
        y_target = round(y0 + slope * (x_target - x0))
        if y_min <= y_target <= y_max:
            candidates.append((x_target, y_target))

    # right edge
    if "right" in eligible:
        x_target = x_max
        y_target = round(y0 + slope * (x_target - x0))
        if y_min <= y_target <= y_max:
            candidates.append((x_target, y_target))

    # top edge
    if "top" in eligible:
        y_target = y_max
        x_target = round(x0 + (y_target - y0) / slope)
        if x_min <= x_target <= x_max:
            candidates.append((x_target, y_target))

    # bottom edge
    if "bottom" in eligible:
        y_target = y_min
        x_target = round(x0 + (y_target - y0) / slope)
        if x_min <= x_target <= x_max:
            candidates.append((x_target, y_target))

    if not candidates:
        return (x0, y0)  # fallback

    # pick closest candidate
    return min(candidates, key=lambda pt: math.hypot(pt[0] - x0, pt[1] - y0))



def tpath(lines, N=5, debug=False, debug_dir="debuging_gif"):
    """Compute traversing paths for both ends of each line.
       If debug=True, saves a plot with lines + traversing paths.
    """
    N_first_slope, N_last_slope = N_slope(lines, N)

    traversing_path = {}

    # flatten all points and compute global bounds
    all_points = [pt for line in lines for pt in line]
    if not all_points:
        return traversing_path
    x_min, y_min = map(min, zip(*all_points))
    x_max, y_max = map(max, zip(*all_points))

    if debug and not os.path.exists(debug_dir):
        os.makedirs(debug_dir)

    for i, line in enumerate(lines):
        if not line:
            continue

        # first
        slope = N_first_slope.get(i)
        x0, y0 = line[0]
        x1, y1 = endpoint(x0, y0, slope, x_min, x_max, y_min, y_max, end = 0)
        path_first = bresenham(x0, y0, x1, y1)
        traversing_path[(i, 0)] = path_first        

        # last
        slope = N_last_slope.get(i)
        x0, y0 = line[-1]
        x1, y1 = endpoint(x0, y0, slope, x_min, x_max, y_min, y_max, end = -1)
        path_last = bresenham(x0, y0, x1, y1)
        traversing_path[(i, -1)] = path_last

        # ---------- DEBUG PLOTTING ----------
        if debug:
            # mark both traversing paths for this line
            marked_layers = [path_first, path_last, line]
            plot_layers(
                lines,
                pixels=256,
                saving_directory=debug_dir,
                marked_layers=marked_layers,
                axis_off=True,
                close_option=True
            )

    print("traversing_paths <=> finished")
    return traversing_path



def square_offset(coord, r=1):
    """Return all integer coordinates in a square neighborhood of radius r."""
    x, y = coord
    return [(x+dx, y+dy) for dx in range(-r, r+1) for dy in range(-r, r+1)]




def spstR(lines, N=5, r=1, debug=False, DEBUG_it = 10, debug_tpath = False):
    """Merge lines by traversing from ends and attaching neighbouring line segments.
    Adds debug plotting every 100 iterations and saves progression as GIF.
    """

    # work on a copy
    lines = [list(line) for line in lines]

    lines = sorted(lines, key=lambda line: len(line), reverse=True)
    lines = [line for line in lines if len(line) >= 2]   # still we need 2 because smooth function will have a problem

    coord_set = {pt for line in lines for pt in line}
    coord_map = {pt: i for i, line in enumerate(lines) for pt in line}
    removables = set()

    traversing_path = tpath(lines, N, debug=debug_tpath)
#     return print("test finished")

    # --- DEBUG directory ---
    if debug:
        debug_dir = "debuging_gif"
        os.makedirs(debug_dir, exist_ok=True)
        # clear previous frames
        for f in glob.glob(os.path.join(debug_dir, "*.png")):
            os.remove(f)

    step_counter = 0  # count traversing iterations

    for i in range(len(lines)):
        if i in removables:
            continue
        if not lines[i]:
            continue

        for end in (0, -1):
            path = traversing_path.get((i, end), [])
            if not path: # if lines[i] did not have any path for that end, meaning its end is already a global min max
                continue

            traversing_coord = path[0] if path else None
            while path and traversing_coord  != path[-1]:
                # to be certain the traversing is moving outside line i
                while path and traversing_coord in lines[i] and traversing_coord != path[-1]: # path[-1] marks the global limit
                    if end == 0:
                        traversing_coord = path[0] if path else None
                        path = path[1:]
                    else: 
                        traversing_coord = path[0] if path else None
                        path = path[1:]
                        
                if not path:  # Exit if path becomes empty
                    break                       
                # update traversing_path; useful when no lines[j] is found
                traversing_path[(i, end)] = path
                
                #### start traversing outside lines i
                # adding out going traversing_coord
                if end == 0:
                    lines[i].insert(0, traversing_coord)
                else:
                    lines[i].append(traversing_coord)

                coord_map[traversing_coord] = i
                coord_set.add(traversing_coord)                    
                ###
                
                for sq in square_offset(traversing_coord, r=r):  # square coords check
                    if sq in coord_set and sq not in lines[i]:  # condition for finding lines[j]
                        
                        j = coord_map[sq]
                        
                        if j in removables:  # skip if j was not belonged to any present line in lines
                            continue
                        
                        # square found a new line => extend lines[i] by sq-lines[j] segment
                        # Find position of sq in lines[j] and return sublist from start to sq (inclusive)
                        position = lines[j].index(sq)
                        if end == 0:  # prepend
                            segment = lines[j][:position + 1]
                            split_segment = lines[j][position + 1:]  # remaining part after sq
                            lines[i] = segment + lines[i]
                        else:  # append
                            segment = lines[j][position:]
                            split_segment = lines[j][:position]  # remaining part before sq
                            lines[i] = lines[i] + segment

                        # Update lines[j] with the split_segment (or remove if empty)
                        if not split_segment:
#                             lines[j] = split_segment
#                         else: 
                            # If no remaining segment, mark for removal ; this shows perfect tsp connection!
                            removables.add(j)

                        # remove possible duplicates caused by traversing_coord == sq since we add t.._c.. early on
                        lines[i] = list(dict.fromkeys(lines[i]))
                        
                        # update coord_map for all coords of lines[j]
                        for pt in lines[i]:
                            coord_map[pt] = i

                        # switch traversing path to j
                        path = traversing_path.get((j, end), [])  # new path
#                         print('j, end',j, end)    
                        # update traversing map for all points on line i
                        traversing_path[(i, end)] = traversing_path[(j, end)]
                        break # break from square check
                
                # --- DEBUG periodic plot ---
                if debug and step_counter % DEBUG_it == 0:
                    timestamp = int(time.time() * 1000)
                    plot_layers(
                        lines,
                        pixels=256,
                        saving_directory=debug_dir,
                        marked_layers=[],
                        axis_off=True,
                        close_option=True,
                    )

    merged_lines = [line for k, line in enumerate(lines) if k not in removables]
    print("merging <=> finished")

    # --- Build GIF from saved PNGs ---
    if debug:
        frames = []
        files = sorted(glob.glob(os.path.join(debug_dir, "*.png")), key=os.path.getmtime)
        for fname in files:
            frames.append(Image.open(fname))
        if frames:
            gif_path = os.path.join(debug_dir, "debug_progress.gif")
            frames[0].save(
                gif_path,
                save_all=True,
                append_images=frames[1:],
                duration=200,
                loop=0,
            )
            print(f"Debug GIF saved to {gif_path}")
    
    
#     # Force close all matplotlib figures and clear memory
#     plt.close('all')
#     plt.clf()
#     plt.cla()
#     gc.collect()
    return merged_lines
    
    
##################### below are to improve smsl   #################################

# Thank God!! but transform outperformed by squeeze_lines_to_safe_min()


def transfer_bottom_to_top(original_lines, new_lines, verbose=False):
    """
    Transfer Y values in `new_lines` so the bottoms match `original_lines` while
    points at the top keep their Y (linearly interpolated per x-column).

    Assumptions (based on prior conversation):
      - `original_lines` and `new_lines` are lists of lines.
      - Each line is a list of (x, y) tuples.
      - Structure is aligned: same number of lines, and every line has the same
        x positions in the same order (only y differs). The function checks
        the X sets match; if not, it raises ValueError.
      - Grouping by X is strict equality (no tolerance).

    Args:
      original_lines: list[list[tuple(x, y)]] -- reference lines
      new_lines:      list[list[tuple(x, y)]] -- lines to be shifted
      verbose:        bool -- print small diagnostics (default False)

    Returns:
      transferred_lines: list[list[(x, new_y)]] with same structure as new_lines.
    """
    # collect min from original (orig_bottom_min_ys), min and max from new
    orig_min = {}
    new_min = {}
    new_max = {}

    # helpers to initialize
    def init_if_absent(d, k):
        if k not in d:
            d[k] = None

    # gather from original_lines
    for line in original_lines:
        for x, y in line:
            init_if_absent(orig_min, x)
            if y is None:
                continue
            if orig_min[x] is None or y < orig_min[x]:
                orig_min[x] = y

    # gather from new_lines
    for line in new_lines:
        for x, y in line:
            init_if_absent(new_min, x)
            init_if_absent(new_max, x)
            if y is None:
                continue
            if new_min[x] is None or y < new_min[x]:
                new_min[x] = y
            if new_max[x] is None or y > new_max[x]:
                new_max[x] = y

    # Check keys match (same x-groups)
    orig_xs = set(orig_min.keys())
    new_xs = set(new_min.keys())
    if orig_xs != new_xs:
        # be explicit about mismatch
        missing_in_new = orig_xs - new_xs
        missing_in_orig = new_xs - orig_xs
        raise ValueError(
            "X-groups mismatch between original_lines and new_lines.\n"
            f"Missing in new_lines: {sorted(missing_in_new)}\n"
            f"Missing in original_lines: {sorted(missing_in_orig)}"
        )

    # order x-values consistently (sorted ascending)
    x_values = sorted(orig_min.keys())

    # compute bottom_delta_ys: orig_bottom - new_bottom
    bottom_delta = {x: (orig_min[x] - new_min[x]) for x in x_values}

    if verbose:
        print("x_values:", x_values)
        print("orig_bottom_min_ys:", [orig_min[x] for x in x_values])
        print("new_bottom_min_ys: ", [new_min[x] for x in x_values])
        print("new_top_max_ys:   ", [new_max[x] for x in x_values])
        print("bottom_delta_ys:  ", [bottom_delta[x] for x in x_values])

    # perform transfer on each point of new_lines
    transferred_lines = []
    for line in new_lines:
        new_line = []
        for x, y in line:
            # if y is None, preserve as None
            if y is None:
                new_line.append((x, None))
                continue

            b = new_min[x]   # bottom in new
            t = new_max[x]   # top in new
            d = bottom_delta[x]  # desired delta to push bottom to orig bottom

            # if top==bottom, give full delta to all points (avoid division by zero)
            if t == b:
                transfer = d
            else:
                # linear from bottom (y=b => factor=1) to top (y=t => factor=0)
                factor = (t - y) / (t - b)
                # clamp factor between 0 and 1 to be safe
                if factor < 0:
                    factor = 0.0
                elif factor > 1:
                    factor = 1.0
                transfer = d * factor

            new_y = y + transfer
            new_line.append((x, new_y))
        transferred_lines.append(new_line)

    return transferred_lines



# Thank God!! this exension is good.
def extend_critical_lines(lines):
    """
    Extend lines that start or end at global y_min or y_max toward global x_min/x_max,
    then run bresenham_line() on them.

    Args:
        lines: List of lines, where each line is a list of (x, y) tuples.

    Returns:
        Modified list of lines with extended critical lines processed by bresenham_line().
    """
    if not lines:
        return []

    # --- Find global min/max x, y ---
    all_points = [pt for line in lines for pt in line if pt is not None]
    x_vals, y_vals = zip(*all_points)
    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)

    modified_lines = []
    for line in lines:
        if not line:
            modified_lines.append(line)
            continue

        new_line = line[:]
        first, last = line[0], line[-1]

        # --- Check critical starts/ends ---
        if first[1] in (y_min, y_max):
            critical_end = (x_min, first[1])  # extend left
            new_line = [critical_end] + new_line

        if last[1] in (y_min, y_max):
            critical_end = (x_max, last[1])  # extend right
            new_line = new_line + [critical_end]

        # --- If critical, apply Bresenham ---
        if new_line != line:
            new_line = bresenham_line(new_line)

        modified_lines.append(new_line)

    return modified_lines




# Thank God!! this squeezing, works somehow good  # Thank God!! by correcting global sorting to x-group sorting, now smooth function does not need this

def squeeze_lines_to_safe_min(lines):
    """
    Squeeze any x-group whose minimum y is below a computed safe_min upward so no point
    remains below safe_min, while keeping the top (y_max) of each affected group fixed
    when possible.

    Algorithm (as requested):
    - Find global x_min, x_max and global y_min, y_max (for information).
    - Group points strictly by x value (no tolerance).
    - For each x-group compute sorted_y, local y_min, y_max and local_delta = y_max - y_min.
    - Compute left_delta (x == x_min) and right_delta (x == x_max).
    - border_min = min(left_delta, right_delta); take the x-group (left or right)
      that attains border_min and set safe_min = that group's y_min.
    - For every x-group whose local y_min < safe_min:
        * If local_delta > 0:
            - For a point with y, compute t = (y - y_min) / local_delta  (0..1)
            - delta_y = (1 - t) * (safe_min - y_min)
            - new_y = y + delta_y
            (This makes bottom move up by full amount, top remain unchanged.)
        * If local_delta == 0 (all points in group have same y):
            - move all points to safe_min (degenerate case — necessary to ensure no point below safe_min).
    - Return a new list of lines with modified (x, y) tuples, preserving original line/point order.

    Args:
        lines: list of lines, where each line is a list of (x, y) tuples. x and y may be ints or floats.

    Returns:
        modified_lines: same structure as `lines` but with adjusted y values so no point is below safe_min.
    """
    # Flatten with references to be able to write back to lines structure
    entries = []  # each: (x, y, line_idx, pt_idx)
    for li, line in enumerate(lines):
        for pi, (x, y) in enumerate(line):
            entries.append({"x": x, "y": y, "li": li, "pi": pi})

    if not entries:
        return []  # nothing to do

    # Compute global x bounds
    xs = [e["x"] for e in entries]
    x_min = min(xs)
    x_max = max(xs)

    # Group strictly by x value
    groups = {}
    for e in entries:
        groups.setdefault(e["x"], []).append(e)

    # For each group compute sorted y, y_min, y_max, delta
    group_stats = {}
    for x_val, g in groups.items():
        ys = [e["y"] for e in g]
        y_min = min(ys)
        y_max = max(ys)
        delta = y_max - y_min
        # keep also sorted list (ascending)
        sorted_by_y = sorted(g, key=lambda ee: ee["y"])
        group_stats[x_val] = {
            "members": g,
            "sorted": sorted_by_y,
            "y_min": y_min,
            "y_max": y_max,
            "delta": delta,
        }

    # Ensure border x groups exist (they must, since entries not empty)
    left_stats = group_stats[x_min]
    right_stats = group_stats[x_max]
    left_delta = left_stats["delta"]
    right_delta = right_stats["delta"]

    # Choose border with minimum delta (if equal, choose left by <=)
    if left_delta <= right_delta:
        border_x = x_min
    else:
        border_x = x_max
    safe_min = group_stats[border_x]["y_min"]

    # Prepare a structure to collect new y values
    # copy original structure
    modified_lines = [list(line) for line in lines]  # shallow copy of tuples (we'll replace tuples)

    # For each group, if its y_min < safe_min, shift points
    for x_val, stats in group_stats.items():
        local_y_min = stats["y_min"]
        local_y_max = stats["y_max"]
        local_delta = stats["delta"]
        members = stats["members"]

        if local_y_min >= safe_min:
            # no need to change this group
            continue

        # amount bottom must be raised
        total_raise = safe_min - local_y_min  # positive

        if local_delta == 0:
            # degenerate: all points share same y. move all to safe_min.
            for e in members:
                li, pi = e["li"], e["pi"]
                x = e["x"]
                modified_lines[li][pi] = (x, safe_min)
        else:
            # normal case: interpolate delta per point so top (y_max) shifts by 0 and bottom by total_raise
            # using t = (y - y_min) / local_delta, delta_y = (1 - t) * total_raise
            for e in members:
                li, pi = e["li"], e["pi"]
                x = e["x"]
                y = e["y"]
                t = (y - local_y_min) / local_delta
                # clamp t for numeric safety
                if t < 0:
                    t = 0.0
                elif t > 1:
                    t = 1.0
                delta_y = (1.0 - t) * total_raise
                new_y = y + delta_y
                modified_lines[li][pi] = (x, new_y)

    return modified_lines




# [TG] but seems useless when smsl is so detailed (removes no line)
# Updated function with debugging statements that report counts and intersections.

def remove_lines_for_missing_border_points(original_lines, new_lines, verbose=True):
    """
    Same logic as previous: rescale new_lines to original bbox; remove entire lines from new_scaled_lines
    when left/right border integer-points are present in original but missing in new_scaled_lines.
    
    Adds debugging output when verbose=True:
      - counts of integer-grid points in original and scaled-new
      - count and sample of intersection
      - list of missing border points (that are in original but not in new_scaled)
      - indices of lines marked for removal (removables)
    
    Returns:
      (result_lines, stats_dict)
    """
    def bbox(lines):
        xs, ys = [], []
        for ln in lines:
            for p in ln:
                try:
                    x, y = float(p[0]), float(p[1])
                except Exception:
                    continue
                xs.append(x); ys.append(y)
        if not xs or not ys:
            return None
        return min(xs), max(xs), min(ys), max(ys)
    
    orig_lines = deepcopy(original_lines)
    new_lines_copy = deepcopy(new_lines)
    
    orig_bbox = bbox(orig_lines)
    new_bbox = bbox(new_lines_copy)
    if orig_bbox is None:
        raise ValueError("original_lines contains no valid points")
    oxmin, oxmax, oymin, oymax = orig_bbox
    
    # rescale new_lines to original bbox (same scaling logic)
    if new_bbox is None:
        new_scaled_lines = deepcopy(new_lines_copy)
    else:
        nxmin, nxmax, nymin, nymax = new_bbox
        new_scaled_lines = []
        def scale_val(v, vmin_from, vmax_from, vmin_to, vmax_to):
            if vmax_from == vmin_from:
                return vmin_to
            ratio = (v - vmin_from) / (vmax_from - vmin_from)
            return vmin_to + ratio * (vmax_to - vmin_to)
        for ln in new_lines_copy:
            scaled_ln = []
            for p in ln:
                try:
                    x, y = float(p[0]), float(p[1])
                except Exception:
                    continue
                sx = scale_val(x, nxmin, nxmax, oxmin, oxmax)
                sy = scale_val(y, nymin, nymax, oymin, oymax)
                scaled_ln.append((sx, sy))
            new_scaled_lines.append(scaled_ln)
    
    # build integer-grid border points for left and right sides
    x_left = int(math.floor(oxmin))
    x_right = int(math.ceil(oxmax))
    y_start = int(math.floor(oymin))
    y_end = int(math.ceil(oymax))
    left_border_points = [(x_left, y) for y in range(y_start, y_end + 1)]
    right_border_points = [(x_right, y) for y in range(y_start, y_end + 1)]
    border_points = left_border_points + right_border_points  # validated
    
    def intcoord(pt):
        return (int(round(pt[0])), int(round(pt[1])))
    
    # build sets for original and new_scaled integer coords
    orig_set = set()
    for ln in orig_lines:
        for p in ln:
            try:
                orig_set.add(intcoord(p))
            except Exception:
                pass
    
    new_set = set()
    # mapping from integer coord to list of line indices
    int_to_lines = {}
    # mapping from integer x to set of line indices (lines that have any point with that x after rounding)
    x_to_line_indices = {}
    for i, ln in enumerate(new_scaled_lines):
        for p in ln:
            try:
                ic = intcoord(p)
            except Exception:
                continue
            new_set.add(ic)
            int_to_lines.setdefault(ic, []).append(i)
            x_to_line_indices.setdefault(ic[0], set()).add(i)
    
    # DEBUG: counts and intersection
    intersection = orig_set & new_set
    stats = {
        "orig_count": len(orig_set),
        "new_scaled_count": len(new_set),
        "intersection_count": len(intersection),
        "intersection_sample": list(sorted(intersection))[:10],
        "orig_sample": list(sorted(orig_set))[:10],
        "new_sample": list(sorted(new_set))[:10],
    }
    if verbose:
        print("DEBUG: integer-grid counts -> original:", stats["orig_count"], 
              " new_scaled:", stats["new_scaled_count"], 
              " intersection:", stats["intersection_count"])
        print("DEBUG: intersection sample (up to 10):", stats["intersection_sample"])
        print("DEBUG: first 10 original points (int-rounded):", stats["orig_sample"])
        print("DEBUG: first 10 new_scaled points (int-rounded):", stats["new_sample"])
    
    # find missing border points
    missing_border_points = [bp for bp in border_points if (bp in orig_set and bp not in new_set)]
    if verbose:
        print("DEBUG: total border points:", len(border_points))
        print("DEBUG: missing border points count:", len(missing_border_points))
        if len(missing_border_points) > 0:
            print("DEBUG: missing border points (sample up to 50):", missing_border_points[:50])
    
    removables = set()
    # iterate the left/right border points and mark entire lines for removal when missing
    for bp in border_points:
        if bp in orig_set and bp not in new_set:
            # find lines that have integer-rounded x == bp[0]
            candidate_line_indices = int_to_lines.get(bp, set())
            if verbose:
                print(f"DEBUG: for missing border point {bp}, candidate line indices with x={bp[0]} ->", candidate_line_indices)
            # mark them for removal (no searching for nearest; simple rule)
            for idx in candidate_line_indices:
                removables.add(idx)
    
    if verbose:
        print("DEBUG: final set of line indices marked for removal:", removables, " count:", len(removables))
    
    # build result by filtering out removables
    result = [ln for i, ln in enumerate(new_scaled_lines) if i not in removables]
    
    stats.update({
        "missing_border_points": missing_border_points,
        "removables": sorted(list(removables)),
        "removables_count": len(removables),
        "result_lines_count": len(result)
    })
    
    return result, stats




# this is not good for CNN but makes it a good look

def filter_by_concentration(
    lines,
    safe_threshold=1,
    critical_threshold=3,
    concentration_threshold=0.75,
    coord_rounding: str = "round",   # "round"|"floor"|"ceil"|"trunc"
    return_removed: bool = False,    # if True, also return set of removed indices
):
    """
    Filter lines where another line j is concentrated inside line i's critical zone.
    A line j is removed if at least `concentration_threshold` fraction of its
    own (integerized, unique) points lie inside line i's critical zone.

    Parameters
    ----------
    lines : list[list[(x,y)]]
        Original lines (each a list of numeric (x,y) tuples). Returned filtered
        lines keep original coords unchanged.
    safe_threshold : int
    critical_threshold : int
    concentration_threshold : float (0..1)
    coord_rounding : str
        How to integerize float coords for mapping: "round","floor","ceil","trunc".
    return_removed : bool
        If True, return (filtered_lines, removables_set); else return filtered_lines.
    """
    if coord_rounding not in {"round", "floor", "ceil", "trunc"}:
        raise ValueError("coord_rounding must be one of 'round','floor','ceil','trunc'")

    def _to_int(v):
        vf = float(v)
        if coord_rounding == "round":
            return int(round(vf))
        if coord_rounding == "floor":
            return math.floor(vf)
        if coord_rounding == "ceil":
            return math.ceil(vf)
        return int(vf)  # trunc

    # integerized versions of lines (preserve order)
    int_lines = []
    unique_int_points = []
    for line in lines:
        pts = [(_to_int(x), _to_int(y)) for (x, y) in line]
        int_lines.append(pts)
        unique_int_points.append(set(pts))

    n = len(lines)
    st = int(safe_threshold)
    ct = int(critical_threshold)

    # build point -> set(indices) mapping using integerized coordinates
    point_to_lines = defaultdict(set)
    for idx, pts in enumerate(int_lines):
        for p in pts:
            point_to_lines[p].add(idx)

    # precompute critical zones (integer points) for each line i
    critical_zones = [set() for _ in range(n)]
    for idx, pts in enumerate(int_lines):
        cz = set()
        for (x, y) in pts:
            for yy in range(y - ct, y + ct + 1):
                cz.add((x, yy))
        critical_zones[idx] = cz

    removables = set()

    # iterate reference lines i
    for i in range(n):
        if i in removables:
            continue
        cz_i = critical_zones[i]
        if not cz_i:
            continue

        # count how many points of each candidate j lie inside cz_i
        counts = defaultdict(int)
        for p in cz_i:
            owners = point_to_lines.get(p)
            if not owners:
                continue
            for j in owners:
                if j == i or j in removables:
                    continue
                counts[j] += 1

        # for each candidate j, compute ratio = cnt / len(unique points of j)
        for j, cnt in counts.items():
            len_j = len(unique_int_points[j])
            if len_j == 0:
                continue  # nothing to test
            ratio = cnt / len_j
            if ratio >= concentration_threshold:
                # remove j (but never remove i)
                if j != i and j not in removables:
                    removables.add(j)
                    # update mapping: remove j from all points it occupies
                    for p in unique_int_points[j]:
                        owners = point_to_lines.get(p)
                        if owners:
                            owners.discard(j)
                            if not owners:
                                del point_to_lines[p]
                    # after removing j from mapping, it won't affect later checks

    filtered = [line for idx, line in enumerate(lines) if idx not in removables]
    if return_removed:
        return filtered, removables
    return filtered




def flip_layers_3d(lines):
    """
    Flip 3D lines (x, y, a) vertically around the maximum y value.

    Args:
        lines: list of list of (x, y, a) points, or a single list of (x, y, a).

    Returns:
        Flipped lines with y mirrored vertically.
    """
    if isinstance(lines[0], tuple):
        lines = [lines]

    # Find the maximum y-value across all layers
    max_y = max(y for line in lines for _, y, _ in line)

    # Flip y-values vertically while keeping x and a intact
    flipped_lines = [
        [(x, max_y - y, a) for x, y, a in line]
        for line in lines
    ]

    return flipped_lines if len(flipped_lines) > 1 else flipped_lines[0]



def average_y_per_x(lines):
    """
    For each line (list of (x, y) points) in 'lines':
    - Group points by identical x values
    - Replace each x group with a single point (x, average_y)
    
    Args:
        lines (list[list[tuple[int, int]]]): List of lines with (x, y) coordinates.
        
    Returns:
        list[list[tuple[float, float]]]: Modified lines where each x has one averaged y.
    """
    result = []
    for line in lines:
        if not line:
            result.append([])
            continue
        
        grouped = defaultdict(list)
        for x, y in line:
            grouped[x].append(y)
        
        # Average and reconstruct sorted by x
        averaged_line = [(x, float(np.mean(ys))) for x, ys in grouped.items()]
        averaged_line.sort(key=lambda p: p[0])  # keep x order consistent
        
        result.append(averaged_line)
    return result



def apply_mask_bias(arr, positive_ratio=1.0, negative_ratio=1.0):
    """
    Return binary masks for positive and negative values in 'arr'
    based on biased ratios of their maximum and minimum magnitudes.

    Args:
        arr (np.ndarray): 2D NumPy array with positive and negative values.
        positive_ratio (float): Fraction (0–1). Keep points >= ratio * max_positive.
        negative_ratio (float): Fraction (0–1). Keep points <= ratio * min_negative.

    Returns:
        biased_mask_p (np.ndarray): Binary mask (uint8) where positive bias applies.
        biased_mask_n (np.ndarray): Binary mask (uint8) where negative bias applies.
    """
    if arr.ndim != 2:
        raise ValueError("Input array must be 2D.")
    if not (0 <= positive_ratio <= 1 and 0 <= negative_ratio <= 1):
        raise ValueError("Ratios must be between 0 and 1.")

    # Determine thresholds
    max_pos = np.max(arr[arr > 0]) if np.any(arr > 0) else 0
    min_neg = np.min(arr[arr < 0]) if np.any(arr < 0) else 0

    # Avoid divide-by-zero behavior
    biased_mask_p = np.zeros_like(arr, dtype=np.uint8)
    biased_mask_n = np.zeros_like(arr, dtype=np.uint8)

    if max_pos > 0:
        pos_thresh = positive_ratio * max_pos
        biased_mask_p[arr >= pos_thresh] = 1

    if min_neg < 0:
        neg_thresh = negative_ratio * min_neg  # e.g., 0.5 * (-10) = -5
        biased_mask_n[arr <= neg_thresh] = 1

    return biased_mask_p, biased_mask_n



# Thank God, flood (puddle) method instead of contour

def find_puddles_cv2(mask, connectivity=8):  # obsoleted by recursive puddle technique
    """
    Find clusters (puddles) in a boolean mask using OpenCV.
    
    Args:
        mask (np.ndarray): 2D boolean array (1=True, 0=False)
        connectivity (int): 4 or 8 for connectivity
    
    Returns:
        list of lists: Each puddle is a list of (x, y) coordinates sorted by x
    """
    # Convert mask to uint8
    mask_uint8 = (mask > 0).astype(np.uint8)
    
    # Get connected components
    num_labels, labels = cv2.connectedComponents(mask_uint8, connectivity=connectivity)
    
    puddles = []
    for i in range(1, num_labels):  # label 0 is background
        coords = np.argwhere(labels == i)
        # Convert to (x, y)
        coords = [(int(x), int(y)) for y, x in coords]
        coords.sort(key=lambda p: (p[0], p[1]))
        puddles.append(coords)
    
    return puddles



def find_puddles_recursive(mask, connectivity=8):
    """
    Find all puddles (outer and nested inner puddles) by iteratively inverting
    each puddle mask inside its bounding box and extracting inner components.
    
    Args:
        mask (np.ndarray): 2D boolean or 0/1 array
        connectivity (int): 4 or 8 for cv2.connectedComponents
    
    Returns:
        list of lists: Each puddle is a list of (x, y) coordinates sorted by x then y.
    """
    mask_uint8 = (mask > 0).astype(np.uint8)

    # initial connected components on whole mask
    num_labels, labels = cv2.connectedComponents(mask_uint8, connectivity=connectivity)
    puddles = []

    # Prepare initial stack of puddle masks (full-image sized uint8 masks)
    stack = []
    seen_signatures = set()  # to avoid reprocessing same region

    for lab in range(1, num_labels):
        comp_mask = (labels == lab).astype(np.uint8)
        # signature: bounding box + area
        ys, xs = np.where(comp_mask)
        if ys.size == 0:
            continue
        ymin, ymax = int(ys.min()), int(ys.max())
        xmin, xmax = int(xs.min()), int(xs.max())
        area = int(comp_mask.sum())
        sig = (xmin, xmax, ymin, ymax, area)
        if sig in seen_signatures:
            continue
        seen_signatures.add(sig)
        stack.append((comp_mask, xmin, ymin))  # store mask and its origin (0,0 for full-image masks)

    # process stack
    while stack:
        comp_mask_full, origin_x, origin_y = stack.pop()
        # record coords for this puddle
        ys, xs = np.where(comp_mask_full)
        if ys.size == 0:
            continue
        coords = [(int(x), int(y)) for y, x in zip(ys, xs)]
        coords.sort(key=lambda p: (p[0], p[1]))
        puddles.append(coords)

        # restrict to bounding box to find inner components more cheaply
        ymin, ymax = int(ys.min()), int(ys.max())
        xmin, xmax = int(xs.min()), int(xs.max())

        sub = comp_mask_full[ymin:ymax + 1, xmin:xmax + 1].copy()
        if sub.size == 0:
            continue

        # invert submask: 1->0 and 0->1 inside bounding box
        inverted = (1 - sub).astype(np.uint8)

        # find components in inverted: these are candidate inner puddles (holes/islands)
        n2, labels2 = cv2.connectedComponents(inverted, connectivity=connectivity)

        for lab2 in range(1, n2):
            comp2 = (labels2 == lab2).astype(np.uint8)
            # check if this component touches the border of the bounding box:
            # if it touches border, it's connected to outside the puddle -> ignore
            h, w = comp2.shape
            touches_border = False
            # check top, bottom, left, right borders
            if comp2[0, :].any() or comp2[-1, :].any() or comp2[:, 0].any() or comp2[:, -1].any():
                touches_border = True

            if touches_border:
                continue

            # this comp2 is fully inside the parent's bounding box -> inner puddle
            # translate comp2 coords to full-image coords
            ys2, xs2 = np.where(comp2)
            if ys2.size == 0:
                continue
            comp2_full = np.zeros_like(comp_mask_full, dtype=np.uint8)
            # place comp2 into the corresponding bbox location
            comp2_full[ymin:ymax + 1, xmin:xmax + 1][ys2, xs2] = 1

            # signature to avoid repeats
            area2 = int(comp2.sum())
            ymin2, ymax2 = int(ys2.min()), int(ys2.max())
            xmin2, xmax2 = int(xs2.min()), int(xs2.max())
            # convert to global coords
            sig = (xmin + xmin2, xmin + xmax2, ymin + ymin2, ymin + ymax2, area2)
            if sig in seen_signatures:
                continue
            seen_signatures.add(sig)
            # push this inner puddle to stack for further processing (it may contain inner islands)
            stack.append((comp2_full, 0, 0))

    return puddles



def split_lines_by_y(lines):
    """
    For each line in lines, group points by x, record min and max y values for each x,
    and create new lines for min and max y. If min and max lines differ, mark the
    original line for removal.
    """
    new_lines = []
    removables = set()

    for idx, line in enumerate(lines):
        x_groups = defaultdict(list)
        for x, y in line:
            x_groups[x].append(y)

        x_sorted = sorted(x_groups.keys())
        min_line = [(x, min(x_groups[x])) for x in x_sorted]
        max_line = [(x, max(x_groups[x])) for x in x_sorted]

        if min_line != max_line:
            # Add both new lines
            new_lines.append(min_line)
            new_lines.append(max_line)
            removables.add(idx)
        else:
            # No split needed, keep original line
            new_lines.append(line)

    # Filter out original lines that were split
    final_lines = [line for idx, line in enumerate(lines) if idx not in removables]
    
    # Include newly created lines
    final_lines.extend([line for line in new_lines if line not in lines])
    
    return final_lines


def plot_lines_scatter_fast(lines, point_size=20, show=True):
    """
    Fast scatter plot for large number of points in lines.
    """
    plt.figure(figsize=(8,6))
    
    # Flatten all points into arrays for plotting
    all_x, all_y, all_colors = [], [], []

    cmap = plt.get_cmap('tab10')
    for i, line in enumerate(lines):
        if len(line) == 0:
            continue
        x, y = zip(*line)
        all_x.extend(x)
        all_y.extend(y)
        all_colors.extend([cmap(i % 10)] * len(line))

    # Single scatter call
    plt.scatter(all_x, all_y, c=all_colors, s=point_size)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title("Scatter plot of lines")
    plt.grid(True)
    if show:
        plt.show()


# --- GUI App Class ---

def plot_layers_GUI(layers, marked_layers=False, pixels=512, line_color='red'):
    """Return a matplotlib Figure and Axes with layers and optional marked layers plotted."""
    # Normalize input
    if isinstance(layers[0], tuple):
        layers = [layers]
    if marked_layers and isinstance(marked_layers[0], tuple):
        marked_layers = [marked_layers]

    dpi = 100
    fig_size = pixels / dpi
    fig = Figure(figsize=(fig_size, fig_size), dpi=dpi)
    ax = fig.add_subplot(111)

    all_x_coords = []
    all_y_coords = []

    def process_layer(layer):
        unique_pts = list(set(layer))
        return sorted(unique_pts, key=lambda pt: pt[0])

    # Plot regular layers in black
    for layer in layers:
        cleaned = process_layer(layer)
        if cleaned:
            x = [pt[0] for pt in cleaned]
            y = [pt[1] for pt in cleaned]
            all_x_coords.extend(x)
            all_y_coords.extend(y)
            ax.plot(x, y, c='black', linewidth=1)

    # Plot marked layers in red
    if marked_layers:
        for layer in marked_layers:
            cleaned = process_layer(layer)
            if cleaned:
                x = [pt[0] for pt in cleaned]
                y = [pt[1] for pt in cleaned]
                all_x_coords.extend(x)
                all_y_coords.extend(y)
                ax.plot(x, y, c=line_color, linewidth=2)

    # Set axis limits
    if all_x_coords and all_y_coords:
        ax.set_xlim(min(all_x_coords), max(all_x_coords))
        ax.set_ylim(min(all_y_coords), max(all_y_coords))

    ax.axis('off')
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    return fig, ax


def plot_seismic_GUI(seismic_profile_transposed, marked_points=False, pixels=512, relief_seismic = 99):
    dpi = 100
    fig_size = pixels / dpi
    fig = Figure(figsize=(fig_size, fig_size), dpi=dpi)
    ax = fig.add_subplot(111)

    relief_seismic = np.clip(relief_seismic, 51, 100)
    vm = np.percentile(seismic_profile_transposed, relief_seismic)
    ax.imshow(seismic_profile_transposed, cmap="RdBu", vmin=-vm, vmax=vm, aspect='auto')

    if marked_points:
        if isinstance(marked_points[0], (int, float)):
            marked_points = [marked_points]
        for x, y in marked_points:
            ax.plot(x, y, 'ro', markersize=4)

    ax.axis('off')
    return fig, ax


# Thank God!! finally twt and cdp representation is correct
def plot_seismic_CROPPED(
    seismic_profile,      # MUST be shape (n_samples, n_traces) -> rows=time, cols=traces/CDP
    pixels=512,
    marked_points=False,  # list of (x, y) points in CDP/TWT space (y in ms if twt_spacing in ms)
    axis_off=True,
    seismic_relief=99,
    inline=-999,
    twt_spacing=4,        # ms per sample
    x_start=0,            # CDP index of left edge of this crop
    y_start=0,             # TWT (ms) at the top of this crop (NOT sample index)
    initial_inline=0,
    crossline=0,    # to specify crossline, make this one 1 and change inline/initial_inline as desired crossline's parameters
    save_option=True,
    close_option=True,
    fixed_name="Seismic Profile",
    return_fig=False,
    svg_option=False
):
    """Plot seismic_profile where rows=time samples and cols=CDP/traces.
       x_start is CDP index; y_start is TWT in same units as twt_spacing (e.g., ms).
    """
    # Validate input
    if not isinstance(seismic_profile, np.ndarray) or seismic_profile.ndim != 2:
        raise ValueError("seismic_profile must be a 2D numpy array with shape (n_samples, n_traces)")

    # compute figure size
    dpi = 100
    if isinstance(pixels, int):
        width_px, height_px = pixels, pixels
    elif isinstance(pixels, (tuple, list)) and len(pixels) == 2:
        width_px, height_px = pixels
    else:
        raise ValueError("pixels must be an int or a tuple of (width, height)")
    fig = plt.figure(figsize=(width_px / dpi, height_px / dpi), dpi=dpi)
    ax = fig.add_subplot(1, 1, 1)

    # amplitude clipping
    seismic_relief = np.clip(seismic_relief, 51, 100)
    vm = np.percentile(seismic_profile, seismic_relief)

    # dimensions (no transpose here: rows=time, cols=traces)
    n_traces, n_samples = seismic_profile.shape
    twt_height = n_samples * twt_spacing   # ms covered by this crop
    x_end = x_start + n_traces
    y_end = y_start + twt_height
    print("CDP:", n_traces, "twt samples:", n_samples, "x_start:", x_start, "y_start:", y_start)
    # show image mapped into absolute CDP/TWT coords
    im = ax.imshow(
        seismic_profile.T,
        cmap="RdBu",
        vmin=-vm,
        vmax=vm,
        aspect='auto',
        extent=[x_start, x_end, y_start, y_end],
        origin='upper'   # top-left origin (visual: TWT increases downward)
    )

    # overlay points (already in CDP, TWT(ms))
    if marked_points:
        if isinstance(marked_points[0], (int, float)):
            marked_points = [marked_points]
        xs, ys = zip(*marked_points)
        ax.plot(xs, ys, 'ro', markersize=4)

    # axis / ticks
    if axis_off:
        ax.axis('off')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    else:
        ax.set_xlabel('CDP number')
        ax.set_ylabel('TWT [ms]')
        if crossline == 0:
            inline_plot = inline + initial_inline
            ax.set_title(f"{fixed_name} inline = {inline_plot}")
        else:
            inline_plot = inline + initial_inline
            ax.set_title(f"{fixed_name} crossline = {inline_plot}")            

        num_ticks = int(width_px/85)  # 6 for 512 px
        twt_ticks = np.linspace(y_start, y_end, num_ticks)
        cdp_ticks = np.linspace(x_start, x_end, num_ticks)
        
        # Round to nearest integer and convert to int
        twt_ticks = np.round(twt_ticks).astype(int)
        cdp_ticks = np.round(cdp_ticks).astype(int)        

        ax.set_yticks(twt_ticks)      # TWT on y-axis
        ax.set_xticks(cdp_ticks)      # CDP on x-axis (fixed)

        # show labels inverted top->bottom as you prefer (top label = y_end, bottom = y_start)
        twt_labels = [f"{int(v)}" for v in twt_ticks]
        ax.set_yticklabels(list(reversed(twt_labels)))

    
        
    if save_option:
        inline_plot = inline + initial_inline
        if crossline == 0:
            plt.savefig(f"{fixed_name}_inline_{inline_plot}.png", bbox_inches='tight', pad_inches=0, dpi=dpi)
            if svg_option:
                plt.savefig(f"{fixed_name}_inline_{inline_plot}.svg", format="svg", bbox_inches='tight', pad_inches=0, dpi=dpi)            
                #plt.savefig(f"{fixed_name}_inline_{inline}.eps", format="eps", bbox_inches='tight', pad_inches=0, dpi=dpi)            
        else:
            plt.savefig(f"{fixed_name}_crossline_{inline_plot}.png", bbox_inches='tight', pad_inches=0, dpi=dpi)
            if svg_option:
                plt.savefig(f"{fixed_name}_crossline_{inline_plot}.svg", format="svg", bbox_inches='tight', pad_inches=0, dpi=dpi)            
                #plt.savefig(f"{fixed_name}_crossline_{inline}.eps", format="eps", bbox_inches='tight', pad_inches=0, dpi=dpi)            
            
            
    if return_fig:
        return fig            
            
    if not close_option:
        plt.show()

    plt.close(fig)




    
####affine  => some of functions below need reconsidering, I just mentioned them because GUI might need them
def best_affine_matches(lines, aff_ref, transform_type="rigid", score_type="hausdorff", n_resample=50):
    """
    For each ref_line in aff_ref, find the best-matching line in `lines`
    under the chosen transformation type and scoring method.

    Parameters
    ----------
    lines : list of list of (x,y)
        Candidate lines.
    aff_ref : list of list of (x,y)
        Reference lines.
    transform_type : {"rigid","similarity","affine"}
        Type of transformation to use.
    score_type : {"mse","hausdorff"}
        Error metric to evaluate similarity.
    n_resample : int
        Number of points to resample each line to before scoring (for fairness).

    Returns
    -------
    result : list
        Best-matching line for each ref_line.
    """

    def resample_line(line, n=n_resample):
        """Resample line to n evenly spaced points along arc length."""
        line = np.array(line, dtype=float)
        if len(line) < 2:
            return line
        dists = np.cumsum(np.r_[0, np.linalg.norm(np.diff(line, axis=0), axis=1)])
        uniform_dists = np.linspace(0, dists[-1], n)
        new_points = np.zeros((n, 2))
        for i in range(2):  # x and y separately
            new_points[:, i] = np.interp(uniform_dists, dists, line[:, i])
        return new_points

    def estimate_transform(src, dst, mode="affine"):
        """Estimate transform matrix src->dst."""
        src = np.array(src[:3], dtype=float)
        dst = np.array(dst[:3], dtype=float)

        if mode == "affine":
            A, B = [], []
            for (x, y), (u, v) in zip(src, dst):
                A.append([x, y, 1, 0, 0, 0])
                A.append([0, 0, 0, x, y, 1])
                B.append(u)
                B.append(v)
            A, B = np.array(A), np.array(B)
            params, *_ = np.linalg.lstsq(A, B, rcond=None)
            M = np.array([[params[0], params[1], params[2]],
                          [params[3], params[4], params[5]]])

        elif mode in {"rigid", "similarity"}:
            # Procrustes alignment
            src_mean, dst_mean = src.mean(0), dst.mean(0)
            src_c, dst_c = src - src_mean, dst - dst_mean
            U, _, Vt = np.linalg.svd(src_c.T @ dst_c)
            R = Vt.T @ U.T
            if np.linalg.det(R) < 0:  # reflection fix
                Vt[-1, :] *= -1
                R = Vt.T @ U.T
            if mode == "similarity":
                scale = np.trace((dst_c @ R).T @ src_c) / np.sum(src_c**2)
            else:
                scale = 1.0
            M = np.array([[scale * R[0, 0], scale * R[0, 1], dst_mean[0] - (scale * R @ src_mean)[0]],
                          [scale * R[1, 0], scale * R[1, 1], dst_mean[1] - (scale * R @ src_mean)[1]]])
        else:
            raise ValueError("Invalid mode")
        return M

    def apply_transform(M, pts):
        pts = np.array(pts, dtype=float)
        ones = np.ones((pts.shape[0], 1))
        return np.hstack([pts, ones]) @ M.T

    def score(transformed, target, method="mse"):
        if method == "mse":
            return np.mean(np.linalg.norm(transformed - target, axis=1))
        elif method == "hausdorff":
            return max(directed_hausdorff(transformed, target)[0],
                       directed_hausdorff(target, transformed)[0])
        else:
            raise ValueError("Invalid score method")

    result = []
    for ref_line in aff_ref:
        best_err, best_line = float("inf"), []
        if len(ref_line) < 3:
            result.append([])
            continue
        for line in lines:
            if len(line) < 3:
                continue
            ref_res = resample_line(ref_line)
            line_res = resample_line(line)
            M = estimate_transform(ref_res, line_res, mode=transform_type)
            ref_t = apply_transform(M, ref_res)
            err = score(ref_t, line_res, method=score_type)
            if err < best_err:
                best_err, best_line = err, line
        result.append(best_line)
    return result


def rescale_2d_array(profile, target_size=(512, 512), interpolation='linear'):
    """
    Rescale a 2D numpy array to a target size.

    Parameters
    ----------
    profile : 2D numpy array
        Input array to rescale.
    target_size : tuple of int
        (T1, T2) target shape (rows, cols). Default (512, 512).
    interpolation : {"nearest","linear","cubic"}
        Interpolation method. Default "linear".

    Returns
    -------
    rescaled : 2D numpy array
        Rescaled array with shape target_size.
    """

    if not isinstance(profile, np.ndarray) or profile.ndim != 2:
        raise ValueError("profile must be a 2D numpy array")

    interp_map = {
        'nearest': cv2.INTER_NEAREST,
        'linear': cv2.INTER_LINEAR,
        'cubic': cv2.INTER_CUBIC
    }

    if interpolation not in interp_map:
        raise ValueError(f"interpolation must be one of {list(interp_map.keys())}")

    T2, T1 = target_size[1], target_size[0]  # OpenCV uses (width, height)
    rescaled = cv2.resize(profile, (T2, T1), interpolation=interp_map[interpolation])
    return rescaled



def produce_aff_3d(
    cropped,
    aff_2d,
    clip=False,
    prefer_transpose_check=True,
    verbose=True,
    return_report=False,
    example_limit=10
):
    """
    Robustly produce (x,y,a) lines from a 2D array and list of (x,y) polylines.

    This function tries to detect whether the input coordinates match the array
    indexing convention (a[row, col] == a[y, x]) or are swapped (a[x, y]).
    It returns lines with amplitude values; out-of-bounds amplitude entries
    are set to np.nan (not None).

    Parameters
    ----------
    cropped : np.ndarray
        2D numpy array (shape: rows=height, cols=width).
    aff_2d : list of list of (x, y)
        Polylines: each line is a list of (x, y) floats or ints.
    clip : bool
        If True, clip coordinates to valid integer indices instead of setting np.nan.
    prefer_transpose_check : bool
        If True, the function will test both (x,y)->(col,row) and swapped (y,x)
        mapping and pick the one that yields fewer out-of-bounds points.
    verbose : bool
        Print a short diagnostic summary.
    return_report : bool
        If True return (aff_3d, report). Otherwise return aff_3d only.
    example_limit : int
        How many out-of-bounds examples to include in the report/print.

    Returns
    -------
    aff_3d : list of list of (x, y, a)
        Polylines with amplitude values (a is float or np.nan).
    (optional) report : dict
        Diagnostic report (only if return_report=True).
    """
    # --- input validation
    if not isinstance(cropped, np.ndarray) or cropped.ndim != 2:
        raise ValueError("cropped must be a 2D numpy array")

    if not isinstance(aff_2d, (list, tuple)):
        raise ValueError("aff_2d must be a list of polylines (list of (x,y) tuples)")

    height, width = cropped.shape

    # Flatten coordinate stats for diagnostics
    all_pts = [(x, y) for line in aff_2d for (x, y) in line]
    if len(all_pts) == 0:
        report = {
            "mapping_used": None,
            "total_points": 0,
            "oob_count": 0,
            "oob_examples": [],
            "aff_min_max": None,
            "cropped_shape": (height, width),
            "message": "No points in aff_2d."
        }
        if verbose:
            print("produce_aff_3d_debug: input aff_2d empty.")
        return ([], report) if return_report else []

    xs = np.array([p[0] for p in all_pts], dtype=float)
    ys = np.array([p[1] for p in all_pts], dtype=float)

    aff_min_max = {"x_min": float(xs.min()), "x_max": float(xs.max()),
                   "y_min": float(ys.min()), "y_max": float(ys.max())}

    # --- helper to compute OOB count for a mapping function
    def oob_count_for_mapping(mapping='xy'):
        cnt = 0
        for x, y in all_pts:
            if mapping == 'xy':
                xi, yi = int(round(x)), int(round(y))
            else:  # 'yx' swapped
                xi, yi = int(round(y)), int(round(x))
            if not (0 <= xi < width and 0 <= yi < height):
                cnt += 1
        return cnt

    # Try both mappings if requested
    mapping_used = 'xy'
    if prefer_transpose_check:
        cnt_xy = oob_count_for_mapping('xy')
        cnt_yx = oob_count_for_mapping('yx')
        # choose the mapping with fewer OOBs
        if cnt_yx < cnt_xy:
            mapping_used = 'yx'
            oob_cnt = cnt_yx
        else:
            mapping_used = 'xy'
            oob_cnt = cnt_xy
    else:
        oob_cnt = oob_count_for_mapping('xy')

    # If many OOBs, print a hint
    total_points = len(all_pts)
    oob_fraction = oob_cnt / total_points

    if verbose:
        print(f"produce_aff_3d_debug: cropped.shape={cropped.shape}")
        print(f"produce_aff_3d_debug: aff_2d points: {total_points}, mapping chosen: '{mapping_used}'")
        print(f"produce_aff_3d_debug: coordinate ranges: x[{aff_min_max['x_min']:.2f},{aff_min_max['x_max']:.2f}], "
              f"y[{aff_min_max['y_min']:.2f},{aff_min_max['y_max']:.2f}]")
        print(f"produce_aff_3d_debug: out-of-bounds count={oob_cnt} ({oob_fraction*100:.2f}%)")
        if oob_cnt > 0:
            hint = ("Coordinates may be in a different orientation (x/y swapped), "
                    "or contain values outside the cropped extent. "
                    "Check whether you used cropped.T earlier.")
            print("Hint:", hint)

    # --- build the result and collect examples of OOBs
    aff_3d = []
    oob_examples = []
    for line in aff_2d:
        line_with_values = []
        for x, y in line:
            # choose mapping
            if mapping_used == 'xy':
                xi, yi = int(round(x)), int(round(y))
                x_out, y_out = x, y
            else:
                xi, yi = int(round(y)), int(round(x))
                # keep original semantics in returned tuple (x,y,a): return original x,y
                x_out, y_out = x, y

            # check bounds
            if 0 <= xi < width and 0 <= yi < height:
                a = float(cropped[yi, xi])  # row = yi, col = xi
            else:
                if clip:
                    # clip into bounds and fetch value
                    xi_clipped = min(max(xi, 0), width - 1)
                    yi_clipped = min(max(yi, 0), height - 1)
                    a = float(cropped[yi_clipped, xi_clipped])
                    # keep example if it was OOB before clipping
                    if len(oob_examples) < example_limit:
                        oob_examples.append({
                            "orig": (x, y),
                            "mapped_idx": (xi, yi),
                            "clipped_idx": (xi_clipped, yi_clipped),
                            "value_used": a
                        })
                else:
                    a = np.nan
                    if len(oob_examples) < example_limit:
                        oob_examples.append({
                            "orig": (x, y),
                            "mapped_idx": (xi, yi),
                            "note": "out_of_bounds"
                        })
            line_with_values.append((x_out, y_out, a))
        aff_3d.append(line_with_values)

    report = {
        "mapping_used": mapping_used,
        "total_points": total_points,
        "oob_count": oob_cnt,
        "oob_fraction": oob_fraction,
        "oob_examples": oob_examples,
        "aff_min_max": aff_min_max,
        "cropped_shape": (height, width),
        "clip": clip
    }

    if verbose:
        print("produce_aff_3d_debug: report summary:")
        print(f"  mapping_used = {mapping_used}")
        print(f"  total_points = {total_points}, oob_count = {oob_cnt}")
        if oob_examples:
            print("  example OOBs (up to {}):".format(len(oob_examples)))
            for ex in oob_examples:
                print("   ", ex)
        else:
            print("  no example OOBs to show (all points in bounds)")

    return (aff_3d, report) if return_report else aff_3d



def lines_2d_to_3d(cropped, lines_2d):
    """
    Convert 2D lines to 3D lines with amplitude.

    Parameters
    ----------
    cropped : 2D numpy array
        The seismic data array.
    lines_2d : list of list of (x,y)
        Each sublist is a line of (x,y) coordinates.

    Returns
    -------
    lines_3d : list of list of (x,y,a)
        Each line now contains (x,y,amplitude) points.
    """

    if not isinstance(cropped, np.ndarray) or cropped.ndim != 2:
        raise ValueError("cropped must be a 2D numpy array")

    lines_3d = []
    nrows, ncols = cropped.shape

    for line in lines_2d:
        line_3d = []
        for x, y in line:
            # Clip coordinates to valid indices
            xi = int(round(x))
            yi = int(round(y))
            xi = max(0, min(ncols - 1, xi))
            yi = max(0, min(nrows - 1, yi))
            a = cropped[yi, xi]
            line_3d.append((x, y, a))
        lines_3d.append(line_3d)

    return lines_3d



def best_line_matches_3d(lines_3d, ref_lines, 
                         transform_type="affine", 
                         score_type="mse", 
                         n_resample=50, 
                         amplitude_normalization=False):
    """
    Match 3D reference lines to candidate 3D lines using affine/similarity/rigid transforms.

    Parameters
    ----------
    lines : list of list of (x,y,a)
        Candidate 3D lines.
    ref_lines : list of list of (x,y,a)
        Reference 3D lines.
    transform_type : {"rigid","similarity","affine"}
        Transformation type.
    score_type : {"mse","hausdorff"}
        Error metric.
    n_resample : int
        Number of resampled points for fair comparison.
    amplitude_normalization : bool
        If True, normalize amplitude range to match y-range.

    Returns
    -------
    result : list
        Best-matching line for each reference line.
    """

    def normalize_amplitude(line):
        line = np.array(line, dtype=float)
        if len(line) == 0:
            return line
        x, y, a = line[:,0], line[:,1], line[:,2]
        yrange = np.ptp(y)
        arange = np.ptp(a)
        if amplitude_normalization and arange > 0:
            scale = yrange / arange
            a = a * scale
        return np.column_stack([x, y, a])


    def resample_line(line, n=n_resample):
        """
        Resample a 3D line to n evenly spaced points along arc length,
        guaranteeing that resampled points lie on the original polyline.

        Parameters
        ----------
        line : list of (x, y, a)
            Original 3D line as a sequence of points.
        n : int
            Number of resampled points.

        Returns
        -------
        np.ndarray of shape (n, 3)
            Resampled 3D points lying on the polyline.
        """
        line = np.array(line, dtype=float)
        if len(line) < 2:
            return line

        # Compute cumulative arc-lengths of the polyline
        seg_lengths = np.linalg.norm(np.diff(line, axis=0), axis=1)
        cum_dists = np.r_[0, np.cumsum(seg_lengths)]
        total_length = cum_dists[-1]

        # Uniform target distances
        target_dists = np.linspace(0, total_length, n)

        # Output points
        new_points = []

        # Walk through each target distance and find segment it belongs to
        seg_idx = 0
        for td in target_dists:
            # Move forward until correct segment found
            while td > cum_dists[seg_idx+1]:
                seg_idx += 1

            # Normalize within segment
            t = (td - cum_dists[seg_idx]) / seg_lengths[seg_idx]
            p0, p1 = line[seg_idx], line[seg_idx+1]

            # Linear interpolation along segment (stays on polyline)
            new_point = (1 - t) * p0 + t * p1
            new_points.append(new_point)

        return np.array(new_points)


    def estimate_transform(src, dst, mode="rigid"):
        """Estimate rigid/similarity/affine transform between 3D src and dst lines."""
        # Safety check
        if np.any(~np.isfinite(src)) or np.any(~np.isfinite(dst)):
            return np.eye(src.shape[1] + 1)

        # Center
        src_mean, dst_mean = src.mean(0), dst.mean(0)
        src_c, dst_c = src - src_mean, dst - dst_mean

        # Check rank (avoid collinear or too few points)
        if min(src_c.shape) < 2 or np.linalg.matrix_rank(src_c) < 2:
            return np.eye(src.shape[1] + 1)

        if mode == "rigid":
            # Standard Orthogonal Procrustes
            H = src_c.T @ dst_c
            H += 1e-12 * np.eye(H.shape[0])  # jitter for SVD stability
            try:
                U, _, Vt = np.linalg.svd(H)
            except np.linalg.LinAlgError:
                return np.eye(src.shape[1] + 1)  # fallback: identity

            R = Vt.T @ U.T
            if np.linalg.det(R) < 0:  # reflection fix
                Vt[-1, :] *= -1
                R = Vt.T @ U.T
            t = dst_mean - R @ src_mean
            M = np.eye(src.shape[1] + 1)
            M[:R.shape[0], :R.shape[1]] = R
            M[:-1, -1] = t
            return M

        elif mode == "similarity":
            # Use Umeyama algorithm for similarity transform
            cov = dst_c.T @ src_c / len(src)
            U, S, Vt = np.linalg.svd(cov)
            R = U @ Vt
            if np.linalg.det(R) < 0:
                Vt[-1, :] *= -1
                R = U @ Vt

            var_src = np.var(src_c, axis=0).sum()
            scale = np.trace(np.diag(S)) / var_src
            t = dst_mean - scale * R @ src_mean

            M = np.eye(src.shape[1] + 1)
            M[:R.shape[0], :R.shape[1]] = scale * R
            M[:-1, -1] = t
            return M

        elif mode == "affine":
            A = np.c_[src, np.ones(len(src))]
            B, _, _, _ = np.linalg.lstsq(A, dst, rcond=None)
            M = np.eye(src.shape[1] + 1)
            M[:-1, :] = B.T
            return M

        else:
            raise ValueError(f"Unknown transform mode: {mode}")


    def apply_transform(M, pts):
        pts = np.array(pts, dtype=float)
        ones = np.ones((pts.shape[0], 1))
        transformed = np.hstack([pts, ones]) @ M.T
        transformed = transformed[:, :-1]  # drop homogeneous column
        # Debugging
#         print("apply_transform -> input shape:", pts.shape, "output shape:", transformed.shape)
        return transformed


    def score(transformed, target, method="mse"):
#         print("score() -> transformed:", transformed.shape, "target:", target.shape)
        if method == "mse":
            return np.mean(np.linalg.norm(transformed - target, axis=1))
        elif method == "hausdorff":
            return max(directed_hausdorff(transformed, target)[0],
                       directed_hausdorff(target, transformed)[0])
        else:
            raise ValueError("Invalid score method")


    result = []
    for ref_line in ref_lines:
        ref_line = normalize_amplitude(ref_line)
        if len(ref_line) < 3:
            result.append([])
            continue
        best_err, best_line = float("inf"), []
        for line in lines_3d:
            if len(line) < 3:
                continue
            ref_res = resample_line(ref_line)
            line_res = resample_line(line)
            M = estimate_transform(ref_res, line_res, mode=transform_type)
            ref_t = apply_transform(M, ref_res)
            err = score(ref_t, line_res, method=score_type)
            if err < best_err:
                best_err, best_line = err, line
        result.append(best_line)
    return result



# new lines_3d codes for mapping
# Thank God, overlay function was fixed

# Assuming flip_layers and flip_layers_3d are defined elsewhere
# def flip_layers(lines): ...
# def flip_layers_3d(lines): ...
# Thank God; fixed by Chatgpt
def map_amplitude_to_lines(
    seismic_profile,
    lines,
    pixels=(512, 512),
    seismic_relief=99,
    inline=-990,
    overlay_debug=True,
    twt_spacing=4,
    axis_off=True,
    x_start=0,
    y_start=0,
    initial_inline=0,
    crossline=0,
    line_color='lime',
    return_fig=False,
    close_option=False,
    svg_option=False
):

    dpi = 100

    # -------------------------------
    # Input validation
    # -------------------------------
    if not isinstance(seismic_profile, np.ndarray) or seismic_profile.ndim != 2:
        raise ValueError("seismic_profile must be a 2D numpy array")

    if isinstance(pixels, int):
        width_px, height_px = pixels, pixels
    else:
        width_px, height_px = pixels

    seismic_relief = np.clip(seismic_relief, 51, 100)
    vm = np.percentile(seismic_profile, seismic_relief)

    # seismic geometry
    n_traces, n_samples = seismic_profile.shape
    x_end = x_start + n_traces
    y_end = y_start + n_samples * twt_spacing

    # -------------------------------
    # Figure
    # -------------------------------
    fig, ax = plt.subplots(figsize=(width_px / dpi, height_px / dpi), dpi=dpi)

    ax.imshow(
        seismic_profile.T,
        cmap='RdBu',
        vmin=-vm,
        vmax=vm,
        aspect='auto',
        extent=[x_start, x_end, y_start, y_end],
        origin='upper'
    )

    # -------------------------------
    # Normalize lines
    # -------------------------------
    if isinstance(lines[0], tuple):
        lines = [lines]

    # -------------------------------
    # Scale lines to seismic world space
    # -------------------------------
    scaled_lines = rescale_lines_to_image(
        lines,
        dimensions=(x_end - x_start, y_end - y_start),
        offset=(x_start, y_start)
    )

    if isinstance(scaled_lines[0], tuple):
        scaled_lines = [scaled_lines]

    # -------------------------------
    # Apply SAME transform as plot_layers
    # -------------------------------
    max_y = max(y for ln in scaled_lines for _, y in ln)


#######   transform seemed unnecessary with y_start and y_end, especially flip_layers made it false
    # transformed = []
    # for ln in scaled_lines:
        # tmp = [(x + x_start, max_y - y + y_start) for x, y in ln]
        # transformed.append(tmp)

    #transformed = flip_layers(transformed)
    
    transformed = scaled_lines    

    # -------------------------------
    # Plot lines (ALWAYS)
    # -------------------------------
    for layer in transformed:
        x = [p[0] for p in layer]
        y = [p[1] for p in layer]
        ax.plot(x, y, color=line_color, linewidth=1.5)

    # -------------------------------
    # Axis handling
    # -------------------------------
    if axis_off:
        ax.axis('off')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    else:
        ax.set_xlabel("CDP number")
        ax.set_ylabel("TWT [ms]")

        inline_plot = inline + initial_inline
        title = "inline" if crossline == 0 else "crossline"
        ax.set_title(f"Seismic Profile {title} = {inline_plot}")

        num_ticks = 6
        twt_ticks = np.linspace(y_start, y_end, num_ticks)
        ax.set_yticks(twt_ticks)
        ax.set_yticklabels(list(reversed([f"{int(v)}" for v in twt_ticks])))

    # -------------------------------
    # Save
    # -------------------------------
    plt.savefig(
        f"seismic_plot_mapped_{inline}.png",
        bbox_inches='tight',
        pad_inches=0,
        dpi=dpi
    )
    
    if svg_option:
        plt.savefig(
            f"seismic_plot_mapped_{inline}.svg",
            bbox_inches='tight',
            pad_inches=0,
            dpi=dpi
        )

    if return_fig:
        return fig

    if not close_option:
        plt.show()

    plt.close(fig)





def seismic_mask_to_lines(seismic_profile, pixels=(512, 512)):
    """
    Convert a binary seismic mask into a list of single-point lines [(x, y)]
    for all nonzero pixels.
    Returns:
        lines: list of single-point lines, e.g. [[(x, y)], [(x, y)], ...]
    """
    if not isinstance(seismic_profile, np.ndarray) or seismic_profile.ndim != 2:
        raise ValueError("seismic_profile must be a 2D numpy array")

    height, width = seismic_profile.shape
    width_px, height_px = pixels if isinstance(pixels, (tuple, list)) else (pixels, pixels)

    # Find nonzero (active) coordinates in the binary mask
    y_indices, x_indices = np.nonzero(seismic_profile)

    # Scaling factors to map array indices to pixel coordinates
    x_scale = (width_px - 1) / (width - 1) if width > 1 else 1
    y_scale = (height_px - 1) / (height - 1) if height > 1 else 1

    # Create single-point lines for all active pixels
    lines = [[(int(x * x_scale), int(y * y_scale))] for x, y in zip(x_indices, y_indices)]

    return lines


def get_polarity_alternation_mask(arr):  # produces noisy outcome, obsoleted
    """
    Mimic digofile() logic to detect both positive→negative and negative→positive transitions.
    """
    # Treat zero as positive like arrRB
    arr_pos = np.where(arr >= 0, arr, 0)
    arr_neg = np.where(arr < 0, np.abs(arr), 0)

    # Shift vertically (compare each sample with next)
    arr_pos_next = np.roll(arr_pos, -1, axis=0)
    arr_neg_next = np.roll(arr_neg, -1, axis=0)

    # digofile logic, both sides combined
    mask_neg_to_pos = (arr_neg == 0) & (arr_pos_next > 0)
    mask_pos_to_neg = (arr_neg > 0) & (arr_pos_next == 0)

    alternation_mask = (mask_neg_to_pos | mask_pos_to_neg).astype(np.uint8)
    alternation_mask[-1, :] = 0  # remove rolled wrap-around

    return alternation_mask



def array_to_coordinate_lists(arr):
    """
    Convert a 2D array to a list of single-point coordinate lists (x, y)
    where nonzero pixels exist after scaling to 0–255.
    """
    # Scale to 0–255 and convert to uint8 binary image
    binary_image = (arr * 255).astype(np.uint8)
    
    # Get coordinates of nonzero pixels
    y_indices, x_indices = np.nonzero(binary_image)
    
    # Build list of lists of (x, y)
#     coords = [[(int(x), int(y))] for x, y in zip(x_indices, y_indices)] 
    coords = [[(int(x), int(y)), (int(x)+1, int(y))] for x, y in zip(x_indices, y_indices)]    # two of the same y-point seems better for our later functions
    
    return coords



def magnet_mask(alternation_mask, threshold=10, distance_threshold=1):
    """
    Merge small clusters (nonzero regions) in a binary mask that are within a 
    given distance threshold from each other.

    Args:
        alternation_mask (np.ndarray): 2D binary array where nonzero pixels indicate clusters.
        threshold (int): Maximum cluster size (in pixels) to consider for merging.
        distance_threshold (int): Max pixel distance to merge nearby clusters.

    Returns:
        list of np.ndarray: List of binary masks for merged clusters.
    """
    # Ensure binary
    mask = (alternation_mask > 0).astype(np.uint8)

    # Define neighborhood structure (for connectivity)
    struct = generate_binary_structure(2, 2)

    changed = True
    while changed:
        changed = False

        # Label connected components
        labeled, num = label(mask, structure=struct)
        if num == 0:
            return []

        # Extract small clusters
        cluster_masks = []
        for i in range(1, num + 1):
            cluster = (labeled == i)
            if cluster.sum() < threshold:
                cluster_masks.append(cluster)

        # Try to merge small clusters close to each other
        merged = []
        for i in range(len(cluster_masks)):
            if cluster_masks[i] is None:
                continue

            a = cluster_masks[i]
            expanded_a = np.pad(a, distance_threshold)
            merged_mask = expanded_a.copy()

            for j in range(i + 1, len(cluster_masks)):
                b = cluster_masks[j]
                if b is None:
                    continue

                # Pad b equally for overlap check
                padded_b = np.pad(b, distance_threshold)
                # Check if expansions overlap
                if np.any(expanded_a & padded_b):
                    merged_mask |= padded_b
                    cluster_masks[j] = None
                    changed = True

            merged.append(merged_mask[distance_threshold:-distance_threshold, distance_threshold:-distance_threshold])

        # Rebuild mask for next iteration
        mask = np.zeros_like(mask)
        for m in merged:
            mask |= m.astype(np.uint8)

    # Return final merged clusters
    labeled, num = label(mask, structure=struct)
    return [(labeled == i).astype(np.uint8) for i in range(1, num + 1)]




##### classic method but very detailed contour ==> seems this much detail is bad for tsp ! outlines puddle contours save structure of 
## layers better

def digofile(compacted_image, negative_sign = False, DS_option = False):
    '''negative_sign shows if you would like to detect negative layers in compacted image (negative seismic amplitude)'''
    T,L,C = np.shape(compacted_image)
    digofile = np.zeros((T,L), dtype = int)
    DS = []
    for k in range(L):
        col = []
        for t in range(T-1):
            if (compacted_image[t,k,2] == 0 and compacted_image[t+1,k,0] > 0) and negative_sign: # [t,k,0] is red [t,k,2] is blue
                col += [t]
                digofile[t,k] = 1 #(top blue)
            elif (compacted_image[t,k,2] > 0 and compacted_image[t+1,k,0] == 0) and not negative_sign: # tracking the other side of the blue layer
                col += [t]
                digofile[t,k] = 1 #(top red)

        DS += [col]
        
    return digofile

    # if DS_option:
        # return digofile, DS
    #########

    
def arrRB(arr):
    """
    Convert a 2D array of size FxS into a 3D array of size FxSx3, 
    where the third dimension defines the RGB channels.
    
    - Red channel: positive values or zero from the input array.
    - Green channel: all values are set to 0.
    - Blue channel: absolute value of the negative values from the input array.
    
    Parameters:
    arr (np.array): 2D input array of size FxS.
    
    Returns:
    np.array: 3D array of size FxSx3 representing the RGB channels.
    """
    # Get the dimensions of the input array
    F, S = arr.shape
    
    # Initialize the 3D array with zeros (RGB)
    rgb_array = np.zeros((F, S, 3), dtype=np.float32)
    
    # Assign positive or zero values to the Red channel
    rgb_array[:, :, 0] = np.where(arr >= 0, arr, 0)
    
    # Green channel remains zeros
    
    # Assign absolute values of negative values to the Blue channel
    rgb_array[:, :, 2] = np.where(arr < 0, np.abs(arr), 0)
    
    return rgb_array



####

# intrusive analysis => contour
def exclude_sub_lines(lines_main, lines_sub):
    lines_sub_set = {tuple(map(tuple, line)) for line in lines_sub if line}
    return [line for line in lines_main if tuple(map(tuple, line)) not in lines_sub_set]
    
    
def custome_floodfill_process(stripe_coordinates_n, arr, length_limit):
    stripe_coordinates_n = remove_subset_lines(stripe_coordinates_n)
    stripe_coordinates_n = split_lines_by_y(stripe_coordinates_n)
    stripe_coordinates_n = integer_fitting(stripe_coordinates_n)
    stripe_coordinates_n = flip_layers(stripe_coordinates_n)
    min_value = min(arr.shape[0], arr.shape[1])
    stripe_coordinates_n = filter_lines_by_length(
        stripe_coordinates_n, min(min_value, length_limit)
    )
    return stripe_coordinates_n
    
    
def plot_outlines_and_blobs(outlines,
                            blobs,
                            pixels=512,
                            color1='green',
                            color2='blue',
                            linewidth1 = 1,
                            linestyle1 = '-',
                            linewidth2 = 1,
                            linestyle2 = '-',
                            x_start = 0,
                            y_stary = 0,
                            save_path="outlines_and_blobs"):
    """
    Plot outlines (green) and blobs (blue), scaled to a pixel canvas.
    Save a perfectly 512x512 image with no white margins if save_path is given.
    """

    outlines = rescale_lines_to_image(outlines, (pixels, pixels),  offset=(x_start, y_start))  # if this did not work, use commented old function for rescaling
    blobs = rescale_lines_to_image(blobs, (pixels, pixels),  offset=(x_start, y_start))

    dpi = 100
    fig = plt.figure(figsize=(pixels / dpi, pixels / dpi), dpi=dpi)

    def process_layer(layer):
        unique_pts = list(set(layer))
        return sorted(unique_pts, key=lambda pt: pt[0])

    # --- Plot blobs (blue)
    for blob in blobs:
        cleaned = process_layer(blob)
        if cleaned:
            x = [p[0] for p in cleaned]
            y = [p[1] for p in cleaned]
            plt.plot(x, y, c=color2, linewidth=linewidth2, linestyle=linestyle2)

    # --- Plot outlines (green)
    for outline in outlines:
        cleaned = process_layer(outline)
        if cleaned:
            x = [p[0] for p in cleaned]
            y = [p[1] for p in cleaned]
            plt.plot(x, y, c=color1, linewidth=linewidth1, linestyle=linestyle1)

    # --- Remove margins & make frame pixel-perfect
    plt.axis('off')
    plt.axis('equal')
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.xlim(0, pixels)
    plt.ylim(0, pixels)
    plt.margins(0, 0)
    fig.subplots_adjust(0, 0, 1, 1)
    fig.patch.set_facecolor('white')

    # --- Save version (no white border)
    if save_path:
        fig.savefig(
            save_path,
            dpi=dpi,
            bbox_inches='tight',
            pad_inches=0,
            facecolor=fig.get_facecolor()
        )
        
        
        base = find_base(save_path)
        
        # fig.savefig(
            # base + ".svg",
            # format="svg",
            # dpi=dpi,
            # bbox_inches='tight',
            # pad_inches=0,
            # facecolor=fig.get_facecolor()
        # )
                
        # fig.savefig(
            # base + ".eps",
            # format="eps",
            # dpi=dpi,
            # bbox_inches='tight',
            # pad_inches=0,
            # facecolor=fig.get_facecolor()
        # )
                
        

    plt.show()
    plt.close(fig)
    


def cross_correlation(manual_interpretation, original_lines, x_threshold=10, y_threshold=10):
    """
    For each manual line, find the original_line with the highest R²
    computed on shared x-values, restricting comparisons to overlapping
    bounding boxes. Lines that have no valid match are skipped.

    Parameters
    ----------
    manual_interpretation : list of list of (x, y)
        Manually interpreted true lines.
    original_lines : list of list of (x, y)
        Original predicted or raw lines.
    x_threshold: threshold for x bounding box search
    y_threshold: threshold for y bounding box search

    Returns
    -------
    lateral_interpretation : list of list of (x, y)
        Correlated lines from original_lines (without None entries).
    """
    lateral_interpretation = []
    R2s = []

    for m_idx, man_line in enumerate(manual_interpretation):
        if not man_line:
            continue

        x_man = np.array([p[0] for p in man_line])
        y_man = np.array([p[1] for p in man_line])
        x_min, x_max = x_man.min(), x_man.max()
        y_min, y_max = y_man.min(), y_man.max()

        best_r2 = -np.inf
        best_line = None

        for o_idx, orig_line in enumerate(original_lines):
            if not orig_line:
                continue

            x_orig = np.array([p[0] for p in orig_line])
            y_orig = np.array([p[1] for p in orig_line])

            ox_min, ox_max = x_orig.min(), x_orig.max()
            oy_min, oy_max = y_orig.min(), y_orig.max()

#             # Skip lines with no bounding box overlap # obsoleted because it excludes features R2 recognizes [TG] for R2.
#             if (ox_max < x_min - x_threshold) or (ox_min > x_max + x_threshold) or (oy_max < y_min - y_threshold) or (oy_min > y_max + y_threshold):
#                 continue

            shared_x = np.intersect1d(x_man, x_orig)
            if shared_x.size < 2:
                continue  # Need at least 2 shared points

            sx = np.sort(shared_x)
            man_lookup = {x: y for x, y in man_line}
            orig_lookup = {x: y for x, y in orig_line}
            y_man_filtered = np.array([man_lookup[x] for x in sx])
            y_orig_filtered = np.array([orig_lookup[x] for x in sx])

            try:
                r2 = r2_score(y_man_filtered, y_orig_filtered)
            except Exception:
                continue

            if np.isnan(r2):
                continue

            if r2 > best_r2:
                best_r2 = r2
                best_line = orig_line

        if best_line is not None:
            lateral_interpretation.append(best_line)
            # store tuple (first_y_value_of_best_line, best_r2)
            R2s.append((best_line[0][1], best_r2))

    # Sort R² results by first point’s y-value
    if R2s:
        R2s_sorted = sorted(R2s, key=lambda x: x[0])
        sorted_y = [y for y, _ in R2s_sorted]
        sorted_r2 = [r for _, r in R2s_sorted]

        df = pd.DataFrame(
            [sorted_r2],
            index=["R²"],
            columns=[f"Layer y={y:.2f}" for y in sorted_y]
        )
        print("\nR² Summary Table (ordered by layer first-point y-value):")
        print(df.to_string(index=True))
    else:
        print("\nNo correlated layers found — R² list is empty.")

    return lateral_interpretation

    
    

def stratigraphic_confusion_matrix(
    manual_interpretation,
    lateral_interpretation,
    interval_labels=None
):
    """
    Build a stratigraphic confusion matrix between manual (true) and lateral (predicted)
    interpretations, using global top/bottom boundaries for correct pixel correlation.
    """

    # --- 1) Shared x-values ---
    all_manual_x = set(x for line in manual_interpretation for x, _ in line)
    all_lateral_x = set(x for line in lateral_interpretation for x, _ in line)
    shared_x_all = sorted(list(all_manual_x & all_lateral_x))
    if len(shared_x_all) == 0:
        raise ValueError("No shared x-values between manual and lateral interpretations.")

    # --- helper: filter lines ---
    def filter_lines(lines):
        out = []
        for ln in lines:
            new_ln = [(x, y) for x, y in ln if x in shared_x_all]
            if len(new_ln) >= 2:
                out.append(new_ln)
        return out

    manual = filter_lines(manual_interpretation)
    lateral = filter_lines(lateral_interpretation)

    # --- 2) Sort top-to-bottom by first-point y ---
    manual.sort(key=lambda ln: ln[0][1], reverse=True)
    lateral.sort(key=lambda ln: ln[0][1], reverse=True)

    # --- 3) Mutual filtering ---
    if len(manual) != len(lateral):
        def keep_matching(src, ref):
            kept = []
            for s in src:
                xs = set(x for x, _ in s)
                keep = any(len(xs & set(x for x, _ in r)) >= 2 for r in ref)
                if keep:
                    kept.append(s)
            return kept
        manual = keep_matching(manual, lateral)
        lateral = keep_matching(lateral, manual)

    if len(manual) != len(lateral):
        raise ValueError(
            f"Unequal number of lines after filtering: manual={len(manual)}, lateral={len(lateral)}."
        )

    n_lines = len(manual)
    if n_lines < 1:
        raise ValueError("Need at least one matched layer to define intervals.")

    # === NEW PART: Add global top and bottom boundary lines ===
    all_y = [y for ln in manual + lateral for _, y in ln]
    global_y_max = max(all_y)
    global_y_min = min(all_y)

    def make_horizontal_line(y_val):
        return [(x, y_val) for x in shared_x_all]

    # Add one above and one below
    manual_with_boundaries = [make_horizontal_line(global_y_max)] + manual + [make_horizontal_line(global_y_min)]
    lateral_with_boundaries = [make_horizontal_line(global_y_max)] + lateral + [make_horizontal_line(global_y_min)]

    n_intervals = len(manual_with_boundaries) - 1  # includes two outer intervals
    true_intervals = n_intervals - 2  # exclude topmost and bottommost

    # --- 5) Validate interval labels ---
    if interval_labels is None:
        labels = [f"I{i}" for i in range(true_intervals)]
    else:
        if len(interval_labels) != true_intervals:
            raise ValueError(
                f"interval_labels length {len(interval_labels)} does not match number of true intervals {true_intervals}."
            )
        labels = list(interval_labels)

    # --- 6) Build interval pixel sets (excluding outer boundary intervals) ---
    def build_interval_pixel_sets(lines):
        interval_pixel_sets = []
        for i in range(n_intervals):
            top = lines[i]
            bot = lines[i + 1]

            top_map = {x: y for x, y in top}
            bot_map = {x: y for x, y in bot}

            pixels = set()
            for x in shared_x_all:
                if (x in top_map) and (x in bot_map):
                    y_top = int(np.round(top_map[x]))
                    y_bot = int(np.round(bot_map[x]))
                    y1, y2 = sorted((y_top, y_bot))
                    if y2 - y1 > 1:
                        for y in range(y1 + 1, y2):
                            pixels.add((int(x), int(y)))
            interval_pixel_sets.append(pixels)
        return interval_pixel_sets

    manual_pixel_sets = build_interval_pixel_sets(manual_with_boundaries)[1:-1]  # exclude outer boundaries
    lateral_pixel_sets = build_interval_pixel_sets(lateral_with_boundaries)[1:-1]

    # --- 7) Build confusion matrix ---
    counts = np.zeros((true_intervals, true_intervals), dtype=int)
    for i in range(true_intervals):
        A = manual_pixel_sets[i]
        for j in range(true_intervals):
            B = lateral_pixel_sets[j]
            overlap = len(A & B) if A else 0
            counts[i, j] = overlap

    # --- 8) Normalize per row ---
    perc = np.full_like(counts, np.nan, dtype=float)
    for i in range(true_intervals):
        denom = counts[i, :].sum()
        if denom > 0:
            perc[i, :] = counts[i, :] / denom * 100

    # --- 9) Build DataFrames ---
    df_counts = pd.DataFrame(counts, index=labels, columns=labels)
    df_percent = pd.DataFrame(np.round(perc, 3), index=labels, columns=labels)

    return df_percent, df_counts
    
    
    
def find_new_point_groups(lines_before, lines_after):
    """
    Find groups of new points that appear in lines_after but not in lines_before.
    Groups are based on the line index they belong to, and within each line,
    only consecutive points (by x with no gaps > 1) are grouped together.

    Args:
        lines_before: list of list of (x, y)
        lines_after: list of list of (x, y)

    Returns:
        new_point_groups: list of list of (x, y)
    """
    # Flatten points from lines_before into a set for quick lookup
    before_points = {pt for line in lines_before for pt in line}

    # Step 1: Identify new points and map them to their line index
    mapped_new_points = []
    for idx, line in enumerate(lines_after):
        for pt in line:
            if pt not in before_points:
                mapped_new_points.append((idx, pt))

    # Step 2: Group new points by their line index
    index_groups = {}
    for idx, pt in mapped_new_points:
        index_groups.setdefault(idx, []).append(pt)

    # Step 3: Within each group, split by x gaps > 1
    new_point_groups = []
    for pts in index_groups.values():
        # Sort points by x value
        pts_sorted = sorted(pts, key=lambda p: p[0])

        # Split into subgroups by consecutive x difference
        subgroup = [pts_sorted[0]]
        i = 1
        while i < len(pts_sorted):
            if pts_sorted[i][0] - pts_sorted[i - 1][0] <= 1:
                subgroup.append(pts_sorted[i])
            else:
                # Gap detected → save subgroup if valid
                if len(subgroup) >= 2:
                    new_point_groups.append(subgroup)
                subgroup = [pts_sorted[i]]
            i += 1

        # Add last subgroup if valid
        if len(subgroup) >= 2:
            new_point_groups.append(subgroup)

    return new_point_groups
    