"""
OutburstMapper -- interactive mapping of cometary outbursts and
surface changes on a 3D shape model.

A PyVista + Qt desktop tool built to map outburst source locations
and surface-change regions of interest (ROIs) onto the shape model
of comet 67P/Churyumov-Gerasimenko, using Rosetta SPICE geometry
and OSIRIS/NAVCAM imagery. It works with any triangulated shape
model (.obj/.ply/.stl/.vtk).

See README.md for setup, the data layout, and a guide to every
feature.

Run
---
    python outburst_mapper.py               (auto-loads the default
                                             model if present next
                                             to this script)
    python outburst_mapper.py my_model.obj  (loads that model)
"""

import copy
import json
import os
import re
import sys
import tempfile
import time

import numpy as np
import pvl
import pyvista as pv
import spiceypy as spice
from scipy.spatial import cKDTree
from vtkmodules.util import numpy_support as ns
from vtkmodules.util.vtkConstants import VTK_UNSIGNED_CHAR
from vtkmodules.vtkCommonCore import vtkIdList, vtkPoints
from vtkmodules.vtkFiltersGeneral import vtkOBBTree
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
from vtkmodules.vtkRenderingCore import vtkActor, vtkCellPicker, vtkPolyDataMapper
from pyvistaqt import QtInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import (
    QFileDialog, QColorDialog, QListWidgetItem, QMessageBox,
)

# Set ROI_GUI_TRACE=1 in the environment to log the linked-zoom path, step by
# step, to stderr. Off by default; there to make a freeze reportable instead
# of a mystery -- the last line printed is the step it never came back from.
_TRACE = bool(os.environ.get("ROI_GUI_TRACE"))
_TRACE_T0 = time.time()


def _trace(message):
    if _TRACE:
        sys.stderr.write(f"[{time.time() - _TRACE_T0:8.3f}s] {message}\n")
        sys.stderr.flush()


LINE_WIDTH = 4
FILL_OPACITY = 0.25
PAINT_ALPHA = 180

# (label, view direction, view-up) for the axis-view buttons, matching
# ParaView's +X/-X/+Y/-Y/+Z/-Z toolbar convention.
AXIS_VIEWS = [
    ("+X", (1, 0, 0), (0, 0, 1)),
    ("-X", (-1, 0, 0), (0, 0, 1)),
    ("+Y", (0, 1, 0), (0, 0, 1)),
    ("-Y", (0, -1, 0), (0, 0, 1)),
    ("+Z", (0, 0, 1), (0, 1, 0)),
    ("-Z", (0, 0, -1), (0, 1, 0)),
]

# SPICE instrument NAIF ID + FOV frame name for each camera. NAVCAM had two
# redundant heads (A/B); Rosetta operations almost exclusively used CAM-A,
# so "NAVCAM" here means NAVCAM-A.
SPICE_CAMERAS = {
    "NAC": (-226111, "ROS_OSIRIS_NAC"),
    "WAC": (-226112, "ROS_OSIRIS_WAC"),
    "NAVCAM": (-226170, "ROS_NAVCAM-A"),
}
# WAC and NAVCAM (both heads) archive their raw PDS3 frames mirrored
# left/right relative to NAC -- a quirk of those cameras' optical paths.
_PDS3_MIRROR_NAIF_IDS = {-226112, -226170, -226180}
# The comet's CK-based body-fixed frame -- confirmed (not assumed) to match
# the DLR shape-model family this app loads: the bundled DSK from the same
# DLR photogrammetry pipeline declares its frame as exactly this ID.
SPICE_COMET_FRAME = "67P/C-G_CK"
SPICE_COMET_NAME = "CHURYUMOV-GERASIMENKO"
SPICE_SC_NAME = "ROSETTA"

# Auto-loaded on startup (if present next to this script) when no model
# path is given on the command line; "Load model…" still opens anything.
DEFAULT_MODEL = "cg-dlr_spg-shap7-v1.0_125Kfacets.obj"
# ======================================================================
# Qt helpers
# ======================================================================

class _NoScrollSpinBox(QtWidgets.QDoubleSpinBox):
    """A spin box that ignores the mouse wheel unless it has been clicked
    into first.

    Qt gives spin boxes Qt.WheelFocus by default, so wheeling over one on
    the way past it hands it the focus AND changes its value. Scrolling
    the control panel therefore used to silently edit whatever the cursor
    swept over -- brush radius, border width, lat/lon, the map's centre
    meridian -- which is easy to miss and annoying to undo.

    Note it *ignores* the wheel event rather than swallowing it: an
    ignored event propagates to the parent, so the scroll area still
    scrolls. Swallowing it (an event filter returning True) would stop
    the panel scrolling at all whenever the cursor was over a spin box."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)  # click/tab to focus, not wheel

    def wheelEvent(self, event):
        if self.hasFocus():
            super().wheelEvent(event)
        else:
            event.ignore()


class _ReorderableList(QtWidgets.QListWidget):
    """The ROI list, reorderable by dragging a row to a new position.

    QListWidget performs an internal move as insert-then-remove rather than
    a true row move, so the dragged row ends up as a NEW item and the model
    passes through an intermediate state where the ROI appears twice. The
    order is only trustworthy once dropEvent has finished, which is exactly
    when this reports it."""

    reordered = QtCore.pyqtSignal()

    def dropEvent(self, event):
        super().dropEvent(event)
        self.reordered.emit()


# ======================================================================
# Geometry helpers
# ======================================================================

def build_roi_geometry(mesh, roi):
    """Return (patch, border) polydata for a freeform (painted) ROI."""
    idx = np.asarray(roi["cells"], dtype=np.int64)
    patch = mesh.extract_cells(idx)
    border = patch.extract_feature_edges(
        boundary_edges=True, feature_edges=False,
        manifold_edges=False, non_manifold_edges=False)
    return patch, border


def hex_to_rgb(hexstr):
    h = hexstr.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def spice_camera_view(instid, frame, et):
    """Return (position_km, boresight_dir, up_dir, view_angle_deg) for a
    SPICE instrument at time `et`, expressed in the comet-fixed
    SPICE_COMET_FRAME (the frame the DLR shape models are delivered in).
    Raises spiceypy.utils.exceptions.SpiceyError if there's no coverage
    (outside the mission, or a gap such as the hibernation cruise)."""
    position, _ = spice.spkpos(
        SPICE_SC_NAME, et, SPICE_COMET_FRAME, "NONE", SPICE_COMET_NAME)
    position = np.asarray(position, dtype=float)

    _, _, boresight, _, _ = spice.getfov(instid, 4)
    boresight = np.asarray(boresight, dtype=float)
    boresight /= np.linalg.norm(boresight)
    ref_vector = np.asarray(spice.gdpool(f"INS{instid}_FOV_REF_VECTOR", 0, 3),
                            dtype=float)
    ref_vector /= np.linalg.norm(ref_vector)
    ref_angle = float(spice.gdpool(f"INS{instid}_FOV_REF_ANGLE", 0, 1)[0])
    cross_angle = float(spice.gdpool(f"INS{instid}_FOV_CROSS_ANGLE", 0, 1)[0])

    rot = np.array(spice.pxform(frame, SPICE_COMET_FRAME, et))
    boresight_c = rot @ boresight
    ref_c = rot @ ref_vector
    up = ref_c - np.dot(ref_c, boresight_c) * boresight_c
    up /= np.linalg.norm(up)

    view_angle = ref_angle + cross_angle  # ~ mean of the two full FOV angles
    return position, boresight_c, up, view_angle


def spice_camera_basis(instid, frame, et):
    """Return (position_km, forward, right, up, right_angle_deg, up_angle_deg)
    for a SPICE instrument at time `et`, in the comet-fixed SPICE_COMET_FRAME.
    `up` is the FOV's reference axis (FOV_REF_VECTOR) -- confirmed to be the
    correct "up" by spice_camera_view(), which uses that same axis directly
    as the pyvista camera's up and produces correctly-oriented views -- with
    half-angle FOV_REF_ANGLE; `right` = cross(forward, up), matching VTK's
    own screen-right convention exactly (verified against
    vtkCamera.GetModelViewTransformMatrix()), with half-angle
    FOV_CROSS_ANGLE."""
    position, _ = spice.spkpos(
        SPICE_SC_NAME, et, SPICE_COMET_FRAME, "NONE", SPICE_COMET_NAME)
    position = np.asarray(position, dtype=float)

    _, _, boresight, _, _ = spice.getfov(instid, 4)
    boresight = np.asarray(boresight, dtype=float)
    boresight /= np.linalg.norm(boresight)
    ref_vector = np.asarray(spice.gdpool(f"INS{instid}_FOV_REF_VECTOR", 0, 3),
                            dtype=float)
    ref_vector /= np.linalg.norm(ref_vector)
    ref_angle = float(spice.gdpool(f"INS{instid}_FOV_REF_ANGLE", 0, 1)[0])
    cross_angle = float(spice.gdpool(f"INS{instid}_FOV_CROSS_ANGLE", 0, 1)[0])

    rot = np.array(spice.pxform(frame, SPICE_COMET_FRAME, et))
    forward = rot @ boresight
    up = rot @ ref_vector
    up = up - np.dot(up, forward) * forward
    up /= np.linalg.norm(up)
    right = np.cross(forward, up)

    return position, forward, right, up, cross_angle, ref_angle


# ======================================================================
# ISIS cube (.cub) images -- minimal reader (no full ISIS install needed)
# and projective texturing onto the shape model.
# ======================================================================

# ISIS "Real" (float32) special-pixel bit patterns: any pixel whose raw
# uint32 bits fall at or above this threshold is NULL/saturated, not data.
_ISIS_REAL_SPECIAL_MIN = np.uint32(0xFF7FFFFA)

_ISIS_DTYPES = {
    "UnsignedByte": "u1", "SignedByte": "i1", "UnsignedWord": "u2",
    "SignedWord": "i2", "UnsignedInteger": "u4", "SignedInteger": "i4",
    "Real": "f4", "Double": "f8",
}


def read_isis_cube(path):
    """Read an ISIS cube's label and single-band pixel data (Tile or
    BandSequential organization) without requiring a full ISIS install.
    Returns (band as float32 ndarray [lines, samples], pvl label)."""
    label = pvl.load(path)
    core = label["IsisCube"]["Core"]
    dims = core["Dimensions"]
    samples, lines, bands = dims["Samples"], dims["Lines"], dims["Bands"]
    ptype = str(core["Pixels"]["Type"])
    endian = "<" if core["Pixels"]["ByteOrder"] == "Lsb" else ">"
    dtype = np.dtype(endian + _ISIS_DTYPES[ptype])
    start = core["StartByte"] - 1

    with open(path, "rb") as f:
        f.seek(start)
        if core["Format"] == "Tile":
            ts, tl = core["TileSamples"], core["TileLines"]
            img = np.empty((bands, lines, samples), dtype=dtype)
            for b in range(bands):
                for y0 in range(0, lines, tl):
                    for x0 in range(0, samples, ts):
                        tile = np.fromfile(f, dtype=dtype,
                                           count=ts * tl).reshape(tl, ts)
                        y1, x1 = min(y0 + tl, lines), min(x0 + ts, samples)
                        img[b, y0:y1, x0:x1] = tile[:y1 - y0, :x1 - x0]
        else:
            img = np.fromfile(f, dtype=dtype,
                              count=bands * lines * samples)
            img = img.reshape(bands, lines, samples)

    band0 = img[0].astype(np.float32)
    return band0, label


def isis_valid_mask(band, ptype):
    """Boolean mask of real (non-NULL/non-saturated) pixels."""
    if ptype == "Real":
        return band.view(np.uint32) < _ISIS_REAL_SPECIAL_MIN
    return np.ones(band.shape, dtype=bool)


def stretch_values_to_uint8(band, valid_mask, lo, hi):
    """Linear stretch to a displayable 8-bit image given explicit
    low/high data values; invalid pixels are forced to 0 (black)."""
    if hi <= lo:
        hi = lo + 1e-12
    out = np.clip((band - lo) / (hi - lo), 0.0, 1.0)
    out8 = (out * 255).astype(np.uint8)
    out8[~valid_mask] = 0
    return out8


def stretch_to_uint8(band, valid_mask, lo_pct, hi_pct):
    """Percentile contrast stretch to a displayable 8-bit image; invalid
    pixels are forced to 0 (black)."""
    valid = band[valid_mask]
    if valid.size == 0:
        lo, hi = 0.0, 1.0
    else:
        lo, hi = np.percentile(valid, [lo_pct, hi_pct])
    return stretch_values_to_uint8(band, valid_mask, lo, hi)


def _pds3_dtype(sample_type, sample_bits):
    st = sample_type.upper()
    endian = "<" if ("LSB" in st or st == "PC_REAL") else ">"
    if "REAL" in st:
        kind = "f"
    elif "UNSIGNED" in st:
        kind = "u"
    else:
        kind = "i"
    return np.dtype(f"{endian}{kind}{sample_bits // 8}")


def read_pds3_image(path):
    """Read a PDS3 (.IMG/.LBL) image -- attached label (data + label in one
    file, e.g. an OSIRIS NAC/WAC release) or detached label (a separate
    .LBL pointing at a headerless .IMG, e.g. NAVCAM) -- into a float32
    array plus its parsed pvl label. Only reads the primary ^IMAGE object;
    ignores any auxiliary objects (e.g. a WAC frame's SIGMA_MAP/QUALITY_MAP)."""
    label = pvl.load(path)
    if "IMAGE" not in label:
        raise ValueError(
            "Not a readable PDS3 image label (missing IMAGE object). "
            "If this is a detached-label product (e.g. NAVCAM), open the "
            ".LBL file instead of the .IMG.")
    img_ptr = label["^IMAGE"]
    if isinstance(img_ptr, (list, tuple)):
        # Detached label, fixed-length records: ("FILE.IMG", record_number).
        img_filename, img_record = img_ptr
        img_path = os.path.join(os.path.dirname(os.path.abspath(path)),
                                str(img_filename))
    elif isinstance(img_ptr, str):
        # Detached label, RECORD_TYPE = UNDEFINED (no record structure --
        # e.g. some NAVCAM releases): just a bare filename, data starts at
        # byte 0 of that file.
        img_path = os.path.join(os.path.dirname(os.path.abspath(path)),
                                img_ptr)
        img_record = 1
    else:
        # Attached label: record number into this same file.
        img_path = path
        img_record = img_ptr

    record_type = str(label.get("RECORD_TYPE", "")).upper()
    if record_type == "UNDEFINED":
        start = 0
    else:
        if "RECORD_BYTES" not in label:
            raise ValueError(
                "Not a readable PDS3 image label (missing RECORD_BYTES for "
                "a non-UNDEFINED RECORD_TYPE).")
        record_bytes = int(label["RECORD_BYTES"])
        start = (int(img_record) - 1) * record_bytes

    image_obj = label["IMAGE"]
    samples = int(image_obj["LINE_SAMPLES"])
    lines = int(image_obj["LINES"])
    dtype = _pds3_dtype(str(image_obj["SAMPLE_TYPE"]), int(image_obj["SAMPLE_BITS"]))

    with open(img_path, "rb") as f:
        f.seek(start)
        data = np.fromfile(f, dtype=dtype, count=lines * samples)
    band = data.reshape(lines, samples).astype(np.float32)
    # A raw-storage row-0-first array needs flipping to match this app's
    # "row 0 = +up" convention (the one the projection math was validated
    # against) exactly when LINE_DISPLAY_DIRECTION = "DOWN" -- verified
    # against a WAC frame's own embedded ground-truth surface intercept
    # (IMAGE_POI/SURF_INT_CART_COORD): before this flip the projected
    # look-direction was off by ~1.06 deg from ESA's own answer; after it,
    # ~0.06 deg (consistent with the label's rounding, not a real error).
    # NAVCAM ships LINE_DISPLAY_DIRECTION = "UP" (the opposite flag), so it's
    # left unflipped -- inferred from that flag, not independently verified
    # the way WAC was (no equivalent ground-truth field in its label).
    line_dir = str(image_obj.get("LINE_DISPLAY_DIRECTION", "DOWN")).upper()
    if line_dir == "DOWN":
        band = band[::-1, :]
    return band, label


def pds3_valid_mask(band):
    return np.isfinite(band)


def pds3_naif_id(label):
    """Map a PDS3 label's INSTRUMENT_ID (+ CHANNEL_ID for NAVCAM's two
    redundant heads) to the matching SPICE instrument NAIF ID."""
    inst = str(label["INSTRUMENT_ID"]).strip().upper()
    if inst == "OSINAC":
        return -226111
    if inst == "OSIWAC":
        return -226112
    if inst == "NAVCAM":
        channel = str(label.get("CHANNEL_ID", "CAM1")).strip().upper()
        return -226180 if channel == "CAM2" else -226170
    raise ValueError(f"Unrecognized PDS3 INSTRUMENT_ID: {inst!r}")


def load_image_any(path):
    """Load a .cub (ISIS) or .IMG/.LBL (PDS3) image, returning a uniform
    (band, valid_mask, naif_id, start_time) regardless of source format."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".cub":
        band, label = read_isis_cube(path)
        ptype = str(label["IsisCube"]["Core"]["Pixels"]["Type"])
        valid = isis_valid_mask(band, ptype)
        naif_id = int(label["IsisCube"]["Kernels"]["NaifFrameCode"])
        start_time = label["IsisCube"]["Instrument"]["StartTime"]
    elif ext in (".img", ".lbl"):
        band, label = read_pds3_image(path)
        valid = pds3_valid_mask(band)
        naif_id = pds3_naif_id(label)
        start_time = label["START_TIME"]
        if naif_id in _PDS3_MIRROR_NAIF_IDS:
            band = np.ascontiguousarray(band[:, ::-1])
            valid = np.ascontiguousarray(valid[:, ::-1])
    else:
        raise ValueError(f"Unsupported image file type: {ext}")
    return band, valid, naif_id, start_time


def project_cube_onto_mesh(mesh, cam_pos, forward, right, up,
                           half_ang_h, half_ang_v, obb_tree):
    """Return (u, v, visible) per mesh point for projecting an image taken
    by a pinhole camera (position/orientation/half-angles given) onto the
    mesh, with per-vertex occlusion testing (line-of-sight ray cast) so the
    texture doesn't wrap onto the far side or shadowed terrain."""
    pts = mesh.points
    rel = pts - cam_pos
    depth = rel @ forward
    x_cam = rel @ right
    y_cam = rel @ up
    with np.errstate(divide="ignore", invalid="ignore"):
        ndc_x = np.where(depth > 0, x_cam / depth / np.tan(half_ang_h), np.nan)
        ndc_y = np.where(depth > 0, y_cam / depth / np.tan(half_ang_v), np.nan)
    u = (ndc_x + 1) / 2
    v = 1 - (ndc_y + 1) / 2
    in_frame = (depth > 0) & (u >= 0) & (u <= 1) & (v >= 0) & (v <= 1)

    visible = np.zeros(len(pts), dtype=bool)
    idx = np.nonzero(in_frame)[0]
    pts_out = vtkPoints()
    cell_ids = vtkIdList()
    for i in idx:
        target = cam_pos + (pts[i] - cam_pos) * 0.999
        hit = obb_tree.IntersectWithLine(tuple(cam_pos), tuple(target),
                                         pts_out, cell_ids)
        visible[i] = (hit == 0)
    return u, v, visible


# ======================================================================
# Equirectangular maps -- a whole-body map image (region map, outburst
# location map, ...) draped over the model by lat/lon, as opposed to a
# single frame projected from where a camera actually was.
# ======================================================================

def equirect_texcoords(points, centre_lon=0.0, east_positive=True):
    """Per-point (u, v) texture coordinates for draping an equirectangular
    map over `points`, which must be given in the model's NATIVE
    (as-loaded) frame -- the same planetocentric, east-positive frame the
    lat/lon lookup box uses, so a map and that box always agree.

    The map is taken to span a full 360 deg of longitude and -90..+90 of
    latitude edge to edge, with +90 along its top row and `centre_lon`
    down its middle column. v=1 addresses the file's top row (checked
    against VTK's own texture orientation rather than assumed), so
    lat=+90 -> v=1."""
    r = np.linalg.norm(points, axis=1)
    r = np.where(r == 0, 1.0, r)
    lat = np.degrees(np.arcsin(np.clip(points[:, 2] / r, -1.0, 1.0)))
    lon = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    if not east_positive:
        lon = -lon
    u = ((lon - (centre_lon + 180.0)) % 360.0) / 360.0
    v = (lat + 90.0) / 180.0
    return u, v


def split_texture_seam(points, faces, u, v):
    """Give every cell straddling the map's left/right edge its own copies
    of its low-u corners, carrying u+1, so the texture runs forwards
    across the seam instead of snapping backwards across the entire map
    (which otherwise smears a band of garbage down the anti-meridian).
    Relies on the texture repeating for u>1, which is VTK's default.
    Returns (points, faces, u, v), each possibly lengthened.

    A corner is shifted when it lies more than half a map BEHIND its own
    cell's rightmost corner -- not merely at u<0.5, which would also drag
    a mid-map corner of a wide cell (the near-polar ones span a lot of
    longitude) a full map's width off target."""
    uf = u[faces]
    shift = (uf.max(axis=1)[:, None] - uf) > 0.5
    if not shift.any():
        return points, faces, u, v
    low = np.unique(faces[shift])
    remap = np.zeros(len(points), dtype=np.int64)
    remap[low] = np.arange(len(low), dtype=np.int64) + len(points)
    return (np.vstack([points, points[low]]),
            np.where(shift, remap[faces], faces),
            np.concatenate([u, u[low] + 1.0]),
            np.concatenate([v, v[low]]))


# ======================================================================
# Paint interaction — subclass vtkInteractorStyleTrackballCamera so ordinary navigation is untouched
# (we just call the inherited On*Event handlers) and only left-drag is
# intercepted while paint mode is on. No event-abort tricks needed.
# ======================================================================

class _PaintStyle(vtkInteractorStyleTrackballCamera):

    def __init__(self, manager):
        super().__init__()
        self.manager = manager
        self.left_down = False
        self.last_facet = None
        self.press_pos = None       # where a non-painting left-press started,
                                    # to tell a probe click from a camera drag
        self.middle_pos = None      # linked middle-drag (pan) in progress
        self.right_pos = None       # linked right-drag (zoom) in progress
        self.AddObserver("LeftButtonPressEvent", self._press)
        self.AddObserver("LeftButtonReleaseEvent", self._release)
        self.AddObserver("MiddleButtonPressEvent", self._middle_press)
        self.AddObserver("MiddleButtonReleaseEvent", self._middle_release)
        self.AddObserver("RightButtonPressEvent", self._right_press)
        self.AddObserver("RightButtonReleaseEvent", self._right_release)
        self.AddObserver("MouseMoveEvent", self._move)
        self.AddObserver("MouseWheelForwardEvent", self._wheel)
        self.AddObserver("MouseWheelBackwardEvent", self._wheel)

    def _wheel(self, _obj=None, event=""):
        """While the panels are linked, the wheel changes the visible angular
        window rather than dollying the camera -- dollying would walk the
        camera off the spacecraft's actual position, which is the one thing
        the link depends on."""
        m = self.manager
        forward = "Forward" in str(event)
        if getattr(m, "zoom_link", False):
            _trace(f"wheel: linked zoom {'in' if forward else 'out'} - start")
            m.zoom_sync_view(0.8 if forward else 1.25)
            _trace("wheel: linked zoom - done")
            return
        _trace("wheel: plain dolly")
        if forward:
            self.OnMouseWheelForward()
        else:
            self.OnMouseWheelBackward()

    # While the panels are linked, EVERY camera-moving interaction has to go
    # through the linked-window path (pan_sync_view / zoom_sync_view): the
    # inherited trackball handlers move the camera's actual pose, which walks
    # it off the spacecraft position -- the image stops following and the
    # link is silently broken. Left-drag and the wheel were already
    # intercepted; middle-drag (trackball pan) and right-drag (trackball
    # dolly) fell through to the inherited handlers and did exactly that.
    def _middle_press(self, *_):
        m = self.manager
        if getattr(m, "zoom_link", False):
            self.middle_pos = m.event_position()
            return
        self.OnMiddleButtonDown()

    def _middle_release(self, *_):
        if self.middle_pos is not None:
            self.middle_pos = None
            return
        self.OnMiddleButtonUp()

    def _right_press(self, *_):
        m = self.manager
        if getattr(m, "zoom_link", False):
            self.right_pos = m.event_position()
            return
        self.OnRightButtonDown()

    def _right_release(self, *_):
        if self.right_pos is not None:
            self.right_pos = None
            return
        self.OnRightButtonUp()

    def _press(self, *_):
        m = self.manager
        if not m.paint_mode or m.mesh is None:
            self.press_pos = m.event_position()
            self.OnLeftButtonDown()
            return
        self.left_down = True
        self.last_facet = None
        fid = m._pick_facet()
        if fid >= 0:
            m._stamp(m._brush_facets(fid))
            self.last_facet = fid

    def _release(self, *_):
        m = self.manager
        if not m.paint_mode:
            self.OnLeftButtonUp()
            # A click that didn't travel is a probe; anything that moved was
            # a camera drag and is left well alone.
            if m.probe_mode and self.press_pos is not None:
                x, y = m.event_position()
                if abs(x - self.press_pos[0]) <= 3 and abs(y - self.press_pos[1]) <= 3:
                    m.probe_facet()
            self.press_pos = None
            return
        self.left_down = False
        self.last_facet = None

    def _move(self, *_):
        m = self.manager
        if getattr(m, "zoom_link", False):
            if self.middle_pos is not None:
                x, y = m.event_position()
                m.pan_sync_view(x - self.middle_pos[0], y - self.middle_pos[1])
                self.middle_pos = (x, y)
                return
            if self.right_pos is not None:
                x, y = m.event_position()
                # up = in, matching the trackball dolly's direction
                m.zoom_sync_view(1.01 ** (self.right_pos[1] - y))
                self.right_pos = (x, y)
                return
        else:
            # the link died mid-drag: don't drive an unlinked camera
            self.middle_pos = self.right_pos = None
        if not m.paint_mode:
            if (getattr(m, "zoom_link", False) and self.press_pos is not None):
                x, y = m.event_position()
                m.pan_sync_view(x - self.press_pos[0], y - self.press_pos[1])
                self.press_pos = (x, y)
                return
            self.OnMouseMove()
            return
        if not self.left_down:
            return
        fid = m._pick_facet()
        if fid < 0:
            return
        facets = (m._between_facets(self.last_facet, fid)
                  if self.last_facet is not None else m._brush_facets(fid))
        m._stamp(facets)
        self.last_facet = fid


# ======================================================================
# Raw-image viewer for a loaded .cub, with zoom (mouse wheel) and an
# adjustable percentile stretch. "Project onto shape model" hands off to
# the main window, which owns the 3D scene.
# ======================================================================

class _ZoomableView(QtWidgets.QGraphicsView):
    """Wheel-zoom + pan by default; a "stretch tool" mode (toggled by the
    dialog) switches left-drag to rubber-band a region instead, emitting
    its bounds (in scene == image-pixel coordinates) on release -- the
    ISIS qview-style "drag to stretch to this region" interaction."""

    regionDragged = QtCore.pyqtSignal(QtCore.QRectF)
    viewChanged = QtCore.pyqtSignal()   # visible area moved or zoomed

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setRenderHints(QtGui.QPainter.Antialiasing
                            | QtGui.QPainter.SmoothPixmapTransform)
        self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        # No visible scrollbars, ever. Panning is hand-drag and framing is
        # programmatic, so they add nothing -- and under ScrollBarAsNeeded
        # they are actively destructive while the panels are linked: a
        # show_rect() that lands near the fits/doesn't-fit threshold makes a
        # scrollbar appear in a DEFERRED layout pass (after _quiet has been
        # dropped), which shrinks the viewport, which fires valueChanged,
        # which sync_3d_to_image reads as a panel resize and answers with
        # another show_rect() at the new viewport size -- hiding the
        # scrollbar again. The two panels then reframe each other forever
        # (one "image reframed" per cycle) and the event loop never idles.
        # Hidden scrollbars still scroll internally; only the flapping stops.
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.stretch_mode = False
        self.user_zoomed = False
        self._rubber_band = None
        self._origin = None
        self._quiet = False         # suppress viewChanged while being driven
        for bar in (self.horizontalScrollBar(), self.verticalScrollBar()):
            bar.valueChanged.connect(self._announce)   # panning

    def _announce(self, *_):
        if not self._quiet:
            self.viewChanged.emit()

    def visible_rect(self):
        """The image-pixel rectangle currently on screen."""
        return self.mapToScene(self.viewport().rect()).boundingRect()

    def show_rect(self, rect, quiet=True):
        """Frame this image-pixel rectangle, without reporting back (so the
        two panels can drive each other without ringing).

        Deliberately not fitInView(): that leaves a hard-coded two-pixel
        margin on every side, so a rectangle handed back and forth between
        two linked views would grow by about a percent each time and the
        framing would slowly zoom itself out. Setting the transform is
        exact."""
        was, self._quiet = self._quiet, quiet
        try:
            self.user_zoomed = True
            viewport = self.viewport().rect()
            if (rect.width() <= 0 or rect.height() <= 0
                    or viewport.width() <= 0 or viewport.height() <= 0):
                return
            scale = min(viewport.width() / rect.width(),
                        viewport.height() / rect.height())
            self.setTransform(QtGui.QTransform.fromScale(scale, scale))
            self.centerOn(rect.center())
        finally:
            self._quiet = was

    def set_roam(self, rect):
        """Give the view room to scroll well past the image while linked.

        QGraphicsView will not scroll past its sceneRect, and centerOn()
        silently pins to it rather than failing. So a linked framing that
        runs off the image's edge -- which the link's clamps deliberately
        allow, and which the whole-frame view (zero scroll range) hits
        immediately -- simply did not happen: the 3D camera panned and the
        image panel sat still, out of sync. Widening the VIEW's own
        sceneRect restores the room; the SCENE's sceneRect stays at the
        pixmap's bounds, so fit_to_content is unaffected. Pass None to
        snap back to tracking the scene's bounds when the link ends."""
        was, self._quiet = self._quiet, True   # a range change may scroll
        try:
            if rect is None:
                self.setSceneRect(QtCore.QRectF())  # track the scene again
            else:
                # the clamps allow the framing centre 1.5 frames out and
                # half-extents of 2 frames, so 3 frames of margin a side
                # covers every reachable framing
                mw, mh = 3.0 * rect.width(), 3.0 * rect.height()
                self.setSceneRect(rect.adjusted(-mw, -mh, mw, mh))
        finally:
            self._quiet = was

    def set_stretch_mode(self, enabled):
        self.stretch_mode = enabled
        self.setDragMode(QtWidgets.QGraphicsView.NoDrag if enabled
                         else QtWidgets.QGraphicsView.ScrollHandDrag)

    def fit_to_content(self):
        scene_items = self.scene().items()
        if scene_items:
            self.fitInView(self.scene().sceneRect(), QtCore.Qt.KeepAspectRatio)

    def wheelEvent(self, event):
        self.user_zoomed = True
        # Scale by how far the wheel actually turned. A mouse notch is 120;
        # a trackpad sends a burst of much smaller deltas, and treating each
        # as a full notch (the old fixed 1.25) meant one flick compounded
        # into an enormous zoom.
        steps = event.angleDelta().y() / 120.0
        factor = float(np.clip(1.25 ** steps, 0.2, 5.0))
        self.scale(factor, factor)
        self._announce()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if not self.user_zoomed:
            self.fit_to_content()
        self._announce()

    def mousePressEvent(self, event):
        if self.stretch_mode and event.button() == QtCore.Qt.LeftButton:
            self._origin = event.pos()
            if self._rubber_band is None:
                self._rubber_band = QtWidgets.QRubberBand(
                    QtWidgets.QRubberBand.Rectangle, self)
            self._rubber_band.setGeometry(QtCore.QRect(self._origin, QtCore.QSize()))
            self._rubber_band.show()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self.stretch_mode and self._origin is not None:
            rect = QtCore.QRect(self._origin, event.pos()).normalized()
            self._rubber_band.setGeometry(rect)
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self.stretch_mode and self._origin is not None:
            rect = QtCore.QRect(self._origin, event.pos()).normalized()
            self._rubber_band.hide()
            self._origin = None
            scene_rect = self.mapToScene(rect).boundingRect()
            if scene_rect.width() > 1 and scene_rect.height() > 1:
                self.regionDragged.emit(scene_rect)
            return
        super().mouseReleaseEvent(event)


class ImagePanel(QtCore.QObject):
    """The loaded image, shown beside the 3D view.

    Deliberately NOT one widget containing both its controls and its display:
    the controls live in a toolbar that spans the full width of BOTH views,
    leaving the image's display area exactly the same size as the 3D view's.
    Stacking the controls above the image instead would make its display area
    shorter than the 3D view's, so the two would have different aspect ratios
    and a linked zoom could not frame the same ground in both."""

    def __init__(self, manager):
        super().__init__(manager)
        self.manager = manager
        self.band = None
        self.valid_mask = None
        self.naif_id = None
        self.start_time = None
        self.title = ""
        self._img8 = None
        self._explicit_range = None  # (lo, hi) from drag-stretch, else None

        # --- the toolbar, spanning both views ---
        self.toolbar = QtWidgets.QWidget()
        bar = QtWidgets.QHBoxLayout(self.toolbar)
        bar.setContentsMargins(4, 2, 4, 2)

        self.lbl_title = QtWidgets.QLabel("No image loaded")
        bar.addWidget(self.lbl_title)
        bar.addSpacing(12)

        bar.addWidget(QtWidgets.QLabel("Stretch %ile"))
        self.sp_lo = _NoScrollSpinBox()
        self.sp_lo.setRange(0.0, 100.0); self.sp_lo.setDecimals(2)
        self.sp_lo.setValue(1.0); self.sp_lo.setMaximumWidth(80)
        bar.addWidget(self.sp_lo)
        self.sp_hi = _NoScrollSpinBox()
        self.sp_hi.setRange(0.0, 100.0); self.sp_hi.setDecimals(2)
        self.sp_hi.setValue(99.5); self.sp_hi.setMaximumWidth(80)
        bar.addWidget(self.sp_hi)
        btn_stretch = QtWidgets.QPushButton("Apply")
        btn_stretch.clicked.connect(self._apply_percentile_stretch)
        bar.addWidget(btn_stretch)

        self.btn_stretch_tool = QtWidgets.QPushButton("Stretch tool (drag)")
        self.btn_stretch_tool.setCheckable(True)
        self.btn_stretch_tool.setToolTip(
            "Off: drag pans, wheel zooms.  On: drag a box to stretch to that "
            "region's min/max (like ISIS qview).")
        self.btn_stretch_tool.toggled.connect(self._toggle_stretch_tool)
        bar.addWidget(self.btn_stretch_tool)
        bar.addSpacing(12)

        btn_view_sc = QtWidgets.QPushButton("View from spacecraft")
        btn_view_sc.setToolTip(
            "Put the 3D camera where this image was taken from, at its own "
            "acquisition time, and link the two panels.")
        btn_view_sc.clicked.connect(self.view_from_spacecraft)
        bar.addWidget(btn_view_sc)
        btn_project = QtWidgets.QPushButton("Project onto model")
        btn_project.clicked.connect(self._project)
        bar.addWidget(btn_project)
        btn_remove = QtWidgets.QPushButton("Remove projection")
        btn_remove.clicked.connect(self.manager.remove_cube_projection)
        bar.addWidget(btn_remove)
        bar.addSpacing(12)

        self.chk_link = QtWidgets.QCheckBox("Link zoom")
        self.chk_link.setEnabled(False)
        self.chk_link.toggled.connect(self.manager.set_zoom_link)
        bar.addWidget(self.chk_link)
        btn_fit = QtWidgets.QPushButton("Fit")
        btn_fit.clicked.connect(self.fit)
        bar.addWidget(btn_fit)
        bar.addStretch(1)
        btn_close = QtWidgets.QPushButton("Close image")
        btn_close.setToolTip(
            "Hide the image panel and go back to the single 3D view. "
            "Unlinks the panels (navigation returns to normal); a "
            "projection already on the model stays until "
            "\"Remove projection\".")
        btn_close.clicked.connect(self.manager.close_image_panel)
        bar.addWidget(btn_close)

        # --- the display, sized to match the 3D view exactly ---
        self.scene = QtWidgets.QGraphicsScene(self)
        self.view = _ZoomableView()
        self.view.setScene(self.scene)
        self.view.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                QtWidgets.QSizePolicy.Expanding)
        # no frame: a QGraphicsView's default border would make its viewport
        # two pixels smaller than the 3D view beside it, and the two panels
        # would then differ slightly in aspect
        self.view.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.view.regionDragged.connect(self._apply_drag_stretch)
        self.view.viewChanged.connect(self.manager.sync_3d_to_image)
        self.pixmap_item = None

    def set_visible(self, visible):
        self.toolbar.setVisible(visible)
        self.view.setVisible(visible)

    def clear(self):
        """Forget the loaded image (the panel is about to be hidden)."""
        self.band = None
        self.valid_mask = None
        self.naif_id = None
        self.start_time = None
        self.title = ""
        self._img8 = None
        self._explicit_range = None
        self.scene.clear()
        self.pixmap_item = None
        self.view.user_zoomed = False
        self.lbl_title.setText("No image loaded")
        self.chk_link.blockSignals(True)
        self.chk_link.setChecked(False)
        self.chk_link.blockSignals(False)
        self.chk_link.setEnabled(False)

    def set_image(self, band, valid_mask, naif_id, start_time, title):
        self.band = band
        self.valid_mask = valid_mask
        self.naif_id = naif_id
        self.start_time = start_time
        self.title = title
        self._explicit_range = None
        self.lbl_title.setText(f"{title} — {self.utc()} UTC")
        self.chk_link.setEnabled(True)
        self._update_pixmap()
        self.fit()

    def utc(self):
        """The image's own acquisition time, as an ISO string SPICE parses."""
        if self.start_time is None:
            return ""
        if hasattr(self.start_time, "strftime"):
            return self.start_time.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]
        return str(self.start_time)

    def fit(self):
        self.view.user_zoomed = False
        self.view.fit_to_content()

    def view_from_spacecraft(self):
        """Jump the 3D camera to this image's own viewpoint and moment. No
        projection needed first -- the instrument and time come off the
        image's own label, which is all SPICE needs."""
        if self.naif_id is None:
            return
        self.manager.view_from_image(self.naif_id, self.start_time)

    def _toggle_stretch_tool(self, checked):
        self.view.set_stretch_mode(checked)

    def _apply_percentile_stretch(self):
        self._explicit_range = None
        self._update_pixmap()

    def _apply_drag_stretch(self, scene_rect):
        x0 = max(0, int(scene_rect.left()))
        y0 = max(0, int(scene_rect.top()))
        x1 = min(self.band.shape[1], int(np.ceil(scene_rect.right())))
        y1 = min(self.band.shape[0], int(np.ceil(scene_rect.bottom())))
        if x1 <= x0 or y1 <= y0:
            return
        region = self.band[y0:y1, x0:x1]
        region_valid = self.valid_mask[y0:y1, x0:x1]
        valid_vals = region[region_valid]
        if valid_vals.size == 0:
            return
        self._explicit_range = (float(valid_vals.min()), float(valid_vals.max()))
        self._update_pixmap()

    def _update_pixmap(self):
        if self.band is None:
            return
        if self._explicit_range is not None:
            lo, hi = self._explicit_range
            img8 = stretch_values_to_uint8(self.band, self.valid_mask, lo, hi)
        else:
            img8 = stretch_to_uint8(self.band, self.valid_mask,
                                    self.sp_lo.value(), self.sp_hi.value())
        self._img8 = np.ascontiguousarray(img8)
        h, w = self._img8.shape
        qimg = QtGui.QImage(self._img8.data, w, h, w,
                            QtGui.QImage.Format_Grayscale8).copy()
        pixmap = QtGui.QPixmap.fromImage(qimg)
        # keep where the user is looking: restretching shouldn't re-frame
        keep = self.view.visible_rect() if self.view.user_zoomed else None
        self.scene.clear()
        self.pixmap_item = self.scene.addPixmap(pixmap)
        self.scene.setSceneRect(QtCore.QRectF(pixmap.rect()))
        if keep is not None and keep.width() > 1:
            self.view.show_rect(keep)
        else:
            self.view.fit_to_content()

    def _project(self):
        if self._img8 is not None:
            self.manager.project_cube(self._img8, self.naif_id, self.start_time)


# ======================================================================
# Main window
# ======================================================================

class ROIManager(QtWidgets.QMainWindow):

    def __init__(self, model_path=None):
        super().__init__()
        self.setWindowTitle("OutburstMapper")
        self._size_to_screen()

        self.mesh = None            # displayed (possibly aligned) mesh
        self.mesh_normals = None    # same mesh with point normals
        self._mesh_original_points = None  # as-loaded points, pre-alignment
        self._latlon_marker_actor = None
        self._axes_actors = {}
        self.model_path = None
        self.transform = None       # 4x4 alignment matrix or None
        self.rois = []              # list of dicts
        self.actors = {}            # roi id -> {"border": actor, "fill": actor}
        self._undo_stack = []       # snapshots of self.rois, newest last,
                                    # taken just before each saved change
        self._counter = 0
        self._updating_ui = False
        self._model_actor = None

        # facet probe state
        self.probe_mode = False
        self._probe_actor = None

        # paint tool state
        self.paint_mode = False
        self._paint_cells = set()   # cell ids in the in-progress ROI
        self._face_centers = None
        self._face_tree = None
        self._picker = vtkCellPicker()
        self._picker.SetTolerance(0.0005)
        self._picker.PickFromListOn()
        self._paint_overlay_actor = None
        self._paint_overlay_array = None
        self._paint_rgba = None
        self._style = None
        self._editing_roi_id = None  # id of the ROI being edited in place,
                                     # or None for "painting a new one"

        self._spice_loaded = False

        self._cube_projection_actor = None
        self._cube_cam_params = None

        # image panel <-> 3D camera link
        self._sync_geom = None      # image/camera geometry from SPICE
        self.zoom_link = False
        self._syncing = False       # guards the two panels driving each other
        self._last_image_viewport = None
        self._render_pending = False

        self._map_actor = None
        self._map_texture = None
        self._map_path = None

        self._session_path = None  # set on manual save/load; auto-save
                                    # then keeps writing to that same file

        self._build_ui()
        self._auto_load_spice()
        if model_path:
            self.load_model(model_path)

    def _size_to_screen(self):
        """Size/center the window to the screen it's launched on, instead
        of a fixed pixel size that may be too big (small laptop) or too
        small (large external monitor)."""
        screen = QtWidgets.QApplication.primaryScreen()
        if screen is None:
            self.resize(1400, 900)
            return
        avail = screen.availableGeometry()
        w = int(avail.width() * 0.9)
        h = int(avail.height() * 0.75)
        self.resize(w, h)
        self.move(avail.x() + (avail.width() - w) // 2,
                 avail.y() + (avail.height() - h) // 2)

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        h = QtWidgets.QHBoxLayout(central)

        # --- left control panel (scrollable -- there are a lot of controls,
        # and the window itself is deliberately not full screen height) ---
        panel_widget = QtWidgets.QWidget()
        panel = QtWidgets.QVBoxLayout(panel_widget)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(panel_widget)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        h.addWidget(scroll, 0)

        btn_load = QtWidgets.QPushButton("Load model…")
        btn_load.clicked.connect(self.dialog_load_model)
        panel.addWidget(btn_load)

        panel.addWidget(QtWidgets.QLabel("View along axis:"))
        view_grid = QtWidgets.QGridLayout()
        for i, (label, vec, up) in enumerate(AXIS_VIEWS):
            b = QtWidgets.QPushButton(label)
            b.clicked.connect(lambda _=False, v=vec, u=up: self._set_view(v, u))
            view_grid.addWidget(b, i // 3, i % 3)
        panel.addLayout(view_grid)

        panel.addWidget(self._hline())
        panel.addWidget(QtWidgets.QLabel("SPICE camera view (Rosetta):"))
        btn_spice_load = QtWidgets.QPushButton("Load SPICE kernels…")
        btn_spice_load.clicked.connect(lambda: self.load_spice_kernels(None))
        panel.addWidget(btn_spice_load)
        self.lbl_spice_status = QtWidgets.QLabel("SPICE: not loaded")
        panel.addWidget(self.lbl_spice_status)

        panel.addWidget(QtWidgets.QLabel("Date/time (UTC) — paste or type:"))
        self.ed_spice_time = QtWidgets.QLineEdit("2014-08-06T10:00:00")
        self.ed_spice_time.setPlaceholderText(
            "e.g. 2015-08-23T02:08:19.862")
        panel.addWidget(self.ed_spice_time)

        cam_row = QtWidgets.QHBoxLayout()
        self._spice_cam_group = QtWidgets.QButtonGroup(self)
        self._spice_cam_group.setExclusive(True)
        for name in ("NAC", "WAC", "NAVCAM"):
            b = QtWidgets.QPushButton(name)
            b.setCheckable(True)
            b.clicked.connect(lambda _=False, n=name: self._apply_camera_view(n))
            self._spice_cam_group.addButton(b)
            cam_row.addWidget(b)
        panel.addLayout(cam_row)
        self.lbl_spice_result = QtWidgets.QLabel("")
        self.lbl_spice_result.setWordWrap(True)
        panel.addWidget(self.lbl_spice_result)

        btn_load_cube = QtWidgets.QPushButton("Load image (.cub/.IMG/.LBL)…")
        btn_load_cube.clicked.connect(self.load_cube_image)
        panel.addWidget(btn_load_cube)

        self.chk_show_projection = QtWidgets.QCheckBox("Show projected image")
        self.chk_show_projection.setEnabled(False)
        self.chk_show_projection.stateChanged.connect(self._toggle_cube_projection)
        panel.addWidget(self.chk_show_projection)

        panel.addWidget(self._hline())
        panel.addWidget(QtWidgets.QLabel("Equirectangular map:"))
        btn_load_map = QtWidgets.QPushButton("Load map…")
        btn_load_map.clicked.connect(self.load_map)
        panel.addWidget(btn_load_map)

        map_row = QtWidgets.QHBoxLayout()
        self.chk_show_map = QtWidgets.QCheckBox("Show map")
        self.chk_show_map.setEnabled(False)
        self.chk_show_map.stateChanged.connect(self._toggle_map)
        map_row.addWidget(self.chk_show_map)
        btn_remove_map = QtWidgets.QPushButton("Remove map")
        btn_remove_map.clicked.connect(self.remove_map)
        map_row.addWidget(btn_remove_map)
        panel.addLayout(map_row)

        map_row2 = QtWidgets.QHBoxLayout()
        map_row2.addWidget(QtWidgets.QLabel("Lon at centre"))
        self.sp_map_lon = _NoScrollSpinBox()
        self.sp_map_lon.setRange(-360.0, 360.0)
        self.sp_map_lon.setDecimals(1); self.sp_map_lon.setSingleStep(10.0)
        self.sp_map_lon.setValue(0.0)
        self.sp_map_lon.valueChanged.connect(self._reproject_map)
        map_row2.addWidget(self.sp_map_lon)
        self.chk_map_west = QtWidgets.QCheckBox("West-positive")
        self.chk_map_west.stateChanged.connect(self._reproject_map)
        map_row2.addWidget(self.chk_map_west)
        panel.addLayout(map_row2)
        self.lbl_map = QtWidgets.QLabel("No map loaded")
        self.lbl_map.setWordWrap(True)
        panel.addWidget(self.lbl_map)

        panel.addWidget(self._hline())
        self.chk_edges = QtWidgets.QCheckBox("Show plate edges")
        self.chk_edges.setChecked(False)
        self.chk_edges.stateChanged.connect(self._restyle_model)
        panel.addWidget(self.chk_edges)

        self.chk_axes = QtWidgets.QCheckBox("Show axes (R=X, G=Y, B=Z)")
        self.chk_axes.setChecked(True)
        self.chk_axes.stateChanged.connect(self._toggle_axes)
        panel.addWidget(self.chk_axes)

        panel.addWidget(self._hline())
        panel.addWidget(QtWidgets.QLabel("Find a point by lat/lon:"))
        latlon_row = QtWidgets.QHBoxLayout()
        latlon_row.addWidget(QtWidgets.QLabel("Lat"))
        self.sp_lat = _NoScrollSpinBox()
        self.sp_lat.setRange(-90.0, 90.0)
        self.sp_lat.setDecimals(2); self.sp_lat.setSingleStep(1.0)
        latlon_row.addWidget(self.sp_lat)
        latlon_row.addWidget(QtWidgets.QLabel("Lon"))
        self.sp_lon = _NoScrollSpinBox()
        self.sp_lon.setRange(0.0, 360.0)
        self.sp_lon.setDecimals(2); self.sp_lon.setSingleStep(1.0)
        latlon_row.addWidget(self.sp_lon)
        btn_latlon = QtWidgets.QPushButton("Go")
        btn_latlon.clicked.connect(self.goto_latlon)
        latlon_row.addWidget(btn_latlon)
        panel.addLayout(latlon_row)
        self.chk_marker = QtWidgets.QCheckBox("Show marker")
        self.chk_marker.setChecked(True)
        self.chk_marker.stateChanged.connect(self._toggle_marker_visibility)
        panel.addWidget(self.chk_marker)
        panel.addWidget(QtWidgets.QLabel(
            "Planetocentric, east-positive (0-360°), in the\n"
            "model's native (as-loaded) frame — unaffected by alignment."))

        panel.addWidget(self._hline())
        panel.addWidget(QtWidgets.QLabel("Probe a facet:"))
        probe_row = QtWidgets.QHBoxLayout()
        self.btn_probe = QtWidgets.QPushButton("Probe mode")
        self.btn_probe.setCheckable(True)
        self.btn_probe.toggled.connect(self._toggle_probe_mode)
        probe_row.addWidget(self.btn_probe)
        btn_probe_clear = QtWidgets.QPushButton("Clear normal")
        btn_probe_clear.clicked.connect(self.clear_probe)
        probe_row.addWidget(btn_probe_clear)
        panel.addLayout(probe_row)
        self.lbl_probe = QtWidgets.QLabel("No facet probed yet.")
        self.lbl_probe.setWordWrap(True)
        panel.addWidget(self.lbl_probe)
        panel.addWidget(QtWidgets.QLabel(
            "Click a facet for its lat/lon/radius (same frame as above)\n"
            "and a long line out along its normal. Dragging still rotates\n"
            "the view, and the lat/lon boxes follow each click."))

        panel.addWidget(self._hline())
        panel.addWidget(QtWidgets.QLabel("Paint an ROI (freeform):"))
        paint_row = QtWidgets.QHBoxLayout()
        self.btn_paint = QtWidgets.QPushButton("Paint mode")
        self.btn_paint.setCheckable(True)
        self.btn_paint.toggled.connect(self._toggle_paint_mode)
        paint_row.addWidget(self.btn_paint)
        self.chk_erase = QtWidgets.QCheckBox("Erase")
        self.chk_erase.setEnabled(False)
        paint_row.addWidget(self.chk_erase)
        panel.addLayout(paint_row)

        paint_row2 = QtWidgets.QHBoxLayout()
        paint_row2.addWidget(QtWidgets.QLabel("Brush r"))
        self.sp_brush = self._spin(0.1)
        paint_row2.addWidget(self.sp_brush)
        self.btn_paint_clear = QtWidgets.QPushButton("Clear paint")
        self.btn_paint_clear.setEnabled(False)
        self.btn_paint_clear.clicked.connect(self._clear_paint)
        paint_row2.addWidget(self.btn_paint_clear)
        panel.addLayout(paint_row2)
        panel.addWidget(QtWidgets.QLabel(
            "Drag to paint cells (multiple strokes OK).\n"
            "Toggle Paint mode off to save as a new ROI."))

        edit_row = QtWidgets.QHBoxLayout()
        self.btn_edit = QtWidgets.QPushButton("Edit selected ROI")
        self.btn_edit.clicked.connect(self.edit_selected)
        edit_row.addWidget(self.btn_edit)
        self.btn_cancel_edit = QtWidgets.QPushButton("Cancel edit")
        self.btn_cancel_edit.setEnabled(False)
        self.btn_cancel_edit.clicked.connect(self.cancel_edit)
        edit_row.addWidget(self.btn_cancel_edit)
        panel.addLayout(edit_row)
        panel.addWidget(QtWidgets.QLabel(
            "Reopens an existing ROI (from this or a loaded session) in\n"
            "the paint tool, its cells already painted, so you can add to\n"
            "or erase parts of it. Toggle Paint mode off to save the edit\n"
            "in place — the ROI keeps its id and position in the list."))

        # ROI parameter controls
        form = QtWidgets.QFormLayout()
        self.ed_name = QtWidgets.QLineEdit()
        form.addRow("Name", self.ed_name)

        self.sp_lw = _NoScrollSpinBox()
        self.sp_lw.setDecimals(1); self.sp_lw.setRange(0.5, 20.0)
        self.sp_lw.setSingleStep(0.5); self.sp_lw.setValue(LINE_WIDTH)
        form.addRow("Border width", self.sp_lw)
        panel.addLayout(form)

        row = QtWidgets.QHBoxLayout()
        self.btn_color = QtWidgets.QPushButton("Color…")
        self.btn_color.clicked.connect(self.choose_color)
        self._color = "#ff2020"
        self._style_color_button()
        row.addWidget(self.btn_color)
        btn_apply = QtWidgets.QPushButton("Apply to selected")
        btn_apply.clicked.connect(self.apply_to_selected)
        row.addWidget(btn_apply)
        panel.addLayout(row)

        self.chk_fill = QtWidgets.QCheckBox("Show translucent fills")
        self.chk_fill.stateChanged.connect(self.refresh_visibility)
        panel.addWidget(self.chk_fill)

        panel.addWidget(self._hline())
        panel.addWidget(QtWidgets.QLabel("ROIs (checkbox = visible):"))
        self.listw = _ReorderableList()
        self.listw.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
        self.listw.setDefaultDropAction(QtCore.Qt.MoveAction)
        self.listw.itemChanged.connect(self.refresh_visibility)
        self.listw.currentItemChanged.connect(self.roi_selected)
        self.listw.reordered.connect(self._sync_roi_order)
        panel.addWidget(self.listw, 1)

        move_row = QtWidgets.QHBoxLayout()
        btn_up = QtWidgets.QPushButton("Move up")
        btn_up.clicked.connect(lambda: self.move_selected(-1))
        move_row.addWidget(btn_up)
        btn_down = QtWidgets.QPushButton("Move down")
        btn_down.clicked.connect(lambda: self.move_selected(1))
        move_row.addWidget(btn_down)
        panel.addLayout(move_row)
        panel.addWidget(QtWidgets.QLabel(
            "Drag a row to move it anywhere in the list; the buttons nudge\n"
            "it one place. The order is saved with the session."))

        vis_row = QtWidgets.QHBoxLayout()
        btn_show_all = QtWidgets.QPushButton("Select all")
        btn_show_all.clicked.connect(lambda: self._set_all_visible(True))
        vis_row.addWidget(btn_show_all)
        btn_hide_all = QtWidgets.QPushButton("Deselect all")
        btn_hide_all.clicked.connect(lambda: self._set_all_visible(False))
        vis_row.addWidget(btn_hide_all)
        panel.addLayout(vis_row)

        del_row = QtWidgets.QHBoxLayout()
        btn_del = QtWidgets.QPushButton("Delete selected ROI")
        btn_del.clicked.connect(self.delete_selected)
        del_row.addWidget(btn_del)
        self.btn_undo = QtWidgets.QPushButton("Undo")
        self.btn_undo.setEnabled(False)
        self.btn_undo.setToolTip(
            "Undo the last ROI change (new ROI, finished edit, apply, or "
            "delete). Restores the list and auto-saves the restored state.")
        self.btn_undo.clicked.connect(self.undo)
        del_row.addWidget(self.btn_undo)
        panel.addLayout(del_row)

        panel.addWidget(self._hline())
        row2 = QtWidgets.QHBoxLayout()
        b_save = QtWidgets.QPushButton("Save session…")
        b_save.clicked.connect(self.save_session)
        b_open = QtWidgets.QPushButton("Load session…")
        b_open.clicked.connect(self.load_session)
        row2.addWidget(b_save); row2.addWidget(b_open)
        panel.addLayout(row2)
        self.lbl_autosave = QtWidgets.QLabel(
            "New ROIs auto-save the session (to ROI_data/ by default).")
        self.lbl_autosave.setWordWrap(True)
        panel.addWidget(self.lbl_autosave)

        b_shot = QtWidgets.QPushButton("Save screenshot…")
        b_shot.clicked.connect(self.screenshot)
        panel.addWidget(b_shot)

        # --- 3D view, with the image beside it and one toolbar over both ---
        viewing = QtWidgets.QWidget()
        stack = QtWidgets.QVBoxLayout(viewing)
        stack.setContentsMargins(0, 0, 0, 0)
        stack.setSpacing(2)
        self.image_panel = ImagePanel(self)
        stack.addWidget(self.image_panel.toolbar)

        self.views_row = QtWidgets.QHBoxLayout()
        self.views_row.setContentsMargins(0, 0, 0, 0)
        self.views_row.setSpacing(2)
        self.plotter = QtInteractor(viewing)
        self.plotter.set_background((0.12, 0.12, 0.14))
        self.plotter.interactor.setSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        # equal stretch and no splitter between them: the two display areas
        # must stay exactly the same size for a linked zoom to frame the same
        # ground in both
        self.views_row.addWidget(self.plotter.interactor, 1)
        self.views_row.addWidget(self.image_panel.view, 1)
        stack.addLayout(self.views_row, 1)
        h.addWidget(viewing, 1)
        self.image_panel.set_visible(False)   # appears when an image loads
        self._style = _PaintStyle(self)
        self._style.SetCurrentRenderer(self.plotter.renderer)
        # Set it through pyvista's own `style` property (not just the raw
        # interactor's SetInteractorStyle) so pyvista's internal bookkeeping
        # (self.plotter.iren._style_class) actually points at our style.
        # Otherwise pyvista's own double-click detector -- always active,
        # completely unrelated to anything we do -- calls its internal
        # update_style() on the next quick double-click-shaped input (e.g.
        # two paint strokes started close together in time and position)
        # and silently reverts the interactor to its stale default style,
        # after which our paint handlers never fire again.
        self.plotter.iren.style = self._style

    def _spin(self, val):
        s = _NoScrollSpinBox()
        s.setDecimals(4); s.setRange(1e-4, 1e6)
        s.setSingleStep(0.01); s.setValue(val)
        return s

    @staticmethod
    def _hline():
        f = QtWidgets.QFrame(); f.setFrameShape(QtWidgets.QFrame.HLine)
        return f

    def _style_color_button(self):
        self.btn_color.setStyleSheet(f"background-color: {self._color};")

    # ------------------------------------------------------------------
    # Model handling
    # ------------------------------------------------------------------
    def dialog_load_model(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load shape model", "",
            "Meshes (*.obj *.ply *.stl *.vtk *.vtp);;All files (*)")
        if path:
            self.load_model(path)

    def load_model(self, path, transform=None):
        self._discard_paint_session()
        self.model_path = path
        self.mesh = pv.read(path).extract_surface().triangulate()
        self._mesh_original_points = self.mesh.points.copy()
        self.transform = None
        if transform is not None:
            self.mesh.transform(np.asarray(transform), inplace=True)
            self.transform = np.asarray(transform)
        self.mesh_normals = self.mesh.compute_normals(
            point_normals=True, cell_normals=False, auto_orient_normals=True)
        self._clear_rois()
        # old snapshots index cells of the previous model, so they can't be
        # restored onto this one
        self._undo_stack = []
        self.btn_undo.setEnabled(False)
        self._session_path = None  # a new model starts a fresh session file
        self.lbl_autosave.setText(
            "New ROIs auto-save the session (to ROI_data/ by default).")
        self.plotter.clear()
        self._latlon_marker_actor = None
        self._probe_actor = None
        self.lbl_probe.setText("No facet probed yet.")
        self._cube_projection_actor = None
        self._cube_cam_params = None
        self._sync_geom = None      # the link is tied to a model + a moment
        self.set_zoom_link(False)
        # a map is body-specific, so it doesn't follow to a new model
        self._map_actor = None
        self._map_texture = None
        self._map_path = None
        self.chk_show_map.blockSignals(True)
        self.chk_show_map.setChecked(False)
        self.chk_show_map.blockSignals(False)
        self.chk_show_map.setEnabled(False)
        self.lbl_map.setText("No map loaded")
        self.chk_show_projection.blockSignals(True)
        self.chk_show_projection.setChecked(False)
        self.chk_show_projection.blockSignals(False)
        self.chk_show_projection.setEnabled(False)
        # Renderer.clear() strips ALL lights and never re-adds any, which
        # otherwise leaves the model rendering as a flat, unlit silhouette.
        self.plotter.enable_lightkit()
        self._display_mesh()
        self._rebuild_axes()
        self._rebuild_face_tree()
        self._ensure_paint_overlay()
        self.plotter.reset_camera()

    def _rebuild_axes(self):
        """(Re)draw world-frame X/Y/Z axis lines through the origin,
        colored red/green/blue, sized to the current mesh's extent."""
        for name in ("axis_x", "axis_y", "axis_z"):
            self.plotter.remove_actor(name, render=False)
        self._axes_actors = {}
        if self.mesh is None:
            return
        xmin, xmax, ymin, ymax, zmin, zmax = self.mesh.bounds
        margin = 1.15
        ex = max(abs(xmin), abs(xmax)) * margin
        ey = max(abs(ymin), abs(ymax)) * margin
        ez = max(abs(zmin), abs(zmax)) * margin
        specs = [
            ("axis_x", (-ex, 0, 0), (ex, 0, 0), "red"),
            ("axis_y", (0, -ey, 0), (0, ey, 0), "green"),
            ("axis_z", (0, 0, -ez), (0, 0, ez), "blue"),
        ]
        show = self.chk_axes.isChecked()
        for name, p0, p1, color in specs:
            line = pv.Line(p0, p1)
            actor = self.plotter.add_mesh(
                line, color=color, line_width=3, pickable=False, name=name)
            actor.SetVisibility(show)
            self._axes_actors[name] = actor

    def _toggle_axes(self, *_):
        show = self.chk_axes.isChecked()
        for actor in self._axes_actors.values():
            actor.SetVisibility(show)
        self.plotter.render()

    def _display_mesh(self):
        """(Re)add the model actor as a plain ParaView-style "Surface"
        view: Gouraud-shaded, no extra shininess, edges optional."""
        if self.mesh is None:
            return
        self._model_actor = self.plotter.add_mesh(
            self.mesh, color="lightgray",
            show_edges=self.chk_edges.isChecked(),
            edge_color="dimgray", line_width=1,
            name="model")
        prop = self._model_actor.GetProperty()
        prop.SetInterpolationToGouraud()
        prop.SetAmbient(0.0)
        prop.SetDiffuse(1.0)
        prop.SetSpecular(0.0)
        self._picker.InitializePickList()
        self._picker.AddPickList(self._model_actor)

    def _restyle_model(self, *_):
        self._display_mesh()
        self.plotter.render()

    def _set_view(self, vector, viewup):
        if self.mesh is None:
            return
        self.plotter.view_vector(vector, viewup=viewup)
        self.plotter.reset_camera()
        self.plotter.render()

    def goto_latlon(self):
        """Highlight the surface point nearest the given planetocentric
        lat/lon (defined in the as-loaded frame, so it stays meaningful
        regardless of any later principal-axis alignment) and point the
        camera straight down at it along the local surface normal.
        Painting is still done by hand from there."""
        if self.mesh is None or self._mesh_original_points is None:
            return
        lat = np.radians(self.sp_lat.value())
        lon = np.radians(self.sp_lon.value())
        target = np.array([np.cos(lat) * np.cos(lon),
                           np.cos(lat) * np.sin(lon),
                           np.sin(lat)])
        dirs = self._mesh_original_points
        norms = np.linalg.norm(dirs, axis=1)
        norms[norms == 0] = 1.0
        dots = (dirs / norms[:, None]) @ target
        idx = int(np.argmax(dots))
        point = self.mesh.points[idx]

        # Fixed size relative to the mesh's own scale -- NOT the paint-brush
        # radius (self.sp_brush), which the user changes freely for painting
        # and previously made this marker's size swing wildly along with it.
        radius = max(self.mesh.length * 0.006, 1e-3)
        marker = pv.Sphere(radius=radius, center=point,
                           theta_resolution=16, phi_resolution=16)
        self._latlon_marker_actor = self.plotter.add_mesh(
            marker, color="yellow", pickable=False, name="latlon_marker")
        self._latlon_marker_actor.SetVisibility(self.chk_marker.isChecked())

        # "view from above" = camera pulled back along the local surface
        # normal at this vertex (mesh_normals shares point indexing with
        # mesh since compute_normals() is called with split_vertices=False)
        normal = self.mesh_normals["Normals"][idx]
        normal = normal / np.linalg.norm(normal)
        up_ref = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(up_ref, normal)) > 0.9:
            up_ref = np.array([0.0, 1.0, 0.0])
        up = up_ref - np.dot(up_ref, normal) * normal
        up /= np.linalg.norm(up)

        cam = self.plotter.camera
        distance = cam.distance
        cam.focal_point = tuple(point)
        cam.position = tuple(point + normal * distance)
        cam.up = tuple(up)
        self.plotter.reset_camera_clipping_range()
        self.plotter.render()

    def _toggle_marker_visibility(self, *_):
        if self._latlon_marker_actor is not None:
            self._latlon_marker_actor.SetVisibility(self.chk_marker.isChecked())
            self.plotter.render()

    # ------------------------------------------------------------------
    # Facet probe — click a facet for its lat/lon/radius and its normal
    # ------------------------------------------------------------------
    def _toggle_probe_mode(self, checked):
        if checked and self.paint_mode:
            QMessageBox.information(
                self, "Paint mode is on",
                "Painting already uses the left mouse button. Toggle Paint "
                "mode off first, then turn Probe mode on.")
            self.btn_probe.blockSignals(True)
            self.btn_probe.setChecked(False)
            self.btn_probe.blockSignals(False)
            return
        self.probe_mode = checked

    def event_position(self):
        return self.plotter.iren.interactor.GetEventPosition()

    def probe_facet(self):
        """Report the clicked facet's lat/lon/radius, and draw a long line
        out along its outward normal (like the X/Y/Z axis lines)."""
        if self.mesh is None:
            return
        fid = self._pick_facet()
        if fid < 0:
            return
        ids = self.mesh.faces.reshape(-1, 4)[fid, 1:4]

        # lat/lon/radius come from the AS-LOADED points, so they agree with
        # the lat/lon box and the map drape whatever alignment is applied.
        native = (self._mesh_original_points
                  if self._mesh_original_points is not None else self.mesh.points)
        p = native[ids].mean(axis=0)
        r = float(np.linalg.norm(p))
        lat = np.degrees(np.arcsin(np.clip(p[2] / r, -1.0, 1.0))) if r else 0.0
        lon = float(np.degrees(np.arctan2(p[1], p[0])) % 360.0)

        # the line itself is drawn on the DISPLAYED mesh. Averaging the three
        # corner normals rather than taking the cross product of the winding
        # gets the outward side for free -- load_model computed those with
        # auto_orient_normals.
        centre = self.mesh.points[ids].mean(axis=0)
        normal = self.mesh_normals["Normals"][ids].mean(axis=0)
        normal = normal / np.linalg.norm(normal)
        length = self.mesh.length
        line = pv.Line(centre - normal * 0.05 * length,
                       centre + normal * 1.20 * length)
        self._probe_actor = self.plotter.add_mesh(
            line, color="magenta", line_width=3, pickable=False,
            name="probe_normal")

        self.lbl_probe.setText(
            f"Facet {fid}:  lat {lat:.2f}°,  lon {lon:.2f}°,  r {r:.3f} km\n"
            f"normal (model frame): "
            f"[{normal[0]:.3f}, {normal[1]:.3f}, {normal[2]:.3f}]")
        # so "Go" re-centres on what was just clicked, and a probed point can
        # be handed straight to the lat/lon workflow
        self.sp_lat.setValue(lat)
        self.sp_lon.setValue(lon)
        self.plotter.render()

    def clear_probe(self):
        self.plotter.remove_actor("probe_normal", render=False)
        self._probe_actor = None
        self.lbl_probe.setText("No facet probed yet.")
        self.plotter.render()

    # ------------------------------------------------------------------
    # SPICE camera view
    # ------------------------------------------------------------------
    def _auto_load_spice(self):
        here = os.path.dirname(os.path.abspath(__file__))
        candidate = os.path.join(here, "ROSETTA", "kernels", "mk", "ROS_OPS.TM")
        if os.path.isfile(candidate):
            self.load_spice_kernels(candidate)

    def load_spice_kernels(self, path):
        if not path:
            path, _ = QFileDialog.getOpenFileName(
                self, "Load SPICE metakernel", "",
                "Metakernels (*.tm *.TM);;All files (*)")
            if not path:
                return
        try:
            mk_dir = os.path.dirname(os.path.abspath(path))
            kernels_root = os.path.dirname(mk_dir)
            with open(path) as f:
                text = f.read()
            # PATH_VALUES in the metakernel is resolved relative to the
            # CALLER's working directory (not the metakernel's own
            # location), so rewrite it to an absolute path -- otherwise
            # loading only works if the app happens to be launched with
            # the mk/ folder as the current directory.
            text = re.sub(r"PATH_VALUES\s*=\s*\(\s*'[^']*'\s*\)",
                          f"PATH_VALUES = ( '{kernels_root}' )", text)
            fd, tmp_path = tempfile.mkstemp(suffix=".tm")
            try:
                with os.fdopen(fd, "w") as f:
                    f.write(text)
                spice.kclear()
                spice.furnsh(tmp_path)
            finally:
                os.remove(tmp_path)
        except Exception as e:
            self._spice_loaded = False
            self.lbl_spice_status.setText("SPICE: not loaded")
            QMessageBox.warning(self, "SPICE error",
                                f"Could not load metakernel:\n{path}\n\n{e}")
            return
        self._spice_loaded = True
        self.lbl_spice_status.setText(f"SPICE: loaded ({os.path.basename(path)})")

    def _apply_camera_view(self, name):
        if self.mesh is None:
            return
        if not self._spice_loaded:
            QMessageBox.warning(
                self, "SPICE not loaded",
                "Load the Rosetta SPICE kernels first "
                "(\"Load SPICE kernels…\" button).")
            return
        utc = self.ed_spice_time.text().strip()
        try:
            et = spice.str2et(utc)
        except Exception as e:
            QMessageBox.warning(
                self, "Unrecognized date/time",
                f"Couldn't parse \"{utc}\" as a UTC date/time.\n\n"
                "Try an ISO 8601-style string, e.g.\n"
                "2015-08-23T02:08:19.862\n\n"
                f"Details: {e}")
            return
        try:
            instid, frame = SPICE_CAMERAS[name]
            pos, boresight, up, view_angle = spice_camera_view(instid, frame, et)
        except Exception as e:
            QMessageBox.warning(
                self, "No SPICE data for that time",
                f"Couldn't compute the {name} view at {utc} UTC.\n\n"
                "This time is likely outside the mission or in a coverage "
                "gap (e.g. the 2011-2014 hibernation cruise).\n\n"
                f"Details: {e}")
            return

        dist = float(np.linalg.norm(pos))
        cam = self.plotter.camera
        cam.position = tuple(pos)
        cam.focal_point = tuple(pos + boresight * dist)
        cam.up = tuple(up)
        cam.view_angle = view_angle
        self.plotter.reset_camera_clipping_range()
        self.plotter.render()
        self.lbl_spice_result.setText(
            f"{name} @ {utc} UTC — {dist:.2f} km from comet center")

    # ------------------------------------------------------------------
    # ISIS cube (.cub) image loading + projection onto the shape model
    # ------------------------------------------------------------------
    def load_cube_image(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load image", "",
            "Images (*.cub *.IMG *.LBL);;ISIS cubes (*.cub);;"
            "PDS3 images (*.IMG *.LBL);;All files (*)")
        if not path:
            return
        try:
            band, valid, naif_id, start_time = load_image_any(path)
        except Exception as e:
            QMessageBox.warning(self, "Could not read image",
                                f"{path}\n\n{e}")
            return
        self.set_zoom_link(False)
        self.image_panel.set_image(band, valid, naif_id, start_time,
                                   os.path.basename(path))
        self.image_panel.set_visible(True)
        # so the SPICE camera-view buttons default to this image's own time
        self.ed_spice_time.setText(self.image_panel.utc())

    def close_image_panel(self):
        """Back to the plain single-panel view, as before an image was
        loaded. Unlinks first (which restores normal navigation); a
        projection already draped on the model is left alone -- "Remove
        projection" exists for that."""
        self.set_zoom_link(False)
        self._sync_geom = None
        self._last_image_viewport = None
        self.image_panel.set_visible(False)
        self.image_panel.clear()

    def project_cube(self, img8, naif_id, start_time):
        if self.mesh is None:
            return
        if not self._spice_loaded:
            QMessageBox.warning(
                self, "SPICE not loaded",
                "Load the Rosetta SPICE kernels first "
                "(\"Load SPICE kernels…\" button).")
            return
        try:
            utc = start_time.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]
            et = spice.str2et(utc)
            _, frame, _, _, _ = spice.getfov(naif_id, 4)
            pos, fwd, right, up, right_ang, up_ang = spice_camera_basis(
                naif_id, frame, et)
        except Exception as e:
            QMessageBox.warning(
                self, "No SPICE data for this image",
                f"Couldn't compute camera geometry for this image "
                f"(instrument/time from its label).\n\n{e}")
            return

        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            tree = vtkOBBTree()
            tree.SetDataSet(self.mesh)
            tree.BuildLocator()
            u, v, visible = project_cube_onto_mesh(
                self.mesh, pos, fwd, right, up,
                np.radians(right_ang), np.radians(up_ang), tree)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()

        faces = self.mesh.faces.reshape(-1, 4)[:, 1:4]
        cell_ok = visible[faces].all(axis=1)
        if not cell_ok.any():
            QMessageBox.information(
                self, "No overlap",
                "This image's footprint doesn't fall on the currently "
                "loaded model at its acquisition time.")
            return

        patch = self.mesh.extract_cells(np.nonzero(cell_ok)[0])
        orig_ids = patch["vtkOriginalPointIds"]
        # numpy_to_texture maps texcoord v=1 to array row 0, so flip v here
        # (project_cube_onto_mesh's v follows normal image-row convention:
        # v=0 at the top/+up row).
        patch.active_texture_coordinates = np.column_stack(
            [u[orig_ids], 1 - v[orig_ids]])

        texture = pv.numpy_to_texture(np.dstack([img8, img8, img8]))
        self.plotter.remove_actor("cube_projection", render=False)
        self._cube_projection_actor = self.plotter.add_mesh(
            patch, texture=texture, pickable=False, name="cube_projection")
        self.chk_show_projection.setEnabled(True)
        self.chk_show_projection.blockSignals(True)
        self.chk_show_projection.setChecked(True)
        self.chk_show_projection.blockSignals(False)

        view_angle = right_ang + up_ang
        self._cube_cam_params = (pos, fwd, up, view_angle)
        # spacecraft view is the default on project, and links the panels
        self.view_from_image(naif_id, start_time)

    def _toggle_cube_projection(self, *_):
        # Just a visibility flip on the already-computed actor/texture, so
        # it's instantaneous no matter how long the projection itself took.
        if self._cube_projection_actor is not None:
            self._cube_projection_actor.SetVisibility(
                self.chk_show_projection.isChecked())
            self.plotter.render()

    def view_from_image(self, naif_id, start_time):
        """Put the 3D camera exactly where the loaded image was taken from,
        at the image's own acquisition time.

        Everything needed comes from the image's own label -- instrument and
        time -- so this no longer waits on a projection first; projecting is
        a much heavier, and quite separate, thing to want."""
        if self.mesh is None:
            return
        if not self._spice_loaded:
            QMessageBox.warning(
                self, "SPICE not loaded",
                "Load the Rosetta SPICE kernels first "
                "(\"Load SPICE kernels…\" button).")
            return
        utc = self.image_panel.utc()
        try:
            et = spice.str2et(utc)
            _, frame, _, _, _ = spice.getfov(naif_id, 4)
            pos, fwd, right, up, right_ang, up_ang = spice_camera_basis(
                naif_id, frame, et)
        except Exception as e:
            QMessageBox.warning(
                self, "No SPICE data for this image",
                f"Couldn't work out where {SPICE_SC_NAME} was at {utc} UTC.\n\n{e}")
            return

        # the geometry that ties image pixels to camera angles, kept so the
        # two panels can be driven from one another
        band = self.image_panel.band
        self._sync_geom = {
            "pos": pos, "forward": fwd, "right": right, "up": up,
            "tan_h": float(np.tan(np.radians(right_ang))),
            "tan_v": float(np.tan(np.radians(up_ang))),
            "width": int(band.shape[1]), "height": int(band.shape[0]),
        }
        self.ed_spice_time.setText(utc)
        self.lbl_spice_result.setText(
            f"{self.image_panel.title} @ {utc} UTC — "
            f"{float(np.linalg.norm(pos)):.2f} km from comet center")
        self._place_sync_camera()
        self.image_panel.fit()
        self.set_zoom_link(True)
        self.image_panel.chk_link.setChecked(True)

    def view_from_spacecraft(self):
        """The projected image's viewpoint (kept for the projection flow)."""
        if self._cube_cam_params is None:
            if self.image_panel.naif_id is not None:
                self.image_panel.view_from_spacecraft()
            return
        pos, fwd, up, view_angle = self._cube_cam_params
        dist = float(np.linalg.norm(pos))
        cam = self.plotter.camera
        cam.position = tuple(pos)
        cam.focal_point = tuple(pos + fwd * dist)
        cam.up = tuple(up)
        cam.view_angle = view_angle
        self.plotter.reset_camera_clipping_range()
        self.plotter.render()

    # ---- linked zoom between the image panel and the 3D view -----------
    def _place_sync_camera(self):
        """Stand the camera at the spacecraft, boresight down the image."""
        g = self._sync_geom
        cam = self.plotter.camera
        dist = float(np.linalg.norm(g["pos"]))
        cam.position = tuple(g["pos"])
        cam.focal_point = tuple(g["pos"] + g["forward"] * dist)
        cam.up = tuple(g["up"])
        # through the same path as every other framing, so the whole frame is
        # the one starting point and the limits apply uniformly
        self._set_camera_link_state(0.0, 0.0, g["tan_v"])

    def set_zoom_link(self, on):
        was = getattr(self, "zoom_link", False)
        self.zoom_link = bool(on) and self._sync_geom is not None
        if hasattr(self, "image_panel"):
            self.image_panel.chk_link.blockSignals(True)
            self.image_panel.chk_link.setChecked(self.zoom_link)
            self.image_panel.chk_link.blockSignals(False)
        if self.zoom_link:
            g = self._sync_geom
            self.image_panel.view.set_roam(
                QtCore.QRectF(0, 0, g["width"], g["height"]))
            self.sync_image_to_3d()
        else:
            if hasattr(self, "image_panel"):
                self.image_panel.view.set_roam(None)
            if was:
                self._restore_free_camera()

    def _restore_free_camera(self):
        """Fold the link's off-axis framing back into an ordinary camera.

        While linked, all navigation happens through view_angle and
        WindowCenter with the camera pose frozen at the spacecraft.
        Unlinking used to leave both behind: a non-zero WindowCenter keeps
        shearing the view, and the focal point -- the trackball's centre
        for orbit and pan -- is wherever the boresight happened to reach,
        kilometres from the comet, which makes ordinary navigation feel
        broken. So: re-aim the view direction at what is currently centred
        on screen (the WindowCenter pan folded into a real rotation, so
        nothing visibly jumps), zero the WindowCenter, and put the focal
        point at the comet centre's depth along that view ray."""
        if self.mesh is None:
            return
        cam = self.plotter.camera
        pos = np.asarray(cam.position, dtype=float)
        forward = np.asarray(cam.focal_point, dtype=float) - pos
        forward /= np.linalg.norm(forward)
        up = np.asarray(cam.up, dtype=float)
        up = up - np.dot(up, forward) * forward
        up /= np.linalg.norm(up)
        right = np.cross(forward, up)
        wcx, wcy = cam.GetWindowCenter()
        half_v = np.tan(np.radians(cam.view_angle) / 2.0)
        half_h = half_v * self._viewport_aspect()
        forward = forward + right * (wcx * half_h) + up * (wcy * half_v)
        forward /= np.linalg.norm(forward)
        depth = float(np.dot(np.asarray(self.mesh.center) - pos, forward))
        depth = max(depth, 0.05 * self.mesh.length)
        new_up = up - np.dot(up, forward) * forward
        cam.SetWindowCenter(0.0, 0.0)
        cam.focal_point = tuple(pos + forward * depth)
        cam.up = tuple(new_up / np.linalg.norm(new_up))
        self.plotter.reset_camera_clipping_range()
        self._request_render()

    def _viewport_aspect(self):
        w, h = self.plotter.renderer.GetSize()
        return (w / h) if h else 1.0

    def _camera_tangent_rect(self):
        """What the 3D camera sees, in tangent-of-angle off the boresight."""
        cam = self.plotter.camera
        half_v = np.tan(np.radians(cam.view_angle) / 2.0)
        half_h = half_v * self._viewport_aspect()
        cx, cy = cam.GetWindowCenter()
        return (cx * half_h - half_h, cx * half_h + half_h,
                cy * half_v - half_v, cy * half_v + half_v)

    # How far the linked zoom may go, in image pixels across the view and in
    # whole frames. Without limits every wheel event multiplies the last, so
    # one trackpad flick (a burst of dozens of events) drove the view angle to
    # 1e-6 degrees and the image's scale factor to 1e5 -- at which point the
    # scroll range approached the 2^31 integer limit and repainting a
    # hugely-magnified pixmap froze the app for seconds.
    LINK_MIN_PIXELS = 8.0
    LINK_MAX_FRAMES = 2.0

    def _clamp_link_state(self, cx, cy, half_v):
        g = self._sync_geom
        if not g:
            return cx, cy, max(float(half_v), 1e-9)
        low = g["tan_v"] * self.LINK_MIN_PIXELS / g["height"]
        high = g["tan_v"] * self.LINK_MAX_FRAMES
        half_v = float(np.clip(half_v, low, high))
        # and don't let it wander off the frame entirely
        return (float(np.clip(cx, -2.0 * g["tan_h"], 2.0 * g["tan_h"])),
                float(np.clip(cy, -2.0 * g["tan_v"], 2.0 * g["tan_v"])),
                half_v)

    def _set_camera_link_state(self, cx, cy, half_v):
        """Frame the view at this centre and half-height, in tangent space.
        WindowCenter shifts the frustum off-axis, so an off-centre crop is
        framed exactly rather than by pointing at its middle and hoping."""
        cx, cy, half_v = self._clamp_link_state(cx, cy, half_v)
        aspect = self._viewport_aspect()
        cam = self.plotter.camera
        cam.view_angle = float(np.degrees(2.0 * np.arctan(half_v)))
        cam.SetWindowCenter(cx / (half_v * aspect), cy / half_v)
        self.plotter.reset_camera_clipping_range()
        self._request_render()

    def _request_render(self):
        """Ask for a redraw without rendering here and now.

        A 3D render costs 10-25 ms and blocks, and a wheel gesture delivers a
        burst of events before the loop gets a chance to repaint, so rendering
        inside each handler froze the window for seconds.

        This asks Qt to repaint the widget instead of calling render()
        directly. Qt already collapses repeated update() calls into a single
        paintEvent, and the render then happens inside Qt's own paint cycle
        with the GL context current -- which on macOS is the safe place for
        it. Driving VTK's render from a timer callback outside that cycle is
        the classic way to get a wedged window on this platform."""
        _trace("    request render")
        widget = getattr(self.plotter, "interactor", None)
        if widget is not None:
            widget.update()
        else:
            self.plotter.render()

    def _flush_render(self):     # kept for tests that count redraw requests
        self._render_pending = False
        self._request_render()

    @staticmethod
    def _link_state_of(x0, x1, y0, y1):
        return (x0 + x1) / 2.0, (y0 + y1) / 2.0, abs(y1 - y0) / 2.0

    def _tangent_to_pixels(self, x0, x1, y0, y1):
        g = self._sync_geom
        px0 = (x0 / g["tan_h"] + 1.0) / 2.0 * g["width"]
        px1 = (x1 / g["tan_h"] + 1.0) / 2.0 * g["width"]
        # image rows run downwards, +up runs up, so the top row is max y
        py0 = (1.0 - y1 / g["tan_v"]) / 2.0 * g["height"]
        py1 = (1.0 - y0 / g["tan_v"]) / 2.0 * g["height"]
        return QtCore.QRectF(px0, py0, px1 - px0, py1 - py0)

    def _pixels_to_tangent(self, rect):
        g = self._sync_geom
        x0 = (2.0 * rect.left() / g["width"] - 1.0) * g["tan_h"]
        x1 = (2.0 * rect.right() / g["width"] - 1.0) * g["tan_h"]
        y1 = (1.0 - 2.0 * rect.top() / g["height"]) * g["tan_v"]
        y0 = (1.0 - 2.0 * rect.bottom() / g["height"]) * g["tan_v"]
        return x0, x1, y0, y1

    def _show_image_link_state(self, cx, cy, half_v):
        """Frame the image panel at the same centre and half-height, built at
        the PANEL's own aspect ratio so fitInView has nothing left to expand.

        Only the centre and height travel between the panels. The two have
        different aspect ratios, so their widths cannot both match; handing
        over each panel's fitted rectangle verbatim would let the other widen
        it to suit itself, and the framing would creep outwards a little on
        every hand-off until the view had zoomed itself out."""
        g = self._sync_geom
        view = self.image_panel.view
        vp = view.viewport().rect()
        aspect = (vp.width() / vp.height()) if vp.height() else 1.0
        half_py = half_v / g["tan_v"] * g["height"] / 2.0
        half_px = half_py * aspect
        centre_px = (cx / g["tan_h"] + 1.0) / 2.0 * g["width"]
        centre_py = (1.0 - cy / g["tan_v"]) / 2.0 * g["height"]
        view.show_rect(QtCore.QRectF(centre_px - half_px, centre_py - half_py,
                                     2 * half_px, 2 * half_py))

    def sync_3d_to_image(self):
        """The image panel moved: match the 3D camera to it."""
        if not getattr(self, "zoom_link", False) or self._syncing:
            return
        size = self.image_panel.view.viewport().size()
        if size != self._last_image_viewport:
            # The panel changed SIZE rather than the user re-framing it: it
            # now shows a bit more or less at the same scale. Resizing the
            # window shouldn't drag the camera around, so the camera stays
            # authoritative and the image is re-framed from it instead.
            self._last_image_viewport = size
            self.sync_image_to_3d()
            return
        self._syncing = True
        try:
            self._set_camera_link_state(*self._link_state_of(
                *self._pixels_to_tangent(self.image_panel.view.visible_rect())))
            # The camera clamps how far the link will zoom; mirror the clamped
            # result straight back, or the image panel could keep zooming on
            # its own past where the 3D view is willing to follow.
            self._show_image_link_state(
                *self._link_state_of(*self._camera_tangent_rect()))
            self._last_image_viewport = self.image_panel.view.viewport().size()
        finally:
            self._syncing = False

    def sync_image_to_3d(self):
        """The 3D camera moved: match the image panel to it."""
        if not getattr(self, "zoom_link", False) or self._syncing:
            return
        self._syncing = True
        try:
            self._show_image_link_state(
                *self._link_state_of(*self._camera_tangent_rect()))
            self._last_image_viewport = self.image_panel.view.viewport().size()
            _trace("    image reframed")
        finally:
            self._syncing = False

    def zoom_sync_view(self, factor):
        """Wheel zoom in the 3D panel while linked: scale the visible angular
        window about its centre, then pull the image panel along."""
        cx, cy, half_v = self._link_state_of(*self._camera_tangent_rect())
        _trace(f"  zoom_sync_view: half_v {half_v:.6g} -> {half_v*factor:.6g}")
        self._set_camera_link_state(cx, cy, half_v * factor)
        _trace("  camera set")
        self.sync_image_to_3d()
        _trace("  image panel synced")

    def pan_sync_view(self, dx, dy):
        """Left-drag in the 3D panel while linked: slide the visible window
        instead of orbiting, since orbiting would leave the spacecraft's
        viewpoint and there would be nothing left to link to."""
        w, h = self.plotter.renderer.GetSize()
        x0, x1, y0, y1 = self._camera_tangent_rect()
        cx, cy, half_v = self._link_state_of(x0, x1, y0, y1)
        self._set_camera_link_state(cx - (x1 - x0) * dx / max(w, 1),
                                    cy - (y1 - y0) * dy / max(h, 1), half_v)
        self.sync_image_to_3d()

    def remove_cube_projection(self):
        self.plotter.remove_actor("cube_projection", render=False)
        self._cube_projection_actor = None
        self._cube_cam_params = None
        self.chk_show_projection.blockSignals(True)
        self.chk_show_projection.setChecked(False)
        self.chk_show_projection.blockSignals(False)
        self.chk_show_projection.setEnabled(False)
        self.plotter.render()

    # ------------------------------------------------------------------
    # Equirectangular map draping
    # ------------------------------------------------------------------
    def load_map(self):
        if self.mesh is None:
            QMessageBox.information(
                self, "No model", "Load a shape model first.")
            return
        here = os.path.dirname(os.path.abspath(__file__))
        maps_dir = os.path.join(here, "maps")
        path, _ = QFileDialog.getOpenFileName(
            self, "Load equirectangular map",
            maps_dir if os.path.isdir(maps_dir) else "",
            "Images (*.png *.jpg *.jpeg *.tif *.tiff *.bmp);;All files (*)")
        if not path:
            return
        try:
            texture = pv.read_texture(path)
        except Exception as e:
            QMessageBox.warning(self, "Could not read map",
                                f"{path}\n\n{e}")
            return
        self._map_texture = texture
        self._map_path = path
        self._apply_map()

    def _apply_map(self):
        """(Re)drape the loaded map over the model at the current
        longitude convention."""
        if self._map_texture is None or self.mesh is None:
            return
        # lat/lon come from the AS-LOADED points (like the lat/lon box), so
        # the drape is unaffected by any principal-axis alignment; the
        # geometry drawn is the displayed mesh. The two share indexing.
        native = (self._mesh_original_points
                  if self._mesh_original_points is not None else self.mesh.points)
        u, v = equirect_texcoords(native, self.sp_map_lon.value(),
                                  east_positive=not self.chk_map_west.isChecked())
        faces = self.mesh.faces.reshape(-1, 4)[:, 1:4]
        pts, faces, u, v = split_texture_seam(self.mesh.points, faces, u, v)
        draped = pv.PolyData(
            pts, np.hstack([np.full((len(faces), 1), 3, np.int64),
                            faces]).ravel())
        draped.active_texture_coordinates = np.column_stack([u, v])

        self.plotter.remove_actor("map_projection", render=False)
        # Unlit: these maps already carry their own shading, and an unlit
        # drape stays readable on terrain angled away from the light.
        self._map_actor = self.plotter.add_mesh(
            draped, texture=self._map_texture, pickable=False,
            lighting=False, name="map_projection")
        # shares geometry with the model actor, so nudge it forward out of
        # z-fighting -- but by less than the paint overlay's -2, so brush
        # strokes still draw on top of the map.
        self._map_actor.GetMapper() \
            .SetRelativeCoincidentTopologyPolygonOffsetParameters(-1, -1)
        self.chk_show_map.setEnabled(True)
        self.chk_show_map.blockSignals(True)
        self.chk_show_map.setChecked(True)
        self.chk_show_map.blockSignals(False)
        self.lbl_map.setText(f"Map: {os.path.basename(self._map_path)}")
        self.plotter.render()

    def _reproject_map(self, *_):
        if self._map_texture is not None:
            self._apply_map()

    def _toggle_map(self, *_):
        if self._map_actor is not None:
            self._map_actor.SetVisibility(self.chk_show_map.isChecked())
            self.plotter.render()

    def remove_map(self):
        self.plotter.remove_actor("map_projection", render=False)
        self._map_actor = None
        self._map_texture = None
        self._map_path = None
        self.chk_show_map.blockSignals(True)
        self.chk_show_map.setChecked(False)
        self.chk_show_map.blockSignals(False)
        self.chk_show_map.setEnabled(False)
        self.lbl_map.setText("No map loaded")
        self.plotter.render()

    # ------------------------------------------------------------------
    # Paint tool (freeform cell painting)
    # ------------------------------------------------------------------
    def _rebuild_face_tree(self):
        self._face_centers = self.mesh.cell_centers().points
        self._face_tree = cKDTree(self._face_centers)

    def _ensure_paint_overlay(self):
        """Translucent overlay sharing the model's geometry, driven by a
        per-cell RGBA array so brush strokes update live with a single
        Modified() call instead of rebuilding geometry every stroke."""
        n = self.mesh.n_cells
        self._paint_rgba = np.zeros((n, 4), np.uint8)
        varr = ns.numpy_to_vtk(self._paint_rgba, deep=0,
                               array_type=VTK_UNSIGNED_CHAR)
        varr.SetNumberOfComponents(4)
        varr.SetName("paint_overlay")
        # must be the ACTIVE cell scalars (not just an extra named array)
        # for the mapper below to pick it up as direct per-cell colour.
        self.mesh.GetCellData().SetScalars(varr)

        mapper = vtkPolyDataMapper()
        mapper.SetInputData(self.mesh)
        mapper.SetScalarModeToUseCellData()
        mapper.ScalarVisibilityOn()
        # the overlay shares exact geometry with the base "model" actor;
        # nudge it forward in depth so it isn't lost to z-fighting.
        mapper.SetRelativeCoincidentTopologyPolygonOffsetParameters(-2, -2)
        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.PickableOff()
        actor.GetProperty().SetInterpolationToFlat()
        self.plotter.renderer.AddActor(actor)
        self._paint_overlay_actor = actor
        self._paint_overlay_array = varr

    def _pick_facet(self):
        if self.mesh is None or self._model_actor is None:
            return -1
        x, y = self.plotter.iren.interactor.GetEventPosition()
        if self._picker.Pick(x, y, 0, self.plotter.renderer):
            return self._picker.GetCellId()
        return -1

    def _brush_facets(self, fid):
        r = self.sp_brush.value()
        hits = self._face_tree.query_ball_point(self._face_centers[fid], r)
        return np.asarray(hits, dtype=np.int64)

    def _between_facets(self, a, b):
        """Facets swept between two consecutive picks, so a fast drag
        doesn't leave gaps between brush stamps."""
        ca, cb = self._face_centers[a], self._face_centers[b]
        r = self.sp_brush.value()
        if np.linalg.norm(cb - ca) > 4.0 * r:
            # picks too far apart to be one continuous stroke (e.g. a fast
            # drag that jumped across the neck) — stamp only the new point
            return self._brush_facets(b)
        pts = ca + np.linspace(0, 1, 8)[:, None] * (cb - ca)
        hits = self._face_tree.query_ball_point(pts, r)
        flat = [np.asarray(hit, np.int64) for hit in hits if len(hit)]
        return np.unique(np.concatenate(flat)) if flat else np.empty(0, np.int64)

    def _stamp(self, facets):
        if len(facets) == 0 or self._paint_rgba is None:
            return
        if self.chk_erase.isChecked():
            self._paint_cells.difference_update(facets.tolist())
            self._paint_rgba[facets] = 0
        else:
            rgb = hex_to_rgb(self._color)
            self._paint_cells.update(facets.tolist())
            self._paint_rgba[facets, 0] = rgb[0]
            self._paint_rgba[facets, 1] = rgb[1]
            self._paint_rgba[facets, 2] = rgb[2]
            self._paint_rgba[facets, 3] = PAINT_ALPHA
        self._paint_overlay_array.Modified()
        self.plotter.render()

    def _toggle_paint_mode(self, checked):
        self.paint_mode = checked
        self.chk_erase.setEnabled(checked)
        self.btn_paint_clear.setEnabled(checked)
        if checked and self.btn_probe.isChecked():
            self.btn_probe.setChecked(False)  # left-click can only mean one thing
        if checked:
            # When editing, edit_selected() has already preloaded the ROI's
            # cells into the overlay -- only a from-scratch ROI starts blank.
            if self._editing_roi_id is None:
                self._paint_cells = set()
                self._reset_overlay()
        else:
            self._finish_paint_session()

    def _clear_paint(self):
        self._paint_cells = set()
        self._reset_overlay()

    def _reset_overlay(self):
        if self._paint_rgba is not None:
            self._paint_rgba[:] = 0
            self._paint_overlay_array.Modified()
            self.plotter.render()

    def _repaint_overlay(self):
        """Redraw the whole overlay from self._paint_cells, in the current
        colour -- used when an edit session opens with an existing ROI's
        cells already painted."""
        if self._paint_rgba is None:
            return
        self._paint_rgba[:] = 0
        if self._paint_cells:
            idx = np.fromiter(self._paint_cells, dtype=np.int64,
                              count=len(self._paint_cells))
            self._paint_rgba[idx, :3] = hex_to_rgb(self._color)
            self._paint_rgba[idx, 3] = PAINT_ALPHA
        self._paint_overlay_array.Modified()
        self.plotter.render()

    def _discard_paint_session(self):
        """Drop any in-progress (not yet saved) paint stroke, e.g. before
        loading a different model out from under it."""
        self._paint_cells = set()
        if self._editing_roi_id is not None:
            self._end_edit_state()
        if self.btn_paint.isChecked():
            self.btn_paint.blockSignals(True)
            self.btn_paint.setChecked(False)
            self.btn_paint.blockSignals(False)
        self.paint_mode = False
        self.chk_erase.setEnabled(False)
        self.btn_paint_clear.setEnabled(False)

    # ---- editing an existing ROI in place ----------------------------
    def edit_selected(self):
        """Reopen the selected ROI in the paint tool so it can be extended
        or trimmed, instead of having to delete it and paint it again."""
        if self.mesh is None:
            return
        if self.paint_mode:
            QMessageBox.information(
                self, "Paint mode is on",
                "Toggle Paint mode off first (that saves whatever is "
                "currently painted), then select an ROI and press "
                "\"Edit selected ROI\".")
            return
        roi = self._selected_roi()
        if roi is None:
            QMessageBox.information(
                self, "No ROI selected",
                "Select an ROI in the list first, then press "
                "\"Edit selected ROI\".")
            return
        cells = np.asarray(roi.get("cells", []), dtype=np.int64)
        if cells.size and (cells.min() < 0 or cells.max() >= self.mesh.n_cells):
            QMessageBox.warning(
                self, "ROI doesn't fit this model",
                f"'{roi['name']}' refers to cells outside the loaded model "
                f"({self.mesh.n_cells} cells) — it was painted on a "
                "different shape model, so it can't be edited here.")
            return

        self._editing_roi_id = roi["id"]
        self._paint_cells = set(cells.tolist())
        # Paint in the ROI's own colour, and load its name/width into the
        # panel fields, since finishing the edit writes those back to it.
        self._color = roi["color"]
        self._style_color_button()
        self.ed_name.setText(roi["name"])
        self.sp_lw.setValue(roi.get("line_width", LINE_WIDTH))
        self._repaint_overlay()
        # its saved outline/fill would sit on top of the live overlay;
        # refresh_visibility hides it for as long as the edit is open.
        self.refresh_visibility()

        self.listw.setEnabled(False)  # the selection must not move mid-edit
        self.btn_edit.setEnabled(False)
        self.btn_cancel_edit.setEnabled(True)
        self._update_paint_button_label()
        self.btn_paint.setChecked(True)  # -> _toggle_paint_mode(True)

    def cancel_edit(self):
        """Leave the edited ROI exactly as it was on disk."""
        if self._editing_roi_id is None:
            return
        self._end_edit_state()
        self._paint_cells = set()
        self._reset_overlay()
        # _editing_roi_id is already cleared, so the resulting
        # _finish_paint_session() call finds nothing painted and does nothing.
        self.btn_paint.setChecked(False)
        self.refresh_visibility()

    def _end_edit_state(self):
        self._editing_roi_id = None
        self.listw.setEnabled(True)
        self.btn_edit.setEnabled(True)
        self.btn_cancel_edit.setEnabled(False)
        self._update_paint_button_label()

    def _update_paint_button_label(self):
        roi = next((r for r in self.rois if r["id"] == self._editing_roi_id),
                   None)
        self.btn_paint.setText("Paint mode" if roi is None
                               else f"Paint mode — editing “{roi['name']}”")

    def _finish_edit(self):
        """Paint mode toggled off during an edit: write the repainted cells
        (plus any name/colour/width tweaks) back onto the same ROI."""
        roi = next((r for r in self.rois if r["id"] == self._editing_roi_id),
                   None)
        cells = sorted(self._paint_cells)
        self._end_edit_state()
        self._paint_cells = set()
        self._reset_overlay()
        if roi is None:
            return
        if not cells:
            QMessageBox.information(
                self, "Nothing left painted",
                f"Every cell of '{roi['name']}' was erased, so it was left "
                "unchanged. Use \"Delete selected ROI\" if you meant to "
                "remove it.")
            self.refresh_visibility()  # un-hide it again
            return
        self._push_undo()
        roi["cells"] = cells
        roi["name"] = self.ed_name.text().strip() or roi["name"]
        roi["color"] = self._color
        roi["line_width"] = self.sp_lw.value()
        item = self._list_item(roi["id"])
        if item is not None:
            item.setText(roi["name"])
        self.rebuild_roi(roi)
        self.plotter.render()
        self._autosave_session()

    def _finish_paint_session(self):
        if self._editing_roi_id is not None:
            self._finish_edit()
            return
        if not self._paint_cells:
            return
        idx = sorted(self._paint_cells)
        self._paint_cells = set()
        self._reset_overlay()

        self._push_undo()
        self._counter += 1
        name = self.ed_name.text().strip() or f"ROI_{self._counter}"
        roi = {
            "id": self._counter,
            "name": name,
            "shape": "freeform",
            "cells": idx,
            "color": self._color,
            "line_width": self.sp_lw.value(),
            "visible": True,
        }
        self.rois.append(roi)
        self._add_list_item(roi)
        self.rebuild_roi(roi)
        self.plotter.remove_actor("latlon_marker", render=False)
        self._latlon_marker_actor = None
        self.ed_name.clear()
        self.plotter.render()
        self._autosave_session()

    # ------------------------------------------------------------------
    # ROI handling
    # ------------------------------------------------------------------
    def _add_list_item(self, roi):
        item = QListWidgetItem(roi["name"])
        item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
        item.setCheckState(QtCore.Qt.Checked if roi["visible"]
                           else QtCore.Qt.Unchecked)
        item.setData(QtCore.Qt.UserRole, roi["id"])
        self.listw.addItem(item)

    def _list_item(self, rid):
        for i in range(self.listw.count()):
            item = self.listw.item(i)
            if item.data(QtCore.Qt.UserRole) == rid:
                return item
        return None

    def _selected_roi(self):
        item = self.listw.currentItem()
        if item is None:
            return None
        rid = item.data(QtCore.Qt.UserRole)
        return next((r for r in self.rois if r["id"] == rid), None)

    def roi_selected(self, item, _prev=None):
        roi = self._selected_roi()
        if roi is None:
            return
        self._updating_ui = True
        self.ed_name.setText(roi["name"])
        self.sp_lw.setValue(roi.get("line_width", LINE_WIDTH))
        self._color = roi["color"]
        self._style_color_button()
        self._updating_ui = False

    def apply_to_selected(self):
        roi = self._selected_roi()
        if roi is None:
            return
        self._push_undo()
        roi["name"] = self.ed_name.text().strip() or roi["name"]
        roi["color"] = self._color
        roi["line_width"] = self.sp_lw.value()
        self.listw.currentItem().setText(roi["name"])
        self.rebuild_roi(roi)
        self._autosave_session()

    def choose_color(self):
        c = QColorDialog.getColor()
        if c.isValid():
            self._color = c.name()
            self._style_color_button()

    # ---- undo ---------------------------------------------------------
    # Whole-list snapshots rather than per-operation inverses: every saved
    # change (new ROI, finished edit, apply, delete) stashes a deep copy of
    # self.rois first, and Undo swaps the last one back in. ROI cell lists
    # are a few hundred KB at worst, so a short stack of copies is cheap,
    # and one restore path covers every kind of change.
    UNDO_LIMIT = 20

    def _push_undo(self):
        self._undo_stack.append(copy.deepcopy(self.rois))
        del self._undo_stack[:-self.UNDO_LIMIT]
        self.btn_undo.setEnabled(True)

    def undo(self):
        """Restore the ROI list to the last snapshot and auto-save it, so
        the file on disk always matches what's shown."""
        if not self._undo_stack:
            return
        if self.paint_mode or self._editing_roi_id is not None:
            QMessageBox.information(
                self, "Paint mode is on",
                "Finish or cancel the current painting/edit first, then "
                "press Undo.")
            return
        for rid in list(self.actors):
            self._remove_actors(rid)
        self.rois = self._undo_stack.pop()
        self._counter = max((r["id"] for r in self.rois), default=0)
        self.listw.blockSignals(True)
        self.listw.clear()
        for roi in self.rois:
            self._add_list_item(roi)
        self.listw.blockSignals(False)
        self.rebuild_all()
        self.plotter.render()
        self._autosave_session()
        self.btn_undo.setEnabled(bool(self._undo_stack))

    def delete_selected(self):
        if self._editing_roi_id is not None:
            QMessageBox.information(
                self, "An ROI is being edited",
                "Finish the edit (Paint mode off) or press \"Cancel edit\" "
                "before deleting an ROI.")
            return
        roi = self._selected_roi()
        if roi is None:
            return
        self._push_undo()
        self._remove_actors(roi["id"])
        self.rois.remove(roi)
        self.listw.takeItem(self.listw.currentRow())
        self.plotter.render()
        self._autosave_session()

    def _remove_actors(self, rid):
        acts = self.actors.pop(rid, {})
        for a in acts.values():
            self.plotter.remove_actor(a, render=False)

    def _clear_rois(self):
        for rid in list(self.actors):
            self._remove_actors(rid)
        self.rois = []
        self.listw.clear()
        self._counter = 0

    def rebuild_roi(self, roi):
        self._remove_actors(roi["id"])
        try:
            patch, border = build_roi_geometry(self.mesh_normals, roi)
        except Exception as e:
            QMessageBox.warning(self, "ROI error",
                                f"Could not build ROI '{roi['name']}':\n{e}")
            return
        b_act = self.plotter.add_mesh(
            border, color=roi["color"],
            line_width=roi.get("line_width", LINE_WIDTH),
            pickable=False, name=f"border_{roi['id']}")
        f_act = self.plotter.add_mesh(
            patch, color=roi["color"], opacity=FILL_OPACITY,
            pickable=False, name=f"fill_{roi['id']}")
        self.actors[roi["id"]] = {"border": b_act, "fill": f_act}
        self.refresh_visibility()

    def move_selected(self, delta):
        """Nudge the selected ROI one place up (-1) or down (+1)."""
        if self._editing_roi_id is not None:
            QMessageBox.information(
                self, "An ROI is being edited",
                "Finish the edit (Paint mode off) or press \"Cancel edit\" "
                "before reordering the list.")
            return
        row = self.listw.currentRow()
        new_row = row + delta
        if row < 0 or not 0 <= new_row < self.listw.count():
            return
        # takeItem/insertItem move the very same item object, so its check
        # state and its ROI id ride along untouched.
        self.listw.blockSignals(True)
        item = self.listw.takeItem(row)
        self.listw.insertItem(new_row, item)
        self.listw.setCurrentItem(item)
        self.listw.blockSignals(False)
        self._sync_roi_order()

    def _sync_roi_order(self):
        """Reorder self.rois to match the list and save. The list's order is
        the order ROIs are written to (and read back from) the session."""
        by_id = {r["id"]: r for r in self.rois}
        order, seen = [], set()
        for i in range(self.listw.count()):
            rid = self.listw.item(i).data(QtCore.Qt.UserRole)
            roi = by_id.get(rid)
            if roi is not None and rid not in seen:
                order.append(roi)
                seen.add(rid)
        # never drop an ROI that somehow has no row: keep it, at the end
        order.extend(r for r in self.rois if r["id"] not in seen)
        unchanged = [r["id"] for r in order] == [r["id"] for r in self.rois]
        self.rois = order
        if not unchanged:  # a drop that changed nothing needn't rewrite the file
            self._autosave_session()

    def _set_all_visible(self, visible):
        """Tick/untick every ROI at once — mainly for a session loaded back
        in, which comes in with everything that was visible when it was
        saved, and is usually far too much to read at once."""
        state = QtCore.Qt.Checked if visible else QtCore.Qt.Unchecked
        # one refresh at the end, rather than one per item (itemChanged is
        # wired to refresh_visibility, which re-renders the whole scene).
        self.listw.blockSignals(True)
        for i in range(self.listw.count()):
            self.listw.item(i).setCheckState(state)
        self.listw.blockSignals(False)
        self.refresh_visibility()

    def rebuild_all(self):
        for roi in self.rois:
            self.rebuild_roi(roi)

    def refresh_visibility(self, *_):
        show_fill = self.chk_fill.isChecked()
        for i in range(self.listw.count()):
            item = self.listw.item(i)
            rid = item.data(QtCore.Qt.UserRole)
            vis = item.checkState() == QtCore.Qt.Checked
            roi = next((r for r in self.rois if r["id"] == rid), None)
            if roi is not None:
                roi["visible"] = vis
            # the ROI currently open in the paint tool is drawn by the live
            # overlay instead, so its saved outline/fill stays hidden (its
            # stored "visible" flag above is untouched by that).
            shown = vis and rid != self._editing_roi_id
            acts = self.actors.get(rid, {})
            if "border" in acts:
                acts["border"].SetVisibility(shown)
            if "fill" in acts:
                acts["fill"].SetVisibility(shown and show_fill)
        self.plotter.render()

    # ------------------------------------------------------------------
    # Session save / load
    # ------------------------------------------------------------------
    def _session_data(self):
        return {
            "model_path": self.model_path,
            "transform": (self.transform.tolist()
                          if self.transform is not None else None),
            "rois": self.rois,
        }

    def _write_session(self, path):
        with open(path, "w") as f:
            json.dump(self._session_data(), f, indent=2)

    def _default_session_path(self):
        here = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(here, "ROI_data")
        os.makedirs(data_dir, exist_ok=True)
        base = (os.path.splitext(os.path.basename(self.model_path))[0]
                if self.model_path else "session")
        return os.path.join(data_dir, f"{base}_rois.json")

    def _autosave_session(self):
        """Called after every new ROI is added. Writes to whatever file is
        already established (via a manual Save/Load) or, the first time,
        a sensible default under ROI_data/ -- no dialog, so it never
        interrupts painting."""
        if self.mesh is None:
            return
        if self._session_path is None:
            self._session_path = self._default_session_path()
        try:
            self._write_session(self._session_path)
        except Exception as e:
            QMessageBox.warning(self, "Auto-save failed",
                                f"{self._session_path}\n\n{e}")
            return
        self.lbl_autosave.setText(f"Auto-saved to: {self._session_path}")

    def save_session(self):
        if self.mesh is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save session", self._session_path or "rois.json",
            "JSON (*.json)")
        if not path:
            return
        self._session_path = path
        self._write_session(path)
        self.lbl_autosave.setText(f"Auto-saved to: {self._session_path}")

    def load_session(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load session", "", "JSON (*.json)")
        if not path:
            return
        with open(path) as f:
            data = json.load(f)
        model = data["model_path"]
        try:
            self.load_model(model, transform=data.get("transform"))
        except Exception:
            # a session written on another machine stores that machine's
            # path (or a bare filename): try the same filename next to this
            # script before asking
            here = os.path.dirname(os.path.abspath(__file__))
            fallback = os.path.join(here, os.path.basename(str(model)))
            try:
                self.load_model(fallback, transform=data.get("transform"))
            except Exception:
                QMessageBox.information(
                    self, "Model not found",
                    f"Couldn't open:\n{model}\nPlease locate the model file.")
                model, _ = QFileDialog.getOpenFileName(
                    self, "Locate shape model", "", "Meshes (*)")
                if not model:
                    return
                self.load_model(model, transform=data.get("transform"))
        self.rois = data["rois"]
        self._counter = max((r["id"] for r in self.rois), default=0)
        for roi in self.rois:
            self._add_list_item(roi)
        self.rebuild_all()
        # keep updating the file we just loaded from, going forward
        self._session_path = path
        self.lbl_autosave.setText(f"Auto-saved to: {self._session_path}")

    def screenshot(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save screenshot", "view.png", "PNG (*.png)")
        if path:
            self.plotter.screenshot(path)


# ======================================================================

def _enable_hang_diagnostics():
    """Let a wedged window be diagnosed instead of guessed at.

    With this in place, pressing Ctrl-\\ in the terminal that launched the app
    (or `kill -QUIT <pid>`) prints a Python traceback for every thread showing
    exactly where it is stuck, and carries on running. Costs nothing when
    nothing goes wrong."""
    try:
        import faulthandler
        import signal
        faulthandler.enable()
        for sig in ("SIGQUIT", "SIGUSR1"):
            if hasattr(signal, sig):
                faulthandler.register(getattr(signal, sig), chain=True)
    except Exception:
        pass        # diagnostics are a convenience, never a reason to fail


def main():
    _enable_hang_diagnostics()
    app = QtWidgets.QApplication(sys.argv)
    if len(sys.argv) > 1:
        model = sys.argv[1]
    else:
        here = os.path.dirname(os.path.abspath(__file__))
        candidate = os.path.join(here, DEFAULT_MODEL)
        model = candidate if os.path.isfile(candidate) else None
    win = ROIManager(model)
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
