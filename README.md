# OutburstMapper

An interactive 3D tool for mapping cometary outburst source locations and
surface-change regions of interest (ROIs) onto a comet shape model.

OutburstMapper was built to survey the outbursts of comet
67P/Churyumov–Gerasimenko observed by the Rosetta mission and the surface
changes associated with them: it renders the 67P shape model, places the
virtual camera exactly where Rosetta's OSIRIS or NAVCAM cameras were at any
moment of the mission (via SPICE), displays and projects the actual mission
imagery onto the terrain, and lets you paint, edit, and organize named
regions directly on the surface. The tool itself is generic — any
triangulated shape model (`.obj`/`.ply`/`.stl`/`.vtk`) can be loaded — but
the SPICE camera features are wired to the Rosetta mission kernels.

This tool and its data archive are companions to the paper *"How the Comet
Crumbles: Mass Wasting Drives Outbursts on Comet
67P/Churyumov-Gerasimenko"* (DOI to be added upon publication): the
archived ROI session contains every outburst footprint and surface-change
region mapped for that work, and the imagery on which the boundaries were
drawn.

Built with PyVista/VTK and Qt (PyQt5).

## Installation

Requires Python 3.10+.

```bash
git clone https://github.com/abhinav-s-jindal/OutburstMapper.git
cd OutburstMapper
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

If a pinned version fails to install on your platform, drop the `==version`
pins in `requirements.txt` and re-run — the tool does not depend on exact
patch versions.

**macOS note:** if you edit the code, always import VTK classes through
narrow `vtkmodules.*` submodules (as the existing code does), never via a
bare `import vtk`. Environments that carry both PyQt5 (Qt5) and PySide6
(Qt6) will otherwise pull two Qt runtimes into one process and crash.

## Data

The shape model, SPICE kernels, mission images, and our ROI session are too
large for GitHub and are archived on Zenodo:

> **DOI: [10.5281/zenodo.22133944](https://doi.org/10.5281/zenodo.22133944)**

Download the archive and unpack it into the repository root as a `data/`
folder, so the layout becomes:

```
OutburstMapper/
├── outburst_mapper.py
└── data/                        # ← the unpacked Zenodo archive
    ├── sessions/
    │   └── 67P_outburst_session.json    # our mapped ROIs & outburst footprints
    ├── shape_model/
    │   └── cg-dlr_spg-shap7-v1.0_125Kfacets.obj   # SHAP7 67P shape model
    ├── spice_kernels/
    │   └── kernels/mk/ROS_OPS.TM        # Rosetta SPICE kernels (metakernel)
    ├── maps/                            # equirectangular maps
    └── ROI_data/                        # per-ROI mission imagery
        ├── roi1/
        │   ├── outburst_img/    # an image showing this ROI's outburst
        │   └── surface_img/     # an image showing the surface change
        ├── roi2/
        │   └── ...
        └── ...
```

- **`data/sessions/`** holds the ROI session: load it in the app (see
  [Loading the ROI session](#loading-the-roi-session)) to restore every
  outburst footprint and surface-change boundary from the paper.
- **`data/ROI_data/roi<N>/`** holds, for each ROI, one representative
  OSIRIS/NAVCAM image of the outburst (`outburst_img/`) and one of the
  resulting surface change (`surface_img/`) — the frames on which that
  ROI's boundary was drawn, so the basis of each mapping can be inspected
  directly (`.cub`, `.IMG`/`.LBL`; open them via **Load image…**).
- **`data/maps/`** holds two equirectangular maps ready to drape (see
  [Draping an equirectangular map](#draping-an-equirectangular-map)):
  - `67P_regions_equirectangular.jpg` — the geomorphological region map of
    67P, from El-Maarry et al. 2015
    ([doi:10.1051/0004-6361/201525723](https://doi.org/10.1051/0004-6361/201525723))
    and El-Maarry et al. 2016
    ([doi:10.1051/0004-6361/201628634](https://doi.org/10.1051/0004-6361/201628634)).
  - `JB_OB_Locs.png` — the outburst source-location map of Vincent et al.
     2016 ([doi:10.1093/mnras/stw2409](https://doi.org/10.1093/mnras/stw2409));
     its numbered markers correspond to the `roi<N>_ob_loc_jb<M>`
     outburst locations in the shipped session.
- The SPICE kernels can alternatively be fetched from the NAIF archive
  (<https://naif.jpl.nasa.gov/pub/naif/pds/data/ro_rl-e_m_a_c-spice-6-v1.0/rossp_1000/>); the app
  expects the metakernel at `data/spice_kernels/kernels/mk/ROS_OPS.TM`
  next to the script, and rewrites the metakernel's `PATH_VALUES` to an
  absolute path internally, so no manual editing of the kernel files is
  needed.

## Running

```bash
python outburst_mapper.py                # auto-loads the default 67P model
python outburst_mapper.py my_model.obj   # or any other shape model
```

On startup the app auto-loads the shape model and the SPICE metakernel
from the unpacked `data/` folder. Without it you can still load any model
via **Load model…** and any metakernel via **Load SPICE kernels…**.

## Loading the ROI session

Our full set of mapped outburst locations and change regions comes with
the data archive: press **Load session…** and pick
`data/sessions/67P_outburst_session.json`. The session restores every ROI
(name, painted facets, colour, border width, list order, visibility) onto
the model it was painted on — the stored model path resolves automatically
to the model file in `data/shape_model/`.

Things to know once a session is loaded:

- A loaded session usually shows *all* ROIs at once. **Deselect all**, then
  tick just the ones you want to inspect.
- Every change (new ROI, edit, rename, delete, reorder) **auto-saves back
  to the loaded file**. Work on a copy if you want to keep the shipped
  session pristine. **Undo** steps back through the last 20 changes.
- Naming convention in our session: `roi<N>` rows are surface-change
  regions; `roi<N>_ob_loc_*` rows are the outburst source locations
  associated with region `<N>`.

## Feature guide

### Navigating

Left-drag orbits, right-drag/wheel zooms, middle-drag pans. The six axis
buttons (**+X … −Z**) snap the camera along the model's body-fixed axes.
**Show axes** draws the X/Y/Z axis lines (red/green/blue).

### Painting ROIs

Toggle **Paint mode** and drag across the surface to paint facets; set the
brush radius with **Brush r**, tick **Erase** to remove misplaced strokes,
and paint in as many separate strokes (or disconnected patches) as needed.
Toggling Paint mode **off** saves the painted facets as a new ROI with the
current name/colour/width — and auto-saves the session. To modify an
existing ROI, select it and press **Edit selected ROI**: it reopens in the
paint tool with its facets pre-painted, and saving writes back to the same
ROI in place.

### The ROI list

Checkboxes control visibility; **Select all**/**Deselect all** flip
everything at once. Drag rows (or **Move up**/**Move down**) to reorder —
the order is what's saved. **Apply to selected** pushes the current
name/colour/border width onto the selected ROI. **Show translucent fills**
adds a see-through fill inside every visible ROI's outline.

### SPICE camera views

Enter any UTC time of the mission, pick **NAC**, **WAC**, or **NAVCAM**,
and the 3D camera moves to that instrument's exact position, boresight,
and field of view at that moment. Times outside SPICE coverage (e.g. the
2011–2014 hibernation cruise) produce a clear error rather than a wrong
view.

### Loading mission images

**Load image (.cub/.IMG/.LBL)…** opens an ISIS cube or a raw PDS3
OSIRIS/NAVCAM frame (open the `.LBL` for detached-label NAVCAM products)
in a panel beside the 3D view — e.g. any frame from
`data/ROI_data/roi<N>/outburst_img/`. The panel has percentile contrast
stretching, a drag-a-box stretch tool, and wheel zoom. **Close image**
returns to the single-panel view.

### View from spacecraft, and linked zoom

**View from spacecraft** reads the loaded image's own instrument and
acquisition time from its label and places the 3D camera exactly where the
image was taken — then links the two panels: zooming or panning either one
(any mouse button) shows the same ground in the other. While linked,
navigation that would move the camera off the spacecraft's true position
(orbiting, dollying) is remapped to pan/zoom of the shared window; untick
**Link zoom** to navigate freely again.

### Projecting an image onto the model

**Project onto model** drapes the loaded image onto the terrain using the
exact SPICE camera geometry, with per-vertex occlusion testing so the
texture doesn't wrap onto the far side or across shadowing terrain.
**Show projected image** toggles it; **Remove projection** discards it.

### Draping an equirectangular map

**Load map…** wraps a whole-body lat/lon map onto the model by
planetocentric latitude/longitude, in the same frame as the lat/lon lookup
box. The two maps shipped in the archive's `data/maps/` folder both use the
app's default convention — 0° at the central column, east-positive, +90°
at the top row — so they drape correctly with no adjustment:

- **`67P_regions_equirectangular.jpg`** — the geomorphological region map
  (El-Maarry et al. 2015, 2016): drape it to read off which named region
  (Imhotep, Ma'at, Anhur, …) an ROI or outburst location falls in.
- **`JB_OB_Locs.png`** — the Vincent et al. (2016) outburst
  source-location map: its numbered markers line up with the
  `roi<N>_ob_loc_jb<M>` ROIs in the shipped session, useful for
  cross-checking a mapped location against the original survey.

For maps drawn to another convention, **Lon at centre** shifts the central
meridian (use 180 for a map whose left edge is 0°) and **West-positive**
flips the longitude direction; both re-drape immediately. The drape
renders *under* ROI outlines and paint strokes, so you can paint straight
on top of it, and it is dropped when a different shape model loads.
**Show map** toggles it; **Remove map** discards it.

### Finding and probing locations

The **lat/lon box + Go** centres the camera over any planetocentric
lat/lon (east-positive, model's native frame) and drops a marker. **Probe
mode** turns clicks into readouts: the clicked facet's lat/lon/radius and
outward normal, drawn as a line, with the lat/lon boxes updated so **Go**
returns there.

### Screenshots

**Save screenshot…** writes the current 3D view to a PNG.

## License

MIT — see [LICENSE](LICENSE).
