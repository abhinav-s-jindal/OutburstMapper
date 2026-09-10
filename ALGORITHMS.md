# OutburstMapper: algorithm details

This document describes the implementation in
[`outburst_mapper.py`](outburst_mapper.py). Installation, data locations,
and controls are described in the [README](README.md). The source reviewed
for this description is [commit fb16052](https://github.com/abhinav-s-jindal/OutburstMapper/blob/fb16052e1332b3996e99024f03e0619a1a45fd7d/outburst_mapper.py).
This documentation adds no changes to the numerical implementation.

## Purpose, inputs, and outputs

OutburstMapper displays mission images alongside a triangulated shape
model, reconstructs the spacecraft camera geometry with SPICE, projects
images onto visible model facets, and stores regions of interest (ROIs)
painted by the user. Projection supplies a geometric context for locating
an outburst; identifying its source and drawing the boundary require
scientific interpretation. The tool does not automatically detect plumes,
triangulate their three-dimensional structure, or determine a causal link
between an outburst and surface change.

Inputs are a shape model; optionally a Rosetta SPICE metakernel and its
referenced kernels; an ISIS cube or PDS3 mission image and label; an
existing ROI session; and/or an equirectangular map. Outputs are a rendered
view, optional PNG screenshots, and JSON sessions containing model and ROI
information. Projected images are display textures, not newly calibrated
science rasters.

## 1. Shape model and coordinate system

`ROIManager.load_model` reads the model with PyVista, extracts its surface,
and triangulates it. It retains a copy of the native vertex coordinates
and applies a saved model transform, if present. The resulting triangle
indices identify the cells stored in each ROI. A `cKDTree` indexes triangle
centres for painting; a separate VTK `vtkOBBTree` is constructed for
line-of-sight tests during image projection.

SPICE geometry uses the body-fixed frame `67P/C-G_CK`, spacecraft
`ROSETTA`, and target `CHURYUMOV-GERASIMENKO`. `spkpos` returns positions in
kilometres. For mission projection, the loaded mesh must therefore have
compatible coordinates, units, and orientation. Loading an arbitrary mesh
is supported for visualization, but does not automatically register it to
the Rosetta frame or convert its units.

## 2. Image reading and display values

`load_image_any` returns the image array, validity mask, instrument NAIF
identifier, and acquisition start time.

- **ISIS `.cub`:** `read_isis_cube` parses the PVL label, reads the pixel
  type and byte order, reconstructs tiled or band-sequential storage, and
  uses the first band as `float32`. Instrument and time come from
  `NaifFrameCode` and `StartTime`. For ISIS `Real` pixels,
  `isis_valid_mask` rejects special-pixel bit patterns at or above
  `0xFF7FFFFA`; other pixel types do not receive this special-pixel test.
  This is a minimal reader and does not apply the label's Base/Multiplier
  calibration to the stored values.
- **PDS3 `.IMG` / `.LBL`:** `read_pds3_image` follows the primary `^IMAGE`
  pointer in an attached or detached label and reads the declared sample
  type, byte order, lines, and samples. Auxiliary image objects are
  ignored. Rows are reversed when `LINE_DISPLAY_DIRECTION` is `DOWN`
  (also the default when absent). `load_image_any` additionally reverses
  columns for WAC and both NAVCAM heads. The PDS3 validity mask accepts
  finite values; it does not interpret every possible product-specific
  special-pixel flag.

The image panel linearly stretches selected valid-pixel percentiles to
8-bit grayscale, clips to 0-255, and sets invalid pixels to black. A
user-drawn stretch box instead supplies the minimum and maximum valid
values within that box. `project_cube` uses this displayed 8-bit image as
its RGB texture. Contrast stretching is for visual interpretation, not
radiometric measurement.

## 3. SPICE camera geometry

`load_spice_kernels` rewrites the metakernel's single `PATH_VALUES` entry
to the local kernels directory in a temporary file, clears the SPICE
kernel pool, and loads that file. Missing coverage or unresolvable
instrument geometry produces an error rather than an estimated pose.

`spice_camera_basis` constructs a pinhole camera as follows:

1. Convert the image start time to SPICE ephemeris time with `str2et`.
2. Obtain the spacecraft position **C** in the comet frame using
   `spkpos(..., "NONE", ...)`. This uses geometric positions without
   light-time or stellar-aberration corrections.
3. Read the instrument boresight and FOV reference vector, and the FOV
   reference and cross half-angles, from the instrument kernels.
4. Rotate the normalized boresight and reference vector into the comet
   frame using `pxform` at the image time.
5. Set **f** to the rotated boresight. Remove the component of the rotated
   reference vector parallel to **f**, then normalize it to obtain **u**
   (camera up). Set camera right **r** = **f** cross **u**.
6. Use the cross half-angle as horizontal half-angle alpha_h and the
   reference half-angle as vertical half-angle alpha_v.

The mapping assumes this rectangular pinhole geometry. It does not fit a
pose to image landmarks, explicitly model lens distortion, or integrate
camera motion over an exposure.

## 4. Projecting an image onto the mesh

For each mesh vertex **P**, `project_cube_onto_mesh` computes:

```text
q = P - C
z = dot(q, f)
x = dot(q, r)
y = dot(q, u)

image_u = (x / (z * tan(alpha_h)) + 1) / 2
image_v = 1 - (y / (z * tan(alpha_v)) + 1) / 2
```

A vertex is in frame only if `z > 0` and both normalized image coordinates
are in `[0, 1]`. Here `image_v = 0` addresses the top image row.

For every in-frame vertex, the algorithm tests the line segment from the
camera to `C + 0.999 * (P - C)` against the mesh's `vtkOBBTree`. A vertex is
accepted if that shortened segment has no mesh intersection. Shortening
avoids counting the destination surface itself as an obstruction.
`ROIManager.project_cube` retains a triangle only when **all three** of
its vertices are accepted. It assigns VTK texture coordinates
`(image_u, 1 - image_v)` to account for VTK's texture orientation.

This is a camera-visibility test, not a calculation of solar illumination
or physical shadowing. Vertex-based acceptance can omit partly visible
triangles near the image edge and can miss obstructions affecting only a
triangle's interior. The 0.999 segment factor is a numerical tolerance,
not a physical clearance. Spatial accuracy remains limited by mesh
resolution, kernel geometry, image orientation, and the user's mapping.

## 5. Linked spacecraft and image views

`view_from_image` stores the camera basis, FOV tangents, and image width
and height. `_pixels_to_tangent` and `_tangent_to_pixels` convert between
image rectangles and camera angular coordinates. The two panels share a
view centre and vertical extent, with horizontal extents adjusted to each
panel's aspect ratio. Linked zoom and pan preserve the spacecraft position
and boresight; ordinary free orbit is restored when linking is disabled.
This changes the view, not the stored image pixels or painted facet IDs.

## 6. Painting and storing ROIs

A VTK cell picker identifies the triangle under the cursor. For brush
radius `R`, `_brush_facets` selects triangle centres within Euclidean
three-dimensional distance `R` of the picked triangle centre. The radius
uses the loaded model's coordinate units. This is a spatial neighbourhood,
not a geodesic or mesh-connectivity flood fill; nearby disconnected
surfaces can fall inside the same brush volume.

To fill gaps during a drag, `_between_facets` samples eight positions
between consecutive picked centres and unions their radius searches.
When the picked centres are more than `4 * R` apart, it stamps only at the
new centre. Painting adds cell IDs to a set; erasing removes them.

When painting ends, the sorted IDs become a new ROI or replace the cells
of the ROI being edited. `build_roi_geometry` extracts those cells and
then their boundary edges, which supply the displayed outline. The ROI
record stores its ID, name, shape (`freeform`), cell list, colour, line
width, and visibility. Display names and colour choices do not change the
underlying geometry.

`_session_data` stores `model_path`, an optional transform, and the ordered
ROI list. The session does **not** store the loaded image, image time,
SPICE kernels, projection texture, camera pose, or a complete action log.
Preserve the matching model and the image/kernel provenance separately.
Cell IDs require the same triangulated mesh and cell ordering; a different
resolution or reindexed model is not interchangeable. Loading a session
sets that file as the autosave destination, so work from a copy when
preserving an archived session. Undo retains up to 20 ROI-state snapshots.

## 7. Equirectangular maps

`equirect_texcoords` uses native model coordinates to calculate
planetocentric latitude `asin(z / norm(P))` and longitude `atan2(y, x)`.
Longitude is optionally negated for a west-positive map. In degrees:

```text
map_u = ((longitude - (centre_longitude + 180)) mod 360) / 360
map_v = (latitude + 90) / 180
```

`split_texture_seam` duplicates a triangle corner when its `u` coordinate
is more than 0.5 below that triangle's maximum `u`, assigning the duplicate
`u + 1`. This allows the repeating texture to cross the longitude seam
without interpolating backwards across the full map. This mapping is
separate from spacecraft-image projection.

## End-to-end pseudocode

```text
mesh = read model, extract surface, triangulate, apply optional transform
build triangle-centre index
if restoring a session:
    restore ROI facet IDs and presentation properties on the matching mesh

load the Rosetta kernels when mission camera features are required
read image band, validity mask, instrument, and start time
stretch valid values to an 8-bit display image
construct the SPICE camera position, basis, and FOV half-angles

for each mesh vertex:
    calculate normalized image coordinates with the pinhole equations
    reject vertices behind the camera or outside the image
    test remaining vertices for mesh occlusion along the viewing segment
retain triangles with three accepted vertices
texture the retained patch with the displayed image

while the user paints or erases:
    pick a triangle and find neighbouring triangle centres
    update the set of selected cell IDs
on completing an ROI:
    extract its patch and boundary, and update the session JSON
```

## Reproducibility record

Keep the software commit, dependency versions, original model file,
kernel set, exact image products and labels, and ROI JSON together. Record
which image supports each mapped boundary. A screenshot helps communicate
the interpretation, but is not a substitute for those inputs. Reopening a
session reproduces its facet selections; evaluating the geological
interpretation also requires the source imagery and viewing geometry.
