# SWOT Landfast Ice Morphology Index Probe

## Objective

Prototype a candidate stability-relevant morphology index from SWOT KaRIn 2D height/freeboard and backscatter fields for Beaufort Sea landfast-ice analysis.

Core question:

Can local roughness, height-gradient, and backscatter-texture metrics from SWOT KaRIn identify coherent landfast-ice morphology features that are plausibly relevant to mechanical stability?

## Novelty Position

This prototype is not claiming that SWOT ridge mapping itself is new. Existing work has already demonstrated SWOT's utility for sea-ice morphology, freeboard, lead/floe discrimination, velocity, and ridging in Beaufort Sea landfast ice.

The seed idea is that SWOT-derived 2D morphology may become useful as an input to landfast-ice stability, persistence, anchoring, or breakout-susceptibility analysis. The present notebook only develops and stress-tests the morphology index; future work must validate it against SAR persistence, breakout timing, bathymetry, or models.

## Dataset

This prototype uses the published Beaufort Sea SWOT-ICESat-2 dataset from Zenodo:

https://zenodo.org/records/15238132

Primary prototype case:

- SWOT: `UMD_SWOT_024_BeaufortSea_Study_Region_02_20230409.nc`
- SWOT swath: `left`
- ICESat-2: `UMD_IC2_0174_BeaufortSea_Study_Region_02_20230401.nc`
- ICESat-2 beam: `gt2l`

The selected Region 02 SWOT and ICESat-2 files overlap spatially, but the ICESat-2 pass is 8 days earlier than the SWOT scene. The notebook treats that comparison as partial spatial support, not definitive validation.

## Method

The notebook:

1. Downloads and inventories the Zenodo NetCDF files.
2. Opens the selected SWOT swath and ICESat-2 beam without assuming variable names in advance.
3. Builds a provisional ice-like analysis mask from finite height/backscatter and broad height bounds.
4. Infers SWOT grid spacing from neighboring latitude/longitude cells.
5. Computes local height roughness, height gradient magnitude, and local sigma0 texture for 500 m, 1 km, and 2 km windows.
6. Combines standardized metrics into a stability morphology index:

   `SMI = z(local height std) + z(height gradient magnitude) + z(local sigma0 std)`

7. Tests height-only, backscatter-only, and combined SMI variants.
8. Thresholds the SMI into high-index roughness features, treated as candidate stability-relevant morphology only.
9. Collocates the SMI to ICESat-2 points with nearest-neighbor sampling in a local projected coordinate system.
10. Runs algorithmic uncertainty checks for window size, threshold choice, component choice, swath-side differences, and row/column structure.

## Figures

### Main result: SWOT height, sigma0, SMI, and high-index roughness features
![Four-panel SMI map with ICESat-2 overlay](figures/four_panel_smi_region02_left_20230409.png)

### ICESat-2 height cross-section through a high-SMI feature
![ICESat-2 cross-section vs SWOT SMI](figures/icesat2_smi_cross_section_region02_left_gt2l.png)

### SMI distribution: high-index pixels vs background ice-like pixels
![SMI histogram](figures/smi_histogram_region02_left_20230409.png)

### Algorithmic uncertainty — window size and threshold sensitivity
![Window and threshold sensitivity](figures/smi_window_threshold_sensitivity_region02_left_20230409.png)

### Algorithmic uncertainty — component sensitivity (height-only, backscatter-only, combined)
![Component sensitivity](figures/smi_component_sensitivity_region02_left_20230409.png)

### Artifact checks — left/right swath behavior and row/column SMI patterns
![Artifact checks](figures/smi_artifact_checks_region02_20230409.png)

### Multi-date SMI persistence across three Region 02 scenes (April 9, April 10, June 16)
![SMI persistence map](figures/smi_persistence_region02_left.png)

## Outputs

Generated figures are saved in `figures/`:

- `four_panel_smi_region02_left_20230409.png`: SWOT height, sigma0, SMI, and high-index roughness features with ICESat-2 overlay.
- `icesat2_smi_cross_section_region02_left_gt2l.png`: ICESat-2 height cross-section through a high-SMI feature with nearest SWOT SMI.
- `smi_histogram_region02_left_20230409.png`: SMI distributions for high-index pixels versus background ice-like pixels.
- `smi_window_threshold_sensitivity_region02_left_20230409.png`: high-index area and feature count versus window size and threshold.
- `smi_component_sensitivity_region02_left_20230409.png`: high-index area and feature count for height-only, backscatter-only, and combined SMI.
- `smi_artifact_checks_region02_20230409.png`: visual checks for left/right swath behavior and row/column SMI patterns.
- `smi_persistence_region02_left.png`: SMI maps for all three Region 02 dates side-by-side plus a persistence count map (0–3 dates high-index), with pairwise Jaccard similarity table.

Current 1 km-window prototype result:

- Valid SWOT analysis pixels: 159,046.
- Top 10 percent high-index threshold: SMI >= 2.53.
- High-index pixels: 15,886, about 951 km2.
- Connected high-index features after a 4-pixel minimum: 448.
- ICESat-2 matched points: 106,422.
- Top 5 percent sampled SMI points have higher mean ICESat-2 20 m versus 500 m height difference than background sampled points, consistent with enhanced along-track morphology but not proof of grounding, anchoring, thickness, or stability.

## Interpretation

High-index areas are candidate morphology features only. They are not confirmed grounded ridges, anchor points, thickness anomalies, or stable landfast-ice zones.

The method is useful if it produces coherent, robust, interpretable morphology features that can be carried forward into future validation against breakout or persistence data.

## Limitations

- No snow depth, density, radar penetration, local sea-surface reconstruction, or uncertainty propagation is included.
- No thickness product is generated.
- No external landfast-ice mask is used.
- The primary ICESat-2 comparison is spatially overlapping but temporally offset by 8 days.
- SWOT artifacts remain plausible, including swath-side differences, along-track/across-track structure, temporal height errors, calibration issues, and outlier height cells.
- The provisional mask is not a sea-ice or landfast-ice classifier.

## Future Validation Layers

1. Breakout prediction layer: compare SWOT roughness/ridge/anchor metrics before breakup against observed breakup locations or dates.
2. SAR persistence layer: use Sentinel-1 to classify landfast versus mobile ice before and after the SWOT scene, then test whether SWOT morphology predicts which zones remain attached.
3. Bathymetric grounding layer: combine SWOT roughness features with bathymetry to identify candidate grounded ridges or shallow-water anchoring zones, then test whether high-roughness features occur at plausible grounding depths.
4. Model-parameterization layer: convert SWOT-derived roughness, ridge density, or morphology-index statistics into parameters relevant to coastal drag, internal ice strength, or landfast-ice anchoring in models.
5. Algorithmic uncertainty layer: quantify how stable the morphology index is under swath-side differences, window size, threshold choice, score-component choice, and known SWOT sea-ice artifacts.

The current notebook prioritizes the algorithmic uncertainty layer because it is doable immediately with the existing dataset.

## Reproducibility

Run the notebook from this folder:

```bash
jupyter nbconvert --to notebook --execute --inplace swot_landfast_morphology_probe.ipynb --ExecutePreprocessor.timeout=900
```

The notebook caches NetCDF files in `data/raw/` and regenerates PNGs in `figures/`.

Main Python dependencies: `xarray`, `numpy`, `pandas`, `scipy`, `matplotlib`, `pyproj`, `netCDF4`, `scikit-image`, and `requests`.
