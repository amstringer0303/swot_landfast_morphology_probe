"""Utilities for a SWOT landfast-ice morphology-index prototype.

The functions in this module deliberately treat high-index areas as candidate
morphology features, not confirmed grounded ridges, anchor points, thickness
products, or stable landfast-ice zones.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests
import xarray as xr
from netCDF4 import Dataset
from pyproj import CRS, Geod, Transformer
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage import measure


ZENODO_RECORD_ID = "15238132"
ZENODO_API_URL = f"https://zenodo.org/api/records/{ZENODO_RECORD_ID}"
ZENODO_RECORD_URL = f"https://zenodo.org/records/{ZENODO_RECORD_ID}"


def lon_to_180(lon: np.ndarray | xr.DataArray) -> np.ndarray | xr.DataArray:
    """Convert 0-360 longitude to -180-180 longitude."""
    return ((lon + 180.0) % 360.0) - 180.0


def zenodo_file_table(record_id: str = ZENODO_RECORD_ID) -> pd.DataFrame:
    """Return Zenodo file metadata as a DataFrame."""
    url = f"https://zenodo.org/api/records/{record_id}"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    record = response.json()
    rows = []
    for item in record.get("files", []):
        rows.append(
            {
                "key": item["key"],
                "size": int(item.get("size", 0)),
                "download_url": item["links"]["self"],
                "checksum": item.get("checksum", ""),
            }
        )
    return pd.DataFrame(rows).sort_values("key").reset_index(drop=True)


def download_zenodo_netcdfs(
    out_dir: str | Path = "data/raw",
    record_id: str = ZENODO_RECORD_ID,
    min_size: int = 1024,
) -> list[Path]:
    """Download NetCDF files from the Zenodo record, skipping tiny placeholders."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = zenodo_file_table(record_id)
    paths: list[Path] = []
    for _, row in files.iterrows():
        name = row["key"]
        size = int(row["size"])
        if not name.endswith(".nc") or size < min_size:
            continue
        path = out_dir / name
        if path.exists() and path.stat().st_size == size:
            paths.append(path)
            continue
        with requests.get(row["download_url"], stream=True, timeout=60) as response:
            response.raise_for_status()
            with path.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        paths.append(path)
    return paths


def _walk_netcdf_group(group: Dataset, source: Path, prefix: str = "") -> list[dict]:
    rows: list[dict] = []
    for name, variable in group.variables.items():
        rows.append(
            {
                "file": source.name,
                "group": prefix or "/",
                "variable": name,
                "dimensions": ", ".join(variable.dimensions),
                "shape": " x ".join(str(s) for s in variable.shape),
                "dtype": str(variable.dtype),
                "units": getattr(variable, "units", ""),
                "long_name": getattr(variable, "long_name", ""),
            }
        )
    for group_name, subgroup in group.groups.items():
        child_prefix = f"{prefix}/{group_name}" if prefix else f"/{group_name}"
        rows.extend(_walk_netcdf_group(subgroup, source, child_prefix))
    return rows


def netcdf_inventory(paths: str | Path | Iterable[str | Path]) -> pd.DataFrame:
    """Inspect variables, groups, dimensions, units, and long names in NetCDF files."""
    if isinstance(paths, (str, Path)):
        paths = [paths]
    rows: list[dict] = []
    for path in [Path(p) for p in paths]:
        with Dataset(path) as dataset:
            rows.extend(_walk_netcdf_group(dataset, path))
    return pd.DataFrame(rows)


def _netcdf_groups(path: str | Path) -> list[str]:
    with Dataset(path) as dataset:
        return list(dataset.groups.keys())


def open_swot(path: str | Path, swath: str = "left") -> xr.Dataset:
    """Open a SWOT NetCDF file and standardize common variable names."""
    path = Path(path)
    groups = _netcdf_groups(path)
    group = swath if swath in groups else None
    ds = xr.open_dataset(path, group=group, decode_times=False, mask_and_scale=True)

    rename: dict[str, str] = {}
    if "longitude" in ds:
        rename["longitude"] = "lon"
    if "latitude" in ds:
        rename["latitude"] = "lat"
    if "sea_ice_height" in ds:
        rename["sea_ice_height"] = "height"
        height_source = "sea_ice_height"
    elif "freeboard" in ds:
        rename["freeboard"] = "height"
        height_source = "freeboard"
    else:
        height_source = ""
    if "backscatter" in ds:
        rename["backscatter"] = "sigma0"
    if "ice_type_flag" in ds:
        rename["ice_type_flag"] = "ice_type"
    ds = ds.rename(rename)

    if "lon" in ds:
        ds["lon"] = lon_to_180(ds["lon"])
    coord_names = [name for name in ("lon", "lat") if name in ds]
    if coord_names:
        ds = ds.set_coords(coord_names)
    ds.attrs.update(
        {
            "source_file": path.name,
            "source_path": str(path),
            "swath": swath if group else "root",
            "height_source": height_source,
        }
    )
    return ds


def open_icesat2(path: str | Path, beam: str = "gt2l") -> pd.DataFrame:
    """Open one ICESat-2 beam and return a standardized point table."""
    path = Path(path)
    ds = xr.open_dataset(path, group=beam, decode_times=False, mask_and_scale=True)
    data = {
        "lon": lon_to_180(ds["longitude"].values),
        "lat": ds["latitude"].values,
        "datetime_seconds": ds["datetime"].values if "datetime" in ds else np.nan,
        "distance_m": ds["distance"].values if "distance" in ds else np.nan,
    }

    height_20m = None
    height_500m = None
    height_source = ""
    for candidate in ("sea_ice_height_20m", "freeboard_20m"):
        if candidate in ds:
            height_20m = ds[candidate].values
            height_source = candidate.replace("_20m", "")
            break
    for candidate in ("sea_ice_height_500m", "freeboard_500m"):
        if candidate in ds:
            height_500m = ds[candidate].values
            height_source = candidate.replace("_500m", "")
            break
    data["height_20m"] = height_20m if height_20m is not None else np.full(ds.sizes["x"], np.nan)
    data["height_500m"] = (
        height_500m if height_500m is not None else np.full(ds.sizes["x"], np.nan)
    )
    frame = pd.DataFrame(data)
    frame["beam"] = beam
    frame["source_file"] = path.name
    frame["height_source"] = height_source
    for col in ("distance_m", "height_20m", "height_500m"):
        frame.loc[np.abs(frame[col]) > 1.0e35, col] = np.nan
    valid = np.isfinite(frame["lon"]) & np.isfinite(frame["lat"])
    valid &= np.isfinite(frame["height_20m"]) | np.isfinite(frame["height_500m"])
    return frame.loc[valid].reset_index(drop=True)


def infer_pixel_spacing_m(lon: np.ndarray, lat: np.ndarray) -> dict[str, float]:
    """Infer median along/across grid spacing from neighboring lon/lat pairs."""
    lon = np.asarray(lon)
    lat = np.asarray(lat)
    geod = Geod(ellps="WGS84")

    def axis_spacing(axis: int) -> float:
        if axis == 0:
            lon1, lon2 = lon[:-1, :], lon[1:, :]
            lat1, lat2 = lat[:-1, :], lat[1:, :]
        else:
            lon1, lon2 = lon[:, :-1], lon[:, 1:]
            lat1, lat2 = lat[:, :-1], lat[:, 1:]
        valid = np.isfinite(lon1) & np.isfinite(lon2) & np.isfinite(lat1) & np.isfinite(lat2)
        if not valid.any():
            return np.nan
        _, _, distance = geod.inv(lon1[valid], lat1[valid], lon2[valid], lat2[valid])
        return float(np.nanmedian(distance))

    along = axis_spacing(0)
    across = axis_spacing(1)
    return {"along_m": along, "across_m": across, "median_m": float(np.nanmedian([along, across]))}


def provisional_ice_mask(
    height: np.ndarray,
    sigma0: np.ndarray | None = None,
    ice_type: np.ndarray | None = None,
    height_bounds: tuple[float, float] = (-1.0, 3.0),
) -> tuple[np.ndarray, dict[str, str]]:
    """Create a simple validity/ice mask for the exploratory prototype."""
    mask = np.isfinite(height)
    criteria = [f"finite height/freeboard", f"{height_bounds[0]} <= height/freeboard <= {height_bounds[1]} m"]
    mask &= (height >= height_bounds[0]) & (height <= height_bounds[1])
    if sigma0 is not None:
        mask &= np.isfinite(sigma0)
        criteria.append("finite sigma0/backscatter")
    if ice_type is not None:
        mask &= np.isfinite(ice_type) & (ice_type > 0)
        criteria.append("ice_type_flag > 0")
    return mask, {"criteria": "; ".join(criteria)}


def window_size_pixels(window_m: float, spacing: dict[str, float]) -> tuple[int, int]:
    """Convert a target window diameter in meters to along/across pixel counts."""
    along = int(max(2, round(window_m / spacing["along_m"])))
    across = int(max(2, round(window_m / spacing["across_m"])))
    return along, across


def local_nan_mean(values: np.ndarray, size: tuple[int, int], min_fraction: float = 0.5) -> np.ndarray:
    """Local mean that ignores NaNs and requires a minimum valid fraction."""
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    kernel = np.ones(size, dtype=float)
    count = ndimage.convolve(valid.astype(float), kernel, mode="constant", cval=0.0)
    total = ndimage.convolve(np.where(valid, values, 0.0), kernel, mode="constant", cval=0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = total / count
    mean[count < kernel.size * min_fraction] = np.nan
    return mean


def local_nan_std(values: np.ndarray, size: tuple[int, int], min_fraction: float = 0.5) -> np.ndarray:
    """Local standard deviation that ignores NaNs."""
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    kernel = np.ones(size, dtype=float)
    count = ndimage.convolve(valid.astype(float), kernel, mode="constant", cval=0.0)
    total = ndimage.convolve(np.where(valid, values, 0.0), kernel, mode="constant", cval=0.0)
    total2 = ndimage.convolve(np.where(valid, values**2, 0.0), kernel, mode="constant", cval=0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = total / count
        variance = total2 / count - mean**2
    std = np.sqrt(np.maximum(variance, 0.0))
    std[count < kernel.size * min_fraction] = np.nan
    return std


def fill_nan_nearest(values: np.ndarray) -> np.ndarray:
    """Fill NaNs with nearest finite neighbors for gradient calculations."""
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    if valid.all():
        return values.copy()
    if not valid.any():
        return values.copy()
    indices = ndimage.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    return values[tuple(indices)]


def gradient_magnitude(
    values: np.ndarray,
    spacing: dict[str, float],
    valid_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Gradient magnitude of height/freeboard in m per horizontal meter."""
    filled = fill_nan_nearest(values)
    grad_along, grad_across = np.gradient(filled, spacing["along_m"], spacing["across_m"])
    grad = np.hypot(grad_along, grad_across)
    grad[~np.isfinite(values)] = np.nan
    if valid_mask is not None:
        grad[~valid_mask] = np.nan
    return grad


def zscore(values: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """Standard z-score using finite values inside the analysis mask."""
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if mask is not None:
        finite &= mask
    mean = float(np.nanmean(values[finite]))
    std = float(np.nanstd(values[finite]))
    if not np.isfinite(std) or std == 0.0:
        return np.full(values.shape, np.nan)
    out = (values - mean) / std
    out[~np.isfinite(values)] = np.nan
    return out


def compute_morphology(
    swot: xr.Dataset,
    window_m: float = 1000.0,
    valid_mask: np.ndarray | None = None,
) -> xr.Dataset:
    """Compute local roughness, gradient, sigma0 texture, and SMI variants."""
    height = swot["height"].values.astype(float)
    sigma0 = swot["sigma0"].values.astype(float) if "sigma0" in swot else None
    if valid_mask is None:
        ice_type = swot["ice_type"].values if "ice_type" in swot else None
        valid_mask, mask_info = provisional_ice_mask(height, sigma0=sigma0, ice_type=ice_type)
    else:
        mask_info = {"criteria": "user-supplied valid mask"}

    spacing = infer_pixel_spacing_m(swot["lon"].values, swot["lat"].values)
    size = window_size_pixels(window_m, spacing)
    height_for_metrics = np.where(valid_mask, height, np.nan)
    sigma_for_metrics = np.where(valid_mask, sigma0, np.nan) if sigma0 is not None else np.full_like(height, np.nan)

    height_std = local_nan_std(height_for_metrics, size=size)
    height_grad = gradient_magnitude(height_for_metrics, spacing=spacing, valid_mask=valid_mask)
    sigma0_std = local_nan_std(sigma_for_metrics, size=size)
    height_roughness_z = zscore(height_std, valid_mask)
    height_gradient_z = zscore(height_grad, valid_mask)
    sigma0_texture_z = zscore(sigma0_std, valid_mask)

    height_only_smi = height_roughness_z + height_gradient_z
    backscatter_only_smi = sigma0_texture_z.copy()
    stability_morphology_index = height_only_smi + backscatter_only_smi
    for index in (height_only_smi, backscatter_only_smi, stability_morphology_index):
        index[~valid_mask] = np.nan

    out = xr.Dataset(
        data_vars={
            "valid_mask": (swot["height"].dims, valid_mask),
            "height_std": (swot["height"].dims, height_std),
            "height_gradient": (swot["height"].dims, height_grad),
            "sigma0_std": (swot["height"].dims, sigma0_std),
            "height_roughness_z": (swot["height"].dims, height_roughness_z),
            "height_gradient_z": (swot["height"].dims, height_gradient_z),
            "sigma0_texture_z": (swot["height"].dims, sigma0_texture_z),
            "height_only_smi": (swot["height"].dims, height_only_smi),
            "backscatter_only_smi": (swot["height"].dims, backscatter_only_smi),
            "stability_morphology_index": (swot["height"].dims, stability_morphology_index),
        },
        coords={name: swot[name] for name in ("lon", "lat") if name in swot},
        attrs={
            "window_m": float(window_m),
            "window_pixels_along": int(size[0]),
            "window_pixels_across": int(size[1]),
            "spacing_along_m": spacing["along_m"],
            "spacing_across_m": spacing["across_m"],
            "mask_criteria": mask_info["criteria"],
            "smi_definition": "z(local height/freeboard std) + z(height/freeboard gradient magnitude) + z(local sigma0 std)",
            "height_only_smi_definition": "z(local height/freeboard std) + z(height/freeboard gradient magnitude)",
            "backscatter_only_smi_definition": "z(local sigma0 std)",
        },
    )
    return out


def index_mask(
    index: np.ndarray,
    valid_mask: np.ndarray | None = None,
    quantile: float | None = 0.90,
    sigma_threshold: float | None = None,
) -> tuple[np.ndarray, float, str]:
    """Threshold a stability-relevant morphology index."""
    finite = np.isfinite(index)
    if valid_mask is not None:
        finite &= valid_mask
    if quantile is not None:
        threshold = float(np.nanquantile(index[finite], quantile))
        label = f"top {int(round((1.0 - quantile) * 100))}%"
    elif sigma_threshold is not None:
        mean = float(np.nanmean(index[finite]))
        std = float(np.nanstd(index[finite]))
        threshold = mean + sigma_threshold * std
        label = f">{sigma_threshold:g} SD above background"
    else:
        raise ValueError("Set either quantile or sigma_threshold.")
    mask = finite & (index >= threshold)
    return mask, threshold, label


def candidate_mask(
    score: np.ndarray,
    valid_mask: np.ndarray | None = None,
    quantile: float | None = 0.90,
    sigma_threshold: float | None = None,
) -> tuple[np.ndarray, float, str]:
    """Backward-compatible wrapper for thresholding morphology-index arrays."""
    return index_mask(score, valid_mask=valid_mask, quantile=quantile, sigma_threshold=sigma_threshold)


def label_index_features(
    mask: np.ndarray,
    spacing: dict[str, float],
    min_pixels: int = 4,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Label connected high-index morphology features and summarize their area."""
    labels = measure.label(mask.astype(bool), connectivity=2)
    rows = []
    pixel_area_km2 = spacing["along_m"] * spacing["across_m"] / 1.0e6
    for label_id in range(1, labels.max() + 1):
        count = int(np.sum(labels == label_id))
        if count < min_pixels:
            labels[labels == label_id] = 0
            continue
        rows.append(
            {
                "feature_id": label_id,
                "pixel_count": count,
                "area_km2": count * pixel_area_km2,
            }
        )
    relabeled = measure.label(labels > 0, connectivity=2)
    if rows:
        table = pd.DataFrame(rows)
        table["feature_id"] = np.arange(1, len(table) + 1)
    else:
        table = pd.DataFrame(columns=["feature_id", "pixel_count", "area_km2"])
    return relabeled, table


def label_candidate_features(
    mask: np.ndarray,
    spacing: dict[str, float],
    min_pixels: int = 4,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Backward-compatible wrapper for labeling high-index morphology features."""
    return label_index_features(mask, spacing=spacing, min_pixels=min_pixels)


def local_projector(lon: np.ndarray, lat: np.ndarray) -> Transformer:
    """Create a local azimuthal equidistant lon/lat-to-meter transformer."""
    lon0 = float(np.nanmedian(lon))
    lat0 = float(np.nanmedian(lat))
    crs = CRS.from_proj4(f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +datum=WGS84 +units=m +no_defs")
    return Transformer.from_crs("EPSG:4326", crs, always_xy=True)


def sample_grid_nearest(
    lon_grid: np.ndarray,
    lat_grid: np.ndarray,
    value_grid: np.ndarray,
    point_lon: np.ndarray,
    point_lat: np.ndarray,
    max_distance_m: float = 750.0,
    valid_mask: np.ndarray | None = None,
) -> pd.DataFrame:
    """Nearest-neighbor sample a SWOT grid at ICESat-2 point locations."""
    finite = np.isfinite(lon_grid) & np.isfinite(lat_grid) & np.isfinite(value_grid)
    if valid_mask is not None:
        finite &= valid_mask
    if not finite.any():
        raise ValueError("No finite grid cells available for sampling.")

    transformer = local_projector(lon_grid[finite], lat_grid[finite])
    grid_x, grid_y = transformer.transform(lon_grid[finite], lat_grid[finite])
    point_x, point_y = transformer.transform(point_lon, point_lat)
    ij = np.column_stack(np.where(finite))
    values = value_grid[finite]
    tree = cKDTree(np.column_stack([grid_x, grid_y]))
    distance, index = tree.query(np.column_stack([point_x, point_y]), k=1)

    sampled = np.full(point_lon.shape, np.nan, dtype=float)
    grid_i = np.full(point_lon.shape, -1, dtype=int)
    grid_j = np.full(point_lon.shape, -1, dtype=int)
    keep = distance <= max_distance_m
    sampled[keep] = values[index[keep]]
    grid_i[keep] = ij[index[keep], 0]
    grid_j[keep] = ij[index[keep], 1]
    return pd.DataFrame(
        {
            "sampled_value": sampled,
            "sample_distance_m": distance,
            "grid_i": grid_i,
            "grid_j": grid_j,
        }
    )


def choose_cross_section_segment(
    points: pd.DataFrame,
    score_col: str = "swot_smi",
    half_width_m: float = 15_000.0,
    min_points: int = 500,
) -> tuple[pd.DataFrame, float]:
    """Select an ICESat-2 segment centered on the highest sampled SWOT SMI."""
    finite = points[np.isfinite(points[score_col]) & np.isfinite(points["distance_m"])].copy()
    if finite.empty:
        return finite, np.nan
    center_idx = finite[score_col].idxmax()
    center_distance = float(finite.loc[center_idx, "distance_m"])
    segment = finite[
        (finite["distance_m"] >= center_distance - half_width_m)
        & (finite["distance_m"] <= center_distance + half_width_m)
    ].copy()
    if len(segment) < min_points:
        segment = finite.copy()
    return segment.reset_index(drop=True), center_distance


def sensitivity_table(
    metrics_by_window: dict[float, xr.Dataset],
    thresholds: Iterable[float | str] = (0.90, 0.95, "mean+2sd"),
    min_pixels: int = 4,
    index_var: str = "stability_morphology_index",
) -> pd.DataFrame:
    """Summarize high-index area and connected-component count across settings."""
    rows = []
    for window_m, metrics in metrics_by_window.items():
        valid = metrics["valid_mask"].values.astype(bool)
        index = metrics[index_var].values
        spacing = {
            "along_m": float(metrics.attrs["spacing_along_m"]),
            "across_m": float(metrics.attrs["spacing_across_m"]),
        }
        for threshold in thresholds:
            if isinstance(threshold, str):
                mask, threshold_value, label = index_mask(
                    index, valid_mask=valid, quantile=None, sigma_threshold=2.0
                )
            else:
                mask, threshold_value, label = index_mask(index, valid_mask=valid, quantile=float(threshold))
            _, features = label_index_features(mask, spacing=spacing, min_pixels=min_pixels)
            pixel_area_km2 = spacing["along_m"] * spacing["across_m"] / 1.0e6
            rows.append(
                {
                    "index_variant": index_var,
                    "window_m": float(window_m),
                    "window_pixels_along": int(metrics.attrs["window_pixels_along"]),
                    "window_pixels_across": int(metrics.attrs["window_pixels_across"]),
                    "threshold": label,
                    "threshold_value": threshold_value,
                    "high_index_pixels": int(mask.sum()),
                    "high_index_area_km2": float(mask.sum() * pixel_area_km2),
                    "feature_count": int(len(features)),
                }
            )
    return pd.DataFrame(rows)


def component_sensitivity_table(
    metrics: xr.Dataset,
    index_vars: Iterable[str] = (
        "height_only_smi",
        "backscatter_only_smi",
        "stability_morphology_index",
    ),
    thresholds: Iterable[float | str] = (0.90, 0.95),
    min_pixels: int = 4,
) -> pd.DataFrame:
    """Summarize sensitivity to height-only, backscatter-only, and combined SMI."""
    rows = []
    valid = metrics["valid_mask"].values.astype(bool)
    spacing = {
        "along_m": float(metrics.attrs["spacing_along_m"]),
        "across_m": float(metrics.attrs["spacing_across_m"]),
    }
    pixel_area_km2 = spacing["along_m"] * spacing["across_m"] / 1.0e6
    for index_var in index_vars:
        index = metrics[index_var].values
        for threshold in thresholds:
            if isinstance(threshold, str):
                mask, threshold_value, label = index_mask(
                    index, valid_mask=valid, quantile=None, sigma_threshold=2.0
                )
            else:
                mask, threshold_value, label = index_mask(index, valid_mask=valid, quantile=float(threshold))
            _, features = label_index_features(mask, spacing=spacing, min_pixels=min_pixels)
            rows.append(
                {
                    "index_variant": index_var,
                    "threshold": label,
                    "threshold_value": threshold_value,
                    "high_index_pixels": int(mask.sum()),
                    "high_index_area_km2": float(mask.sum() * pixel_area_km2),
                    "feature_count": int(len(features)),
                }
            )
    return pd.DataFrame(rows)


def smi_variogram(
    smi: np.ndarray,
    valid: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    max_distance_m: float = 60_000,
    n_bins: int = 24,
    n_sample: int = 6000,
    seed: int = 42,
) -> pd.DataFrame:
    """Empirical semi-variogram of the SMI in projected metric coordinates.

    Returns a DataFrame with columns: distance_m, semivariance, n_pairs.
    Subsamples up to n_sample valid pixels to keep runtime tractable.
    """
    rng = np.random.default_rng(seed)
    lon_v = lon[valid].ravel()
    lat_v = lat[valid].ravel()
    val_v = smi[valid].ravel()
    finite = np.isfinite(val_v)
    lon_v, lat_v, val_v = lon_v[finite], lat_v[finite], val_v[finite]

    if len(lon_v) > n_sample:
        idx = rng.choice(len(lon_v), n_sample, replace=False)
        lon_v, lat_v, val_v = lon_v[idx], lat_v[idx], val_v[idx]

    transformer = local_projector(lon_v, lat_v)
    x, y = transformer.transform(lon_v, lat_v)
    pts = np.column_stack([x, y])

    from scipy.spatial.distance import pdist
    dists = pdist(pts)
    i_idx, j_idx = np.triu_indices(len(val_v), k=1)
    sq_diff = (val_v[i_idx] - val_v[j_idx]) ** 2

    keep = dists <= max_distance_m
    dists, sq_diff = dists[keep], sq_diff[keep]

    bin_edges = np.linspace(0, max_distance_m, n_bins + 1)
    bin_idx = np.digitize(dists, bin_edges) - 1
    rows = []
    for b in range(n_bins):
        mask = bin_idx == b
        if mask.sum() < 10:
            continue
        rows.append({
            "distance_m": float((bin_edges[b] + bin_edges[b + 1]) / 2),
            "semivariance": float(np.mean(sq_diff[mask]) / 2),
            "n_pairs": int(mask.sum()),
        })
    return pd.DataFrame(rows)


def feature_shape_metrics(
    feature_labels: np.ndarray,
    spacing: dict[str, float],
    min_pixels: int = 10,
) -> pd.DataFrame:
    """Orientation, eccentricity, and elongation for each labelled high-index feature.

    Uses skimage regionprops on the label image. Orientation is in degrees clockwise
    from the across-track (column) axis.
    """
    from skimage.measure import regionprops

    pixel_area_km2 = spacing["along_m"] * spacing["across_m"] / 1e6
    rows = []
    for prop in regionprops(feature_labels):
        if prop.num_pixels < min_pixels:
            continue
        # prop.orientation is in radians, counter-clockwise from column axis
        angle_deg = float(np.degrees(prop.orientation))
        rows.append({
            "feature_id": int(prop.label),
            "area_km2": float(prop.num_pixels * pixel_area_km2),
            "eccentricity": float(prop.eccentricity),
            "orientation_deg": angle_deg,
            "axis_major_km": float(prop.axis_major_length * spacing["along_m"] / 1000),
            "axis_minor_km": float(prop.axis_minor_length * spacing["along_m"] / 1000),
            "elongation": float(
                prop.axis_major_length / prop.axis_minor_length
                if prop.axis_minor_length > 0 else np.nan
            ),
        })
    return pd.DataFrame(rows)


def download_bathymetry(
    lon_min: float,
    lon_max: float,
    lat_min: float,
    lat_max: float,
    cache_dir: str | Path = "data/raw",
    resolution_deg: float = 0.01,
) -> xr.Dataset:
    """Download GEBCO 2023 bathymetry for a bounding box, cached as a NetCDF.

    Falls back to ETOPO 2022 via NOAA ERDDAP if GEBCO is unavailable.
    Elevation values are in metres; negative = below sea level.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{lat_min:.2f}_{lat_max:.2f}_{lon_min:.2f}_{lon_max:.2f}"
    cache_path = cache_dir / f"gebco_bathy_{tag}.nc"
    if cache_path.exists():
        return xr.open_dataset(cache_path)

    # ---- primary: GEBCO 2023 REST download API ----
    try:
        url = "https://download.gebco.net/"
        params = {
            "format": "netCDF4_compressed",
            "north": lat_max,
            "south": lat_min,
            "west": lon_min,
            "east": lon_max,
        }
        r = requests.get(url, params=params, timeout=180)
        r.raise_for_status()
        import io
        ds = xr.open_dataset(io.BytesIO(r.content))
        ds.to_netcdf(cache_path)
        return ds
    except Exception:
        pass

    # ---- fallback: ETOPO 2022 via NOAA ERDDAP ----
    step = max(1, int(round(resolution_deg / (1 / 60))))  # convert deg to arcminute steps
    erddap = (
        "https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo360.nc"
        f"?altitude[({lat_min}):{step}:({lat_max})][({lon_min % 360}):{step}:({lon_max % 360})]"
    )
    r = requests.get(erddap, timeout=180)
    r.raise_for_status()
    import io
    ds = xr.open_dataset(io.BytesIO(r.content))
    # Normalise to 'elevation' variable name and -180/180 longitude.
    if "altitude" in ds:
        ds = ds.rename({"altitude": "elevation"})
    if "longitude" in ds:
        ds["longitude"] = lon_to_180(ds["longitude"])
        ds = ds.rename({"longitude": "lon", "latitude": "lat"})
    ds.to_netcdf(cache_path)
    return ds


def smi_by_depth_bin(
    smi: np.ndarray,
    depth: np.ndarray,
    valid: np.ndarray,
    high_index: np.ndarray,
    bins: list[float] | None = None,
) -> pd.DataFrame:
    """Summarise SMI statistics and high-index fraction by water depth bin.

    Depth should be in metres (negative = below sea level).
    """
    if bins is None:
        bins = [0, -5, -10, -20, -40, -100, -np.inf]
    smi_flat = smi[valid].ravel()
    depth_flat = depth[valid].ravel()
    hi_flat = high_index[valid].ravel()
    finite = np.isfinite(smi_flat) & np.isfinite(depth_flat)
    smi_flat, depth_flat, hi_flat = smi_flat[finite], depth_flat[finite], hi_flat[finite]

    rows = []
    for lo, hi_bound in zip(bins[:-1], bins[1:]):
        # lo >= hi_bound since depth is negative
        mask = (depth_flat <= lo) & (depth_flat > hi_bound)
        if mask.sum() == 0:
            continue
        rows.append({
            "depth_bin": f"{int(hi_bound)} to {int(lo)} m" if np.isfinite(hi_bound) else f"< {int(lo)} m",
            "n_pixels": int(mask.sum()),
            "mean_smi": float(np.nanmean(smi_flat[mask])),
            "p90_smi": float(np.nanpercentile(smi_flat[mask], 90)),
            "high_index_fraction": float(hi_flat[mask].mean()),
        })
    return pd.DataFrame(rows)


def make_common_grid(
    lon_arrays: list[np.ndarray],
    lat_arrays: list[np.ndarray],
    valid_arrays: list[np.ndarray],
    resolution_deg: float = 0.01,
) -> tuple[np.ndarray, np.ndarray]:
    """Regular lon/lat grid covering the union of valid pixels from all swaths."""
    all_lon = np.concatenate([lon[v] for lon, v in zip(lon_arrays, valid_arrays)])
    all_lat = np.concatenate([lat[v] for lat, v in zip(lat_arrays, valid_arrays)])
    lon_centers = np.arange(np.nanmin(all_lon), np.nanmax(all_lon) + resolution_deg, resolution_deg)
    lat_centers = np.arange(np.nanmin(all_lat), np.nanmax(all_lat) + resolution_deg, resolution_deg)
    return lon_centers, lat_centers


def rasterize_smi(
    lon: np.ndarray,
    lat: np.ndarray,
    smi: np.ndarray,
    high_index: np.ndarray,
    valid: np.ndarray,
    lon_centers: np.ndarray,
    lat_centers: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Bin SMI and high-index mask onto a regular grid using mean aggregation.

    Returns (smi_grid, hi_fraction_grid) both shaped (n_lat, n_lon).
    NaN where no valid native pixels fell in that bin.
    """
    from scipy.stats import binned_statistic_2d

    res = float(lon_centers[1] - lon_centers[0]) if len(lon_centers) > 1 else 0.01
    lon_edges = np.append(lon_centers - res / 2, lon_centers[-1] + res / 2)
    lat_edges = np.append(lat_centers - res / 2, lat_centers[-1] + res / 2)

    lon_flat = lon[valid]
    lat_flat = lat[valid]
    smi_flat = np.where(np.isfinite(smi[valid]), smi[valid], np.nan)
    hi_flat = high_index[valid].astype(float)

    def _bin(values: np.ndarray) -> np.ndarray:
        finite = np.isfinite(values)
        if not finite.any():
            return np.full((len(lat_centers), len(lon_centers)), np.nan)
        result, _, _, _ = binned_statistic_2d(
            lon_flat[finite], lat_flat[finite], values[finite],
            statistic="mean", bins=[lon_edges, lat_edges],
        )
        return result.T  # (lat, lon)

    smi_grid = _bin(smi_flat)
    hi_grid = _bin(hi_flat)
    return smi_grid, hi_grid


def persistence_count(
    hi_grids: list[np.ndarray],
    hi_threshold: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Count how many dates each grid cell is high-index.

    Parameters
    ----------
    hi_grids:
        List of hi_fraction arrays from rasterize_smi, float with NaN where no data.
    hi_threshold:
        Mean fraction cutoff to call a cell high-index for a given date.

    Returns
    -------
    count_grid : float array, NaN where no date has coverage.
    coverage_grid : int array, number of dates with data in each cell.
    """
    shape = hi_grids[0].shape
    count = np.zeros(shape, dtype=float)
    coverage = np.zeros(shape, dtype=int)
    for g in hi_grids:
        has_data = np.isfinite(g)
        coverage += has_data.astype(int)
        count += np.where(has_data, (g >= hi_threshold).astype(float), 0.0)
    count_grid = np.where(coverage > 0, count, np.nan)
    return count_grid, coverage


def jaccard_similarity(
    hi_a: np.ndarray,
    hi_b: np.ndarray,
    hi_threshold: float = 0.5,
) -> float:
    """Jaccard similarity between two hi_fraction grids over their shared coverage."""
    both = np.isfinite(hi_a) & np.isfinite(hi_b)
    if not both.any():
        return np.nan
    a = hi_a[both] >= hi_threshold
    b = hi_b[both] >= hi_threshold
    union = int((a | b).sum())
    return float((a & b).sum() / union) if union > 0 else np.nan


def swath_axis_artifact_summary(index: np.ndarray, valid_mask: np.ndarray) -> pd.DataFrame:
    """Simple diagnostic for row/column index structure that could indicate striping."""
    index = np.asarray(index, dtype=float)
    masked = np.where(valid_mask, index, np.nan)
    row_median = np.nanmedian(masked, axis=1)
    col_median = np.nanmedian(masked, axis=0)
    return pd.DataFrame(
        [
            {
                "axis": "along-track rows",
                "median_of_axis_medians": float(np.nanmedian(row_median)),
                "std_of_axis_medians": float(np.nanstd(row_median)),
                "p95_abs_axis_median": float(np.nanpercentile(np.abs(row_median), 95)),
            },
            {
                "axis": "across-track columns",
                "median_of_axis_medians": float(np.nanmedian(col_median)),
                "std_of_axis_medians": float(np.nanstd(col_median)),
                "p95_abs_axis_median": float(np.nanpercentile(np.abs(col_median), 95)),
            },
        ]
    )
