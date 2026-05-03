from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

import swot_landfast_utils as u


def _mesh(shape=(18, 12)):
    row = np.arange(shape[0])[:, None]
    col = np.arange(shape[1])[None, :]
    lon = 220.0 + col * 0.01 + row * 0.002
    lat = 73.0 + row * 0.01 + col * 0.001
    height = 0.35 + 0.06 * np.sin(row / 2.0) + 0.04 * np.cos(col / 3.0)
    sigma0 = 5.0 + 1.2 * np.cos(row / 3.0) + 0.7 * np.sin(col / 2.0)
    return lon.astype("f8"), lat.astype("f8"), height.astype("f8"), sigma0.astype("f8")


def _write_grouped_file(
    path: Path,
    group_name: str,
    lon_name: str,
    lat_name: str,
    height_name: str,
    sigma0_name: str,
):
    lon, lat, height, sigma0 = _mesh()
    with Dataset(path, "w") as nc:
        nc.createDimension("x_along", lon.shape[0])
        nc.createDimension("y_across", lon.shape[1])
        group = nc.createGroup(group_name)
        for name, values in (
            (lon_name, lon),
            (lat_name, lat),
            (height_name, height),
            (sigma0_name, sigma0),
        ):
            var = group.createVariable(name, "f8", ("x_along", "y_across"))
            var[:] = values


def test_open_swot_standardizes_umd_zenodo_group(tmp_path):
    path = tmp_path / "umd_swot.nc"
    _write_grouped_file(
        path,
        group_name="left",
        lon_name="longitude",
        lat_name="latitude",
        height_name="sea_ice_height",
        sigma0_name="backscatter",
    )

    ds = u.open_swot(path, swath="left")

    assert {"lon", "lat", "height", "sigma0"}.issubset(ds.variables)
    assert ds.attrs["product_format"] == "umd_zenodo"
    assert ds.attrs["height_source"] == "sea_ice_height"
    assert ds.attrs["sigma0_source"] == "backscatter"
    assert float(ds["lon"].max()) < 0


def test_open_swot_standardizes_raw_l2_karin_group(tmp_path):
    path = tmp_path / "swot_l2_karin.nc"
    _write_grouped_file(
        path,
        group_name="right",
        lon_name="longitude",
        lat_name="latitude",
        height_name="ssh_karin",
        sigma0_name="sig0_karin",
    )

    ds = u.open_swot(path, swath="right")

    assert {"lon", "lat", "height", "sigma0"}.issubset(ds.variables)
    assert ds.attrs["product_format"] == "swot_l2_karin"
    assert ds.attrs["height_source"] == "ssh_karin"
    assert ds.attrs["sigma0_source"] == "sig0_karin"
    np.testing.assert_allclose(ds["height"].values[0, 0], 0.39)


def test_open_swot_reports_missing_standard_variables(tmp_path):
    path = tmp_path / "bad.nc"
    with Dataset(path, "w") as nc:
        nc.createDimension("x", 2)
        nc.createVariable("longitude", "f8", ("x",))[:] = [220.0, 220.1]

    with pytest.raises(ValueError, match="missing"):
        u.open_swot(path)


def test_compute_morphology_and_sensitivity_tables(tmp_path):
    path = tmp_path / "swot_l2_karin.nc"
    _write_grouped_file(
        path,
        group_name="left",
        lon_name="longitude",
        lat_name="latitude",
        height_name="ssh_karin",
        sigma0_name="sig0_karin",
    )
    swot = u.open_swot(path, swath="left")
    valid, _ = u.provisional_ice_mask(swot["height"].values, sigma0=swot["sigma0"].values)
    metrics = u.compute_morphology(swot, window_m=1000, valid_mask=valid)

    assert "stability_morphology_index" in metrics
    assert "height_only_smi" in metrics
    assert "backscatter_only_smi" in metrics
    assert np.isfinite(metrics["stability_morphology_index"].values[valid]).any()

    high, threshold, label = u.index_mask(
        metrics["stability_morphology_index"].values,
        valid_mask=valid,
        quantile=0.90,
    )
    assert label == "top 10%"
    assert np.isfinite(threshold)
    assert high.sum() > 0

    sensitivity = u.sensitivity_table(
        {500: metrics, 1000: metrics},
        thresholds=(0.90, "mean+2sd"),
    )
    assert {"high_index_pixels", "high_index_area_km2", "feature_count"}.issubset(sensitivity.columns)
    assert len(sensitivity) == 4

    component = u.component_sensitivity_table(metrics, thresholds=(0.90,))
    assert set(component["index_variant"]) == {
        "height_only_smi",
        "backscatter_only_smi",
        "stability_morphology_index",
    }


def test_sampling_labels_and_persistence_helpers():
    lon, lat, _, _ = _mesh(shape=(6, 5))
    values = np.arange(lon.size, dtype=float).reshape(lon.shape)
    points = u.sample_grid_nearest(
        lon,
        lat,
        values,
        np.array([lon[2, 3], lon[4, 1]]),
        np.array([lat[2, 3], lat[4, 1]]),
        max_distance_m=10,
    )

    np.testing.assert_allclose(points["sampled_value"].values, [values[2, 3], values[4, 1]])

    mask = np.zeros((6, 5), dtype=bool)
    mask[1:3, 1:3] = True
    mask[4, 4] = True
    labels, features = u.label_index_features(mask, {"along_m": 250.0, "across_m": 250.0}, min_pixels=2)
    assert labels.max() == 1
    assert len(features) == 1

    count, coverage = u.persistence_count(
        [mask.astype(float), np.where(mask, 1.0, np.nan)],
        hi_threshold=0.5,
    )
    assert np.nanmax(count) == 2
    assert np.nanmax(coverage) == 2
