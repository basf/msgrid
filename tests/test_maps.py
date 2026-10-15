import numpy as np
import pytest

from msgrid import Structure, compute_mep, compute_mlp, run_maps
from msgrid.mlp import METHODS, assign_fi, get_table, load_logp_file


def _grid_points(gm):
    axes = [gm.origin[a] + np.arange(gm.data.shape[a]) * gm.spacing for a in range(3)]
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)


def test_read_pqr_group_charges(pqr_file):
    s = Structure.from_file(pqr_file)
    assert s.has_charges
    pos, q = s.heavy_group_charges()
    assert len(pos) == 14  # 17 atoms - 3 H
    # backbone N of ALA 1 carries its bonded H (1.0 A away)
    assert q[0] == pytest.approx(-0.4157 + 0.2719)
    # CB carries HB1
    assert q[2] == pytest.approx(-0.1825 + 0.0603)


def test_read_pdb_has_no_charges(pdb_file):
    s = Structure.from_file(pdb_file)
    assert not s.has_charges
    assert s.is_hydrogen.sum() == 3
    with pytest.raises(ValueError):
        s.heavy_group_charges()


@pytest.mark.parametrize("method", list(METHODS))
def test_mlp_matches_bruteforce(pdb_file, method):
    s = Structure.from_file(pdb_file)
    xyz, fi = assign_fi(s, get_table("Crippen"))
    gm = compute_mlp(s, method=method, spacing=2.0, cutoff=None)
    pts = _grid_points(gm)
    d = np.linalg.norm(pts[:, None, :] - xyz[None, :, :], axis=2)
    with np.errstate(divide="ignore"):
        expected = (100.0 * fi * METHODS[method](d)).sum(axis=1)
    mask = np.isfinite(expected)
    np.testing.assert_allclose(gm.data.ravel()[mask], expected[mask], rtol=1e-9, atol=1e-9)


def test_mlp_fauchere_cutoff_is_negligible(pdb_file):
    s = Structure.from_file(pdb_file)
    full = compute_mlp(s, spacing=2.0, cutoff=None).data
    cut = compute_mlp(s, spacing=2.0).data
    assert np.abs(full - cut).max() < 1e-6


def test_mlp_residue_alias(pdb_file, tmp_path):
    text = pdb_file.read_text().replace("ASP", "ASH")
    p = tmp_path / "ash.pdb"
    p.write_text(text)
    compute_mlp(Structure.from_file(p), spacing=3.0)


def test_mlp_missing_atom_is_reported(pdb_file, tmp_path):
    p = tmp_path / "bad.pdb"
    p.write_text(pdb_file.read_text().replace(" CB ", " CX "))
    with pytest.raises(KeyError, match="CX"):
        compute_mlp(Structure.from_file(p))


def test_mep_matches_bruteforce(pqr_file):
    s = Structure.from_file(pqr_file)
    pos, q = s.heavy_group_charges()
    gm = compute_mep(s, spacing=1.5, padding=3.0, cutoff_factor=1e3)
    pts = _grid_points(gm)
    r = np.maximum(np.linalg.norm(pts[:, None, :] - pos[None, :, :], axis=2), 0.5)
    expected = 1200.0 * (q * np.exp(-r / 8.0) / (4.0 * r)).sum(axis=1)
    np.testing.assert_allclose(gm.data.ravel(), expected, rtol=1e-4, atol=1e-3)


def test_threads_give_same_result(pqr_file):
    s = Structure.from_file(pqr_file)
    a = compute_mep(s, spacing=1.0, n_jobs=1).data
    b = compute_mep(s, spacing=1.0, n_jobs=3).data
    np.testing.assert_array_equal(a, b)


def test_dx_roundtrip(pqr_file, tmp_path):
    griddata = pytest.importorskip("gridData")
    gm = compute_mep(Structure.from_file(pqr_file), spacing=1.5)
    out = tmp_path / "mep.dx"
    gm.write_dx(out)
    g = griddata.Grid(str(out))
    assert g.grid.shape == gm.data.shape
    np.testing.assert_allclose(g.origin, gm.origin, atol=1e-5)
    np.testing.assert_allclose(g.grid, gm.data, rtol=1e-5, atol=1e-4)


def test_run_maps_both(pqr_file):
    out = run_maps(pqr_file, mlp_options={"spacing": 2.0}, mep_options={"spacing": 1.5})
    assert set(out) == {"MLP", "MEP"}
    assert all(p.exists() and p.stat().st_size > 0 for p in out.values())
    assert out["MLP"].name == "mol_MLP.dx"


def test_run_maps_pdb_needs_no_mep(pdb_file):
    with pytest.raises(ValueError, match="charges"):
        run_maps(pdb_file)
    out = run_maps(pdb_file, mep=False, mlp_options={"spacing": 2.0})
    assert set(out) == {"MLP"}


def test_load_logp_file_does_not_execute(tmp_path):
    marker = tmp_path / "pwned"
    f = tmp_path / "logp.py"
    f.write_text(
        f"open({str(marker)!r}, 'w').write('x')\nCrippen = {{'ALA': {{'CA': -0.2}}}}\n"
    )
    assert load_logp_file(f) == {"ALA": {"CA": -0.2}}
    assert not marker.exists()


@pytest.fixture
def fake_pdb2pqr(tmp_path, pqr_file):
    import shutil
    import sys

    template = tmp_path / "template.pqr.txt"
    shutil.copy(pqr_file, template)
    exe = tmp_path / "fake_pdb2pqr"
    exe.write_text(
        f"#!{sys.executable}\n"
        "import shutil, sys\n"
        f"open({str(tmp_path / 'args.txt')!r}, 'w').write(' '.join(sys.argv[1:-2]))\n"
        f"shutil.copy({str(template)!r}, sys.argv[-1])\n"
    )
    exe.chmod(0o755)
    return exe


def test_run_maps_pdb2pqr_conversion(pdb_file, fake_pdb2pqr, tmp_path):
    out = run_maps(
        pdb_file,
        mlp_options={"spacing": 2.0},
        mep_options={"spacing": 1.5},
        pdb2pqr_options={"executable": str(fake_pdb2pqr), "ph": 6.5},
    )
    assert set(out) == {"PQR", "MLP", "MEP"}
    assert out["PQR"].name == "mol.pqr"
    args = (tmp_path / "args.txt").read_text()
    assert "--ff=AMBER" in args and "--with-ph=6.5" in args and "--titration-state-method=propka" in args


def test_pdb2pqr_failure_is_reported(pdb_file, tmp_path):
    bad = tmp_path / "bad_pdb2pqr"
    bad.write_text("#!/bin/sh\necho boom >&2\nexit 3\n")
    bad.chmod(0o755)
    with pytest.raises(RuntimeError, match="boom"):
        run_maps(pdb_file, pdb2pqr_options={"executable": str(bad)})


def test_pdb2pqr_missing(pdb_file, monkeypatch):
    monkeypatch.setattr("shutil.which", lambda _name: None)
    with pytest.raises(RuntimeError, match="not found"):
        run_maps(pdb_file, pdb2pqr_options={})


def test_pdb_altlocs_keep_first(tmp_path):
    pdb = tmp_path / "alt.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA A   1      11.000  12.000  13.000  1.00  0.00           N\n"
        "ATOM      2  CA AALA A   1      12.000  12.000  13.000  0.50  0.00           C\n"
        "ATOM      3  CA BALA A   1      12.100  12.000  13.000  0.50  0.00           C\n"
        "ATOM      4  C   ALA A   1      13.000  12.000  13.000  1.00  0.00           C\n"
    )
    s = Structure.from_file(pdb)
    assert len(s) == 3
    assert s.xyz[1, 0] == pytest.approx(12.0)


def test_pdb_bad_record_reports_line(tmp_path):
    pdb = tmp_path / "bad.pdb"
    pdb.write_text("ATOM      1  N   ALA A   1      xxxxxx  12.000  13.000  1.00  0.00           N\n")
    with pytest.raises(ValueError, match=r"bad\.pdb:1"):
        Structure.from_file(pdb)


def test_pqr_run_together_coordinates(tmp_path):
    pqr = tmp_path / "merged.pqr"
    pqr.write_text(
        "ATOM      1  N   ALA A   1    -123.456-234.567-345.678 -0.4157 1.8240\n"
    )
    s = Structure.from_file(pqr)
    np.testing.assert_allclose(s.xyz[0], [-123.456, -234.567, -345.678])
    assert s.charge[0] == pytest.approx(-0.4157)
    assert s.radius[0] == pytest.approx(1.824)


def test_pqr_metal_ion_is_not_hydrogen(tmp_path):
    pqr = tmp_path / "ion.pqr"
    pqr.write_text(
        "ATOM      1  N   ALA A   1       1.000   2.000   3.000 -0.4157 1.8240\n"
        "HETATM    2 HG   HG  A   2       5.000   2.000   3.000  2.0000 1.5000\n"
        "ATOM      3  H   ALA A   1       1.000   3.000   3.000  0.2719 0.6000\n"
    )
    s = Structure.from_file(pqr)
    assert s.is_hydrogen.tolist() == [False, False, True]


@pytest.mark.parametrize("bad", [0, -1.0, float("nan")])
def test_invalid_spacing_is_rejected(pdb_file, bad):
    with pytest.raises(ValueError, match="spacing"):
        compute_mlp(Structure.from_file(pdb_file), spacing=bad)


def test_invalid_mep_option_is_rejected(pqr_file):
    with pytest.raises(ValueError, match="eps"):
        compute_mep(Structure.from_file(pqr_file), eps=0)


def test_buckingham_grid_is_finite(tmp_path):
    pdb = tmp_path / "one.pdb"
    pdb.write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\n")
    gm = compute_mlp(Structure.from_file(pdb), method="buckingham", spacing=1.0)
    assert np.isfinite(gm.data).all()


def test_cli_reports_errors_without_traceback(pdb_file, capsys):
    from msgrid.cli import main

    assert main([str(pdb_file), "--no-mlp", "--no-mep"]) == 1
    assert "Nothing to do" in capsys.readouterr().err

    assert main([str(pdb_file)]) == 1  # PDB has no charges, MEP enabled
    assert "charges" in capsys.readouterr().err


def test_cli_runs_both_maps(pqr_file, tmp_path):
    from msgrid.cli import main

    rc = main([str(pqr_file), "-d", str(tmp_path / "out"), "--mlp-spacing", "2", "--mep-spacing", "1.5"])
    assert rc == 0
    assert (tmp_path / "out" / "mol_MLP.dx").exists()
    assert (tmp_path / "out" / "mol_MEP.dx").exists()


def test_standalone_tools_report_errors(pdb_file, tmp_path, capsys):
    from msgrid import mep, mlp

    assert mlp.main(["-f", str(pdb_file), "-o", str(tmp_path / "o.dx"), "-s", "0"]) == 1
    assert "spacing" in capsys.readouterr().err
    assert mep.main([str(pdb_file), str(tmp_path / "o.dx")]) == 1
    assert "charges" in capsys.readouterr().err


@pytest.mark.parametrize("fmt", ["pqr", "auto"])
def test_pdb_content_is_not_read_as_pqr(pdb_file, tmp_path, fmt):
    fake = tmp_path / "actually_pdb.pqr"
    fake.write_text(pdb_file.read_text())
    with pytest.raises(ValueError, match="looks like a PDB"):
        Structure.from_file(fake, fmt=fmt)
