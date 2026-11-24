"""The geNomad database checksum manifest: hashing, fingerprints and verification.

No network and no geNomad: these tests build a small fake database folder on disk.
"""

import pytest

from refmasker import database
from refmasker.cli import main


@pytest.fixture
def fake_db(tmp_path):
    """A DEST folder holding a `genomad_db` with a few files, one in a subfolder."""
    db = tmp_path / "genomad_db"
    (db / "sub").mkdir(parents=True)
    (db / "version.txt").write_text("1.9\n")
    (db / "genomad_db.h3i").write_bytes(b"\x00binary profile data" * 10)
    (db / "sub" / "names.tsv").write_text("a\tb\n")
    return db


def test_hash_tree_covers_every_file(fake_db):
    sums = database.hash_tree(fake_db, show_progress=False)
    assert [f.path for f in sums] == ["genomad_db.h3i", "sub/names.tsv", "version.txt"]
    assert all(len(f.sha256) == 64 for f in sums)
    assert sum(f.size for f in sums) == sum(p.stat().st_size for p in database.db_files(fake_db))


def test_hash_tree_rejects_an_empty_folder(tmp_path):
    (tmp_path / "genomad_db").mkdir()
    with pytest.raises(FileNotFoundError):
        database.hash_tree(tmp_path / "genomad_db", show_progress=False)


def test_fingerprint_is_stable_and_order_independent(fake_db):
    sums = database.hash_tree(fake_db, show_progress=False)
    assert database.fingerprint(sums) == database.fingerprint(list(reversed(sums)))
    assert len(database.fingerprint(sums)) == 64


def test_fingerprint_changes_with_contents_and_with_names(fake_db):
    before = database.fingerprint(database.hash_tree(fake_db, show_progress=False))
    (fake_db / "version.txt").write_text("1.10\n")
    after = database.fingerprint(database.hash_tree(fake_db, show_progress=False))
    assert after != before
    (fake_db / "version.txt").rename(fake_db / "version2.txt")
    assert database.fingerprint(database.hash_tree(fake_db, show_progress=False)) != after


def test_manifest_round_trip(fake_db):
    sums = database.hash_tree(fake_db, show_progress=False)
    path = database.write_manifest(fake_db, sums, "1.9")
    assert path == fake_db.parent / "genomad_db.sha256sums"  # beside the database, not inside it
    recorded = database.read_manifest(path)
    assert recorded.sums == {f.path: f.sha256 for f in sums}
    assert recorded.fingerprint == database.fingerprint(sums)
    assert recorded.db_version == "1.9"
    assert recorded.written_by  # the refmasker version that wrote it


def test_manifest_checksum_lines_are_sha256sum_format(fake_db):
    sums = database.hash_tree(fake_db, show_progress=False)
    path = database.write_manifest(fake_db, sums, "1.9")
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    assert len(lines) == len(sums)
    for line in lines:
        digest, sep, name = line.partition("  ")
        assert sep and len(digest) == 64 and int(digest, 16) >= 0
        assert name and not name.startswith(" ")


def test_read_manifest_rejects_junk(tmp_path):
    bad = tmp_path / "x.sha256sums"
    bad.write_text("# only comments\n")
    with pytest.raises(ValueError, match="no checksum lines"):
        database.read_manifest(bad)
    bad.write_text("not-a-checksum-line\n")
    with pytest.raises(ValueError, match="sha256"):
        database.read_manifest(bad)


def test_verify_passes_on_an_untouched_database(fake_db):
    database.write_manifest(fake_db, database.hash_tree(fake_db, show_progress=False), "1.9")
    result, recorded = database.verify(fake_db, show_progress=False)
    assert result.ok
    assert result.checked == 3
    assert result.fingerprint == result.recorded_fingerprint == recorded.fingerprint


@pytest.mark.parametrize(
    "change,field",
    [
        ("edit", "changed"),
        ("delete", "missing"),
        ("add", "added"),
    ],
)
def test_verify_detects_changes(fake_db, change, field):
    database.write_manifest(fake_db, database.hash_tree(fake_db, show_progress=False), "1.9")
    if change == "edit":
        (fake_db / "version.txt").write_text("tampered\n")
    elif change == "delete":
        (fake_db / "version.txt").unlink()
    else:
        (fake_db / "extra.bin").write_bytes(b"planted")
    result, _ = database.verify(fake_db, show_progress=False)
    assert not result.ok
    assert getattr(result, field)
    if change != "add":
        assert result.fingerprint != result.recorded_fingerprint


def test_verify_without_a_manifest_says_so(fake_db):
    with pytest.raises(FileNotFoundError, match="checksum manifest"):
        database.verify(fake_db, show_progress=False)


def test_resolve_db_dir_accepts_the_parent_or_the_database_itself(fake_db, tmp_path):
    assert database.resolve_db_dir(tmp_path) == fake_db
    assert database.resolve_db_dir(fake_db) == fake_db
    with pytest.raises(FileNotFoundError):
        database.resolve_db_dir(tmp_path / "nowhere")


def test_manifest_path_sits_next_to_the_database(tmp_path):
    assert database.manifest_path(tmp_path / "genomad_db").name == "genomad_db.sha256sums"


def test_human_bytes():
    assert database.human_bytes(512) == "512 B"
    assert database.human_bytes(2 * 1024) == "2.0 kB"
    assert database.human_bytes(3 * 1024**3) == "3.0 GB"


def test_cli_verify_round_trip(fake_db, tmp_path, capsys):
    database.write_manifest(fake_db, database.hash_tree(fake_db, show_progress=False), "1.9")
    assert main(["download-db", str(tmp_path), "--verify"]) == 0
    assert "every checksum matches" in capsys.readouterr().out

    (fake_db / "version.txt").write_text("tampered\n")
    assert main(["download-db", str(tmp_path), "--verify"]) == 2
    out = capsys.readouterr().out
    assert "FAILED" in out and "version.txt" in out


def test_cli_download_refuses_to_overwrite_an_existing_database(fake_db, tmp_path, caplog):
    # Exit code 2 and a message, not a traceback and not a deleted database.
    assert main(["download-db", str(tmp_path)]) == 2
    assert "already exists" in caplog.text
    assert (fake_db / "version.txt").is_file()
