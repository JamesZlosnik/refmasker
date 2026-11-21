import os

import pytest

from refmasker.safeio import open_replacing


def test_replaces_a_planted_symlink_instead_of_following_it(tmp_path):
    victim = tmp_path / "victim.txt"
    victim.write_text("important\n")
    link = tmp_path / "mask.bed"
    link.symlink_to(victim)
    with open_replacing(link) as fh:
        fh.write("chr\t0\t10\n")
    assert victim.read_text() == "important\n"
    assert not link.is_symlink()
    assert link.read_text() == "chr\t0\t10\n"


def test_failed_write_leaves_old_file_and_no_temp(tmp_path):
    out = tmp_path / "summary.tsv"
    out.write_text("old\n")
    with pytest.raises(RuntimeError), open_replacing(out) as fh:
        fh.write("half")
        raise RuntimeError("boom")
    assert out.read_text() == "old\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["summary.tsv"]


def test_permissions_follow_umask(tmp_path):
    old = os.umask(0o027)
    try:
        with open_replacing(tmp_path / "a.txt") as fh:
            fh.write("x")
    finally:
        os.umask(old)
    assert (tmp_path / "a.txt").stat().st_mode & 0o777 == 0o640


def test_binary_mode(tmp_path):
    with open_replacing(tmp_path / "a.bin", "wb") as fh:
        fh.write(b"\x00\x01")
    assert (tmp_path / "a.bin").read_bytes() == b"\x00\x01"
