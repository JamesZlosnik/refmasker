"""Output folder safety: ownership, permissions, symbolic links and clean-up of work/."""

import os

import pytest

from refmasker.pipeline import clean_workdir, prepare_outdir


def test_new_outdir_gives_others_no_access(tmp_path):
    old = os.umask(0o022)
    try:
        prepare_outdir(tmp_path / "out", force=False)
    finally:
        os.umask(old)
    assert (tmp_path / "out").stat().st_mode & 0o777 == 0o750


def test_refuses_folder_of_another_user(tmp_path, monkeypatch):
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr(os, "getuid", lambda: out.stat().st_uid + 1)
    with pytest.raises(PermissionError, match="another user"):
        prepare_outdir(out, force=False)


def test_refuses_world_writable_folder(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    out.chmod(0o777)
    with pytest.raises(PermissionError, match="every user"):
        prepare_outdir(out, force=False)


def test_refuses_a_file(tmp_path):
    out = tmp_path / "out"
    out.write_text("x")
    with pytest.raises(FileExistsError, match="not a folder"):
        prepare_outdir(out, force=True)


def test_force_refuses_symlinks_in_previous_output(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "refmasker.log").write_text("")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (out / "work").symlink_to(elsewhere)
    with pytest.raises(FileExistsError, match="symbolic links"):
        prepare_outdir(out, force=True)


def test_clean_workdir_keeps_tool_logs_after_failure(tmp_path):
    work = tmp_path / "work"
    (work / "isescan").mkdir(parents=True)
    (work / "isescan" / "isescan.log").write_text("error details\n")
    (work / "reference.fna").write_text(">chr\nACGT\n")
    clean_workdir(work, tmp_path / "logs", keep_logs=True)
    assert not work.exists()
    assert (tmp_path / "logs" / "isescan.log").read_text() == "error details\n"


def test_clean_workdir_after_success_keeps_nothing(tmp_path):
    work = tmp_path / "work"
    (work / "isescan").mkdir(parents=True)
    (work / "isescan" / "isescan.log").write_text("")
    clean_workdir(work, tmp_path / "logs", keep_logs=False)
    assert not work.exists()
    assert not (tmp_path / "logs").exists()
