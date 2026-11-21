from refmasker.external import writable_by_others


def test_writable_by_others(tmp_path):
    tool = tmp_path / "bin" / "nucmer"
    tool.parent.mkdir(mode=0o755)
    tool.write_text("")
    tool.chmod(0o755)
    tool.parent.chmod(0o755)
    assert not writable_by_others(tool)
    tool.chmod(0o757)
    assert writable_by_others(tool)
    tool.chmod(0o755)
    tool.parent.chmod(0o777)  # anyone could replace the file
    assert writable_by_others(tool)
    tool.parent.chmod(0o1777)  # sticky, like /tmp: others cannot replace your files
    assert not writable_by_others(tool)
    tool.parent.chmod(0o755)
