from refmasker.detectors.repeats import RepeatHit
from refmasker.regions import Region
from refmasker.seqio import Record
from refmasker.summary import repeat_sensitivity, summarise


def test_total_counts_overlapping_categories_once():
    recs = [Record("a", "", "A" * 1000), Record("b", "", "A" * 1000)]
    regions = [
        Region("a", 0, 100, "repeats"),
        Region("a", 50, 150, "is"),
        Region("b", 0, 10, "is"),
    ]
    rows = {(r.record, r.category): r for r in summarise(regions, recs, ["repeats", "is"])}
    assert rows[("ALL", "repeats")].bp == 100
    assert rows[("ALL", "is")].bp == 110
    assert rows[("ALL", "total")].bp == 160
    assert rows[("ALL", "total")].pct == 8.0
    assert rows[("b", "total")].bp == 10
    assert rows[("a", "total")].regions == 1


def test_sensitivity_marks_chosen_setting():
    hits = [
        RepeatHit("a", 0, 150, "a", 500, 650, 150, 99.0),
        RepeatHit("a", 1000, 1400, "a", 2000, 2400, 400, 85.0),
    ]
    rows = repeat_sensitivity(hits, 10000, 200, 90.0, (100, 200), (80.0, 90.0))
    by = {(r.min_len, r.min_id): r for r in rows}
    assert by[(100, 80.0)].bp == 300 + 800
    assert by[(200, 90.0)].bp == 0
    assert by[(200, 90.0)].chosen
    assert sum(r.chosen for r in rows) == 1
