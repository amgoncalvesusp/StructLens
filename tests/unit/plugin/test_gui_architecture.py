"""Hard architectural boundaries for the report-first Qt plugin."""

from pathlib import Path


def test_every_gui_production_module_stays_within_the_800_line_hard_limit() -> None:
    gui_root = Path(__file__).parents[2].parent / "src" / "structlens" / "plugin" / "gui"
    line_counts = {
        path.relative_to(gui_root).as_posix(): len(path.read_text(encoding="utf-8").splitlines())
        for path in gui_root.rglob("*.py")
    }
    oversized = {name: count for name, count in line_counts.items() if count > 800}

    assert not oversized, f"GUI modules exceed the hard 800-line limit: {oversized}"
