from src.lib import rnastructure_pf as rpf


def test_parse_probability_plot_formats() -> None:
    text = "\n".join(
        [
            "# comment",
            "1 10 0.2",
            "2 9 3.0",
            "3 8 -1.0",
        ]
    )
    pairs = rpf.parse_probability_plot(text)
    assert pairs[(0, 9)] == 0.2
    assert abs(pairs[(1, 8)] - 1e-3) < 1e-9
    assert (2, 7) not in pairs
