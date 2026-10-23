from chrom_core.trim import trim_window


def test_trims_low_quality_flanks():
    # high-quality span is bases 30..129; boundary-straddling windows still
    # pass the mean threshold (e.g. window at 17: (13*10 + 7*40)/20 = 20.5),
    # so the span extends a few bases into the low-quality flanks
    quality = [10] * 30 + [40] * 100 + [10] * 30
    assert trim_window(quality, window=20, min_q=20) == (17, 143)


def test_all_good_is_full_length():
    assert trim_window([30] * 50, window=20, min_q=20) == (0, 50)


def test_all_bad_returns_empty():
    assert trim_window([5] * 50, window=20, min_q=20) == (0, 0)


def test_empty_input():
    assert trim_window([], window=20, min_q=20) == (0, 0)


def test_window_larger_than_read_falls_back_to_whole_read():
    assert trim_window([40, 40, 40], window=20, min_q=20) == (0, 3)


def test_mean_not_median_one_bad_base_keeps_window():
    # window mean (19*40 + 0)/20 = 38 >= 20 -> still good
    quality = [40] * 19 + [0] + [40] * 30
    assert trim_window(quality, window=20, min_q=20) == (0, 50)
