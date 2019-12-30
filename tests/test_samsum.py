
def test_overlapping_intervals():
    from samsum import alignment_utils
    coords_one = (0, 100)
    coords_two = (50, 101)
    coords_three = (101, 200)
    assert alignment_utils.overlapping_intervals(coords_one, coords_two) is True
    assert alignment_utils.overlapping_intervals(coords_two, coords_three) is True
    assert alignment_utils.overlapping_intervals(coords_one, coords_three) is False
