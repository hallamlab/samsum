import unittest


class MyTestCase(unittest.TestCase):
    def test_overlapping_intervals(self):
        from samsum import alignment_utils
        coords_one = (0, 100)
        coords_two = (50, 101)
        coords_three = (101, 200)
        self.assertTrue(alignment_utils.overlapping_intervals(coords_one, coords_two))
        self.assertTrue(alignment_utils.overlapping_intervals(coords_two, coords_three))
        self.assertFalse(alignment_utils.overlapping_intervals(coords_one, coords_three))
        return


if __name__ == '__main__':
    unittest.main()
