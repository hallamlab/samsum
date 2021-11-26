import os
import pytest
import unittest


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.output_dir = "tmp_dir"

    def tearDown(self) -> None:
        from shutil import rmtree
        if os.path.isdir(self.output_dir):
            rmtree(self.output_dir)
        return

    def test_int_to_str(self):
        from samsum import utilities
        self.assertEqual('001', utilities.int_to_str(1, 3))
        self.assertEqual('101', utilities.int_to_str(101, 3))
        with pytest.raises(SystemExit):
            utilities.int_to_str(101, 2)
        return

    def test_make_sure_dir_exists(self):
        from samsum import utilities
        self.assertFalse(os.path.isdir(self.output_dir))
        utilities.make_sure_dir_exists(dir_path=self.output_dir)
        self.assertTrue(os.path.isdir(self.output_dir))
        return


if __name__ == '__main__':
    unittest.main()
