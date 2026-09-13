import unittest

from scripts.workflow_progress import progress_value, unread_lines


class WorkflowProgressTests(unittest.TestCase):
    def test_formats_fraction_as_whole_percent(self):
        value, label = progress_value(["Finished job"] * 11, total=15)
        self.assertAlmostEqual(value, 11 / 15)
        self.assertEqual(label, "73%")

    def test_each_page_can_read_the_complete_existing_log(self):
        lines = ["first", "second"]
        first_page, first_cursor = unread_lines(lines, 0)
        second_page, second_cursor = unread_lines(lines, 0)
        self.assertEqual(first_page, lines)
        self.assertEqual(second_page, lines)
        self.assertEqual(first_cursor, 2)
        self.assertEqual(second_cursor, 2)


if __name__ == "__main__":
    unittest.main()
