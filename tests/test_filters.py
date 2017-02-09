import unittest
from sidekick import filters

class Read(object):
    """
    Dummy version of a pysam Bam read
    """
    def __init__(self, is_unmapped, mate_is_unmapped,
                 is_paired, is_supplementary,
                 mapping_quality):
        self.is_unmapped = is_unmapped
        self.mate_is_unmapped = mate_is_unmapped
        self.is_paired = is_paired
        self.is_supplementary = is_supplementary
        self.mapping_quality = mapping_quality


class HeroTests(unittest.TestCase):

    def setUp(self):
        """
        Construct some reads that should pass or fail
        filterfn1
        """
        self.passing = Read(False, True, True, False, 60)
        self.lowqual = Read(False, True, True, False, 45)
        self.failing = Read(True, False, True, False, 60)
    
    def test_passing_passes(self):
        self.assertTrue(filters.hero(self.passing))

    def test_lowqual_fails(self):
        self.assertFalse(filters.hero(self.lowqual))

    def test_failing_fails(self):
        self.assertFalse(filters.hero(self.failing))


class SidekickTests(unittest.TestCase):

    def setUp(self):
        """
        Construct some reads that should pass or fail
        filterfn1
        """
        self.passing = Read(True, False, True, False, 60)
        self.failing = Read(False, True, True, False, 60)
    
    def test_passing_passes(self):
        self.assertTrue(filters.sidekick(self.passing))

    def test_failing_fails(self):
        self.assertFalse(filters.sidekick(self.failing))


if __name__ == '__main__':
    unittest.main()
