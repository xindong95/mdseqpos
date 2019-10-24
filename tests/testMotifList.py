"""
SYNOPSIS: Tests for the MotifList class
"""
import sys
import unittest

sys.path.append('..')
from lib.motif import MotifList
from lib.motif import Motif
from lib.motif import DuplicateMotif

class AcceptanceTests(unittest.TestCase):
    def test_no_duplicates(self):
        tmp = Motif()
        ml = MotifList()
        ml.append(tmp)
        #adding it again MUST raise DuplicateMotif Exception
        self.assertRaises(DuplicateMotif, ml.append, tmp)

    def test_append(self):
        ml = MotifList()
        self.assertEqual(len(ml), 0)
        ml.append(Motif())
        self.assertEqual(len(ml), 1)

    def test_cull(self):
        """Test that the cull w/ default params works"""
        pass

    def test_cull_maxmotif(self):
        """Test the max motifs are honored by cull"""
        pass

    def test_cull_pval(self):
        """Test that pval is honored by cull"""
        pass

    def test_from_xml_file(self):
        """1. write a sample xml file, 2. try to read it in"""
        pass

    def test_to_xml(self):
        """NOTE: I'm NOT GOING TO TEST THIS FN b/c the output is fragile"""
        pass

    def test_save_to_xml(self):
        """
        1. Read in the sample xml file (from_xml) --> L1
        2. save list xml --> foo.xml
        3. Read in foo.xml --> L2
        4. assertEqual(L1, L2)
        """
        pass

if __name__ == '__main__':
    unittest.main()
