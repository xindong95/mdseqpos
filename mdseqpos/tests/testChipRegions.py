"""
SYNOPSIS: Tests for the Motif class
"""
import sys
import unittest

sys.path.append('..')
from lib.chipregions import ChipRegions
class AcceptanceTests(unittest.TestCase):
    
    def test_init_and_read(self):
        """Test that ChipRegions are read in on initialization"""
        chrom = ['chr4' for i in range(10)]
        chromStart = [3006483, 3011755, 6371870, 6403510, 6608743, 6851518,
                      7099291, 7125877, 7133863, 7616485]
        chromEnd = [3007652, 3012807, 6372866, 6404169, 6609367, 6852478,
                    7100429, 7127047, 7134524, 7617329]

        foo = ChipRegions('fixtures/test.bed', 'mm8')
        pairs = [(chrom, foo.chrom),
                 (chromStart, foo.chromStart),
                 (chromEnd, foo.chromEnd)]

        for p in pairs:
            for (ref, val) in zip(p[0], p[1]):
                #print ref, val
                self.assertEqual(ref, val)

        self.assertTrue(foo.genome, 'mm8')

    def test_read_sequence(self):
        foo = ChipRegions('fixtures/test.bed', 'mm8')
        foo.read_sequence()
        #print foo.sequence

    def test_bad_region(self):
        """When reading in a region whose chr name is bad, it should be 
        ignored, NOT none"""
        foo = ChipRegions('fixtures/bad.bed', 'mm8')
        foo.read_sequence()
        self.assertEqual(len(foo.sequence), 0)

    def test_ext_read_sequence(self):
        pass

    def test_len_read_sequence(self):
        pass

    def test_to_fasta(self):
        pass
    def test_to_bed(self):
        pass
    
    def test_trim(self):
        """Trim method trims the region down to a certain length centered
        on a mid-point, default = 600bps"""
        chromStart = [3006767, 3011981, 6372068, 6403539, 6608755, 6851698,
                      7099560, 7126162, 7133893, 7616607]
        chromEnd = [3007367, 3012581, 6372668, 6404139, 6609355, 6852298,
                    7100160, 7126762, 7134493, 7617207]
        
        foo = ChipRegions('fixtures/test.bed', 'mm8')
        foo.trim()

        pairs = [(chromStart, foo.chromStart),
                 (chromEnd, foo.chromEnd)]

        for p in pairs:
            for (ref, val) in zip(p[0], p[1]):
                #print ref, val
                self.assertEqual(ref, val)

    def test_mdmodule(self):
        pass
    def test_motifscan(self):
        pass

if __name__ == '__main__':
    unittest.main()
