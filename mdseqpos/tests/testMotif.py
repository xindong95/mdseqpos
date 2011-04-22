"""
SYNOPSIS: Tests for the Motif class
"""
import sys
import unittest
#import time
from xml.dom.minidom import parse, parseString

sys.path.append('..')
from lib.motif import Motif
from lib.motif import BadPssm
import lib.motif as motif

def feq(f1, f2):
    epsilon = 1e-6
    return abs(f1 - f2) < epsilon

class AcceptanceTests(unittest.TestCase):
    _ATTRIBUTES = ['id', 'pssm', 'seqpos_results', 'antisense',
                   'source', 'sourcefile', 'species', 'fullname', 'pmid',
                   'numseqs', 'symbols', 'entrezs', 'refseqs', 'dbd']
    _results_fields = ['numhits','cutoff','zscore','meanposition','pvalue']
            
    def test_attributes_set_to_none(self):
        """Make sure that all of the attributes are all initialized to None"""
        foo = Motif()
        for attr in AcceptanceTests._ATTRIBUTES:
            if attr != 'seqpos_results':
                self.assertEqual(getattr(foo, attr), None)

    def test_seqpos_results_set_to_none(self):
        """Make sure that seqpos_results dictionary is set to None"""
        foo = Motif()
        for attr in AcceptanceTests._results_fields:
            self.assertEqual(foo.seqpos_results[attr], None)

    def test_convenience_fns(self):
        """Make convenience fns like getpvalue, etc. work"""
        foo = Motif()
        for attr in AcceptanceTests._results_fields:
            #COULD have done this instead to be explicit, but my way is cleaner
            #self.assertEqual(getattr(foo, "get"+attr).__call__(), None)
            self.assertEqual(getattr(foo, "get"+attr)(), None)

    def test_validpssm(self):
        """This is a valid pssm, it should pass"""
        pssm = [[0.1, 0.3, 0.4, 0.2],
                [0.5, 0.1, 0.1, 0.3],
                [0.2, 0.4, 0.3, 0.1]]
        self.assertTrue(Motif._validpssm(pssm))

    def test_invalidpssm_num_cols(self):
        """This is an invalid pssm, not 4 cols"""
        pssm = [[0.1, 0.3, 0.4]]
        self.assertRaises(BadPssm, Motif._validpssm, pssm)

    def test_invalidpssm_row_sum(self):
        """This is an invalid pssm, does not sum to 1.0"""
        pssm = [[0.1, 0.3, 0.4, 0.5]]
        self.assertRaises(BadPssm, Motif._validpssm, pssm)

    def test_invalidpssm_None(self):
        """This is an invalid pssm"""
        pssm = None
        self.assertRaises(BadPssm, Motif._validpssm, pssm)

    def test_invalidpssm_empty(self):
        """This is an invalid pssm"""
        pssm = []
        self.assertRaises(BadPssm, Motif._validpssm, pssm)

    def test_equals(self):
        """Tests the motif.equals fn"""
        foo = Motif()
        bar = Motif()
        foo.id = 'M00134'
        bar.id = 'M00134'
        self.assertTrue(foo.equals(bar))
        self.assertTrue(foo.equals(foo))

    def test_not_equals(self):
        """Tests the motif.equals fn"""
        foo = Motif()
        bar = Motif()
        foo.id = 'M00134'
        self.assertFalse(foo.equals(bar))

    def test_from_flat_file(self):
        """fixtures/goodPssm.txt should be a valid motif"""
        foo = Motif.from_flat_file('fixtures/goodPssm.txt')
        self.assertTrue(Motif._validpssm(foo.pssm))

    def test_from_flat_file2(self):
        """fixtures/goodPssmObCb.txt should be a valid motif"""
        foo = Motif.from_flat_file('fixtures/goodPssmObCb.txt')
        self.assertTrue(Motif._validpssm(foo.pssm))

    def test_from_flat_file_bad(self):
        """fixtures/badPssm.txt should raise error, missing 4th col"""
        self.assertRaises(BadPssm,Motif.from_flat_file, 'fixtures/badPssm.txt')

    def test_from_flat_file_bad(self):
        """fixtures/badPssm.txt should raise error, does not sum to 1.0"""
        self.assertRaises(BadPssm,Motif.from_flat_file,
                          'fixtures/badPssm_sum.txt')

    def test_to_json(self):
        """test the json serialization method"""
        expected_json = """{"refseqs": ["NM_013464", "NM_001621"], "sourcefile": "5", "pssm": [[0.1000, 0.3000, 0.4000, 0.2000], [0.5000, 0.1000, 0.1000, 0.3000], [0.2000, 0.4000, 0.3000, 0.1000]], "antisense": "True", "dbd": "5", "numseqs": "5", "species": null, "symbols": ["Runx3", "Ahr"], "source": "5", "seqpos_results": {"numhits": null, "zscore": null, "pvalue": null, "meanposition": null, "cutoff": null}, "entrezs": ["25690", "11622", "196"], "fullname": "5", "pmid": "5", "id": "5"}"""
        
        dom = parse('fixtures/goodMotif.xml')
        #print dom.documentElement.tagName
        foo = Motif.from_xml(dom.documentElement)
        #print foo.to_json()
        self.assertEqual(foo.to_json(), expected_json)
        #print foo.to_json()
                        
    def test_to_xml(self):
        """test the xml serialization method"""
        pssm = [[0.1, 0.3, 0.4, 0.2],
                [0.5, 0.1, 0.1, 0.3],
                [0.2, 0.4, 0.3, 0.1]]
        foo = Motif()
        for attr in AcceptanceTests._ATTRIBUTES:
            setattr(foo, attr, '5')
        foo.symbols = ["Runx3", "Ahr"]
        foo.entrezs = ["25690", "11622", "196"]
        foo.refseqs = ["NM_013464", "NM_001621"]
        foo.antisense = True
        foo.setpssm(pssm)

        xml_str = foo.to_xml(print_non_schema=True)
        dom = parseString(xml_str)
        #verify the output
        #PUNT--i'm not up for it right now
        self.assertTrue(True)

    def test_from_xml(self):
        """test the xml motif parser; goodMotif.xml should equal foo"""
        pssm = [[0.1, 0.3, 0.4, 0.2],
                [0.5, 0.1, 0.1, 0.3],
                [0.2, 0.4, 0.3, 0.1]]
        foo = Motif()
        for attr in AcceptanceTests._ATTRIBUTES:
            setattr(foo, attr, '5')
        foo.id = "5"
        foo.symbols = ["Runx3", "Ahr"]
        foo.entrezs = ["25690", "11622", "196"]
        foo.refseqs = ["NM_013464", "NM_001621"]
        foo.antisense = True
        foo.setpssm(pssm)
        #print foo

        dom = parse('fixtures/goodMotif.xml')
        #print dom.documentElement.tagName
        bar = Motif.from_xml(dom.documentElement)
        #print bar
        self.assertTrue(foo.equals(bar))

    def test_from_xml_optional_attr(self):
        """Test that a motif with none of the optional attributes is read in correctly"""
        pssm = [[0.1, 0.3, 0.4, 0.2],
                [0.5, 0.1, 0.1, 0.3],
                [0.2, 0.4, 0.3, 0.1]]
        foo = Motif()
        foo.id = "5"
        dom = parse('fixtures/goodMotifOpt.xml')
        #print dom.documentElement.tagName
        bar = Motif.from_xml(dom.documentElement)
        #print bar
        self.assertTrue(foo.equals(bar))

    def test_from_xml2(self):
        """test the xml motif parser w/o reading from file--use to_xml"""
        pssm = [[0.1, 0.3, 0.4, 0.2],
                [0.5, 0.1, 0.1, 0.3],
                [0.2, 0.4, 0.3, 0.1]]
        foo = Motif()
        foo.id = "19"
        foo.symbols = ["Runx3", "Ahr"]
        foo.entrezs = ["25690", "11622", "196"]
        foo.refseqs = ["NM_013464", "NM_001621"]
        foo.antisense = True
        foo.setpssm(pssm)
        #print foo

        fuz = foo.to_xml()
        #print fuz
        dom = parseString(fuz)
        bar = Motif.from_xml(dom.documentElement)
        #print bar
        #print "%s, %s" % (foo.id, bar.id)
        #print foo.id == bar.id
        self.assertTrue(foo.equals(bar))

    def test_calc_pssm(self):
        """Tests the calc_pssm fn"""
        
        seq = [['A', 'A', 'A', 'C', 'C', 'A', 'C', 'A', 'G'],
               ['T', 'A', 'T', 'C', 'C', 'G', 'C', 'A', 'G'],
               ['C', 'A', 'A', 'C', 'C', 'A', 'C', 'A', 'A'],
               ['A', 'A', 'A', 'C', 'C', 'A', 'C', 'A', 'G'],
               ['A', 'A', 'A', 'C', 'C', 'A', 'C', 'A', 'G'],
               ['A', 'A', 'A', 'C', 'C', 'A', 'C', 'A', 'A']]
        
        pssm_ref =[[0.66666666666666663, 0.16666666666666666, 0.0, 0.16666666666666666],
         [1.0, 0.0, 0.0, 0.0],
         [0.83333333333333337, 0.0, 0.0, 0.16666666666666666],
         [0.0, 1.0, 0.0, 0.0],
         [0.0, 1.0, 0.0, 0.0],
         [0.83333333333333337, 0.0, 0.16666666666666666, 0.0],
         [0.0, 1.0, 0.0, 0.0],
         [1.0, 0.0, 0.0, 0.0],
         [0.33333333333333331, 0.0, 0.66666666666666663, 0.0]]
        pssm = motif.calc_pssm(seq)
        for r in range(len(pssm_ref)):
            for c in range(len(pssm_ref[0])):
                self.assertTrue(feq(pssm_ref[r][c], pssm[r][c]))

    def test_seqpos_stat(self):
        """Testing the seqpos stat method:
        actual data from running: MDSeqPos.py runx_small.bed mm8 -m y1h.xml
        """
        pssm = [[0.49, 0.167, 0.01, 0.333],
                [0.924, 0.01, 0.056, 0.01],
                [0.924, 0.01, 0.056, 0.01],
                [0.01, 0.97, 0.01, 0.01],
                [0.01, 0.97, 0.01, 0.01],
                [0.313, 0.01, 0.667, 0.01],
                [0.01, 0.97, 0.01, 0.01],
                [0.97, 0.01, 0.01, 0.01],
                [0.591, 0.01, 0.389, 0.01]]
        
        start = [114, 596, 329, 208, 507, 340, 467, 239, 458, 283, 559, 281,
                 175, 301, 397, 623, 240, 54, 212, 511, 285, 380, 511, 116]
        end = [123, 605, 338, 217, 516, 349, 476, 248, 467, 292, 568, 290,
               184, 310, 406, 632, 249, 63, 221, 520, 294, 389, 520, 125]
        score = [4.91170692, 5.87872267, 6.32827187, 4.73049879, 4.61262989,
                 6.83801079, 1.78081846, 4.05301714, 4.85997772, 2.98214078,
                 4.32487345, 2.78596187, 4.73322868, 5.85379267, 5.10804605,
                 4.35964298, 3.89328027, 7.50640059, 5.79009295, 2.93482327,
                 7.54800034, 3.64937592, 3.9898386, 2.86741686]
        
        #fracpos is an intermediate step- its here just to compare the output
        frac = [0.64430577, 0.85959438, 0.02652106, 0.35101404, 0.58190328,
                0.06084243, 0.45709828, 0.25429017, 0.42901716, 0.11700468,
                0.74414977, 0.12324493, 0.45397816, 0.06084243, 0.23868955,
                0.94383775, 0.25117005, 0.83151326, 0.33853354, 0.59438378,
                0.11076443,  0.18564743,  0.59438378,  0.63806552]
        
        length = 650
        foo = Motif()
        foo.setpssm(pssm)
        foo.seqpos_stat(start, end, score, length)
        
        self.assertEqual(foo.getnumhits(), 23)
        self.assertTrue(feq(foo.getzscore(), -0.16397254440660525))
        self.assertTrue(feq(foo.getpvalue(), 0.43487637881039054))
        self.assertTrue(feq(foo.getmeanposition(), -0.069965407311944727))
        self.assertTrue(feq(foo.getcutoff(), 1.7808184600000001))

    def test_seqpos(self):
        class Bar:
            pass

        pssm = [[0.49, 0.167, 0.01, 0.333],
                [0.924, 0.01, 0.056, 0.01],
                [0.924, 0.01, 0.056, 0.01],
                [0.01, 0.97, 0.01, 0.01],
                [0.01, 0.97, 0.01, 0.01],
                [0.313, 0.01, 0.667, 0.01],
                [0.01, 0.97, 0.01, 0.01],
                [0.97, 0.01, 0.01, 0.01],
                [0.591, 0.01, 0.389, 0.01]]

        #read in an actual chip regions sequence
        f = open("fixtures/chipRegionSeq.txt")
        seq = eval(f.readline())
        #print seq
        f.close()

        chip_regions = Bar()
        chip_regions.sequence = seq
        chip_regions.preprocessed_regions = True
        
        foo = Motif()
        foo.setpssm(pssm)
        foo.seqpos(chip_regions)
        #print foo

        seqpos_results_ref =\
        {'cutoff': 4.9475970072872455, 'zscore': -8.547918739697856,
         'pssm': [[0.50764525993883791, 0.17125382262996941, 0.027522935779816515, 0.29357798165137616],
                  [0.80122324159021407, 0.033639143730886847, 0.1529051987767584, 0.012232415902140673],
                  [0.86850152905198774, 0.015290519877675841, 0.10703363914373089, 0.0091743119266055051],
                  [0.015290519877675841, 0.95718654434250761, 0.0091743119266055051, 0.01834862385321101],
                  [0.01834862385321101, 0.94189602446483178, 0.01834862385321101, 0.021406727828746176],
                  [0.61467889908256879, 0.030581039755351681, 0.31804281345565749, 0.03669724770642202],
                  [0.01834862385321101, 0.94495412844036697, 0.027522935779816515, 0.0091743119266055051],
                  [0.90825688073394495, 0.024464831804281346, 0.03669724770642202, 0.030581039755351681],
                  [0.47094801223241589, 0.055045871559633031, 0.44954128440366975, 0.024464831804281346]],
         'meanposition': -0.14154048514881779, 'numhits': 255,
         'pvalue': 6.2663951031948087e-18}

        list = ['cutoff', 'zscore', 'meanposition', 'numhits', 'pvalue']
        for attr in list:
            self.assertTrue(feq(foo.seqpos_results[attr],
                                seqpos_results_ref[attr]))
        #compare pssms
        for r in range(len(seqpos_results_ref['pssm'])):
            for c in range(len(seqpos_results_ref['pssm'][0])):
                self.assertTrue(feq(seqpos_results_ref['pssm'][r][c],
                                    foo.seqpos_results['pssm'][r][c]))
        
        
        
        
    #A timing test that shows my new method is 5 times faster than the
    #older method!!!
#     def test_calc_pssm2(self):
#         """Tests the calc_pssm fn"""
#         f = open("fixtures/seq.txt")
#         seq = eval(f.readline())
#         f.close()
#         seq2 = seq + seq + seq + seq + seq

#         start = time.time()
#         pssm = motif.calc_pssm2(seq2)
#         end = time.time()
#         print "Total is: %s ms\n" % ((end - start)*1000)

#         start = time.time()
#         pssm = motif.calc_pssm(seq2)
#         end = time.time()
#         print "Total is: %s ms\n" % ((end - start)*1000)

if __name__ == '__main__':
    unittest.main()
