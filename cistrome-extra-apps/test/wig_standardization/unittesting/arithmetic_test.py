from arithmetic import MyMath
import unittest

class TestSequenceFunctions(unittest.TestCase):
    def setUp(self):
        pass
        
    def tearDown(self):
        pass

    def testsum(self):
        self.assertEquals(MyMath.sum(2,3), 5)

    def testsum2(self):
        self.assertEquals(MyMath.sum(2,-2),0)

    def testsumBAD(self):
        self.assertRaises(Exception, MyMath.sum, 2, "a")
        
    def testsumBAD2(self):
        self.assertRaises(Exception, MyMath.sum, 2, [])

    def testNEWBUG(self):
        self.assertRaises(Exception, MyMath.sum, 2, [])
        
if __name__ == '__main__':
        unittest.main()
        
