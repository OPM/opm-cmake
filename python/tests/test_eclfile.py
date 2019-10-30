import unittest
import sys
import numpy as np

from opm.io.ecl import EclFile


class TestEclFile(unittest.TestCase):
    
    def test_getListOfArrays(self):

        refList=["INTEHEAD","LOGIHEAD","DOUBHEAD","PORV","DEPTH","DX","DY","DZ","PORO",
                 "PERMX","PERMY", "PERMZ","NTG","TRANX","TRANY","TRANZ","TABDIMS","TAB",
                 "ACTNUM","EQLNUM","FIPNUM","PVTNUM","SATNUM","TRANNNC"]
    
        self.assertRaises(ValueError, EclFile, "/file/that/does_not_exists")

        file2uf = EclFile("tests/data/SPE9.INIT")
        self.assertEqual(file2uf.getNumArrays(), 24)

        arrList = file2uf.getListOfArrays()
        self.assertEqual(arrList, refList)

        file2f = EclFile("tests/data/SPE9.FINIT")
        self.assertEqual(file2f.getNumArrays(), 24)

    def test_get_function(self):

        file1 = EclFile("tests/data/SPE9.INIT")

        first = file1.get(0)
        self.assertEqual(len(first), 95)

        # get fourth array in file SPE9.INIT which is PORV
        test1 = file1.get(3)
        test2 = file1.get("PORV")

        for val1, val2 in zip(test1, test2):
            self.assertEqual(val1, val2)

    def test_get_function_float(self):

        file1 = EclFile("tests/data/SPE9.INIT")
         
        dzList=[20.0, 15.0, 26.0, 15.0, 16.0, 14.0, 8.0, 8.0, 18.0, 12.0, 19.0, 18.0, 20.0, 50.0, 100.0]
        poroList = [0.087, 0.097, 0.111, 0.16, 0.13, 0.17, 0.17, 0.08, 0.14, 0.13, 0.12, 0.105, 0.12, 0.116, 0.157]
        ft3_to_bbl = 0.1781076
        
        refporv = []
        
        for poro, dz in zip(dzList, poroList):
            for i in range(0,600):
                refporv.append(300.0*300.0*dz*poro*ft3_to_bbl)
        
        self.assertTrue(file1.hasArray("PORV"))
        porv_np = file1.get("PORV")
        
        self.assertEqual(len(porv_np), 9000)

        self.assertTrue(isinstance(porv_np, np.ndarray))
        self.assertEqual(porv_np.dtype, "float32")

        porv_list = file1.get("PORV")
        
        for val1, val2 in zip(porv_np, refporv):
            self.assertLess(abs(1.0 - val1/val2), 1e-6)

            
    def test_get_function_double(self):

        refTabData=[0.147E+02, 0.2E+21, 0.4E+03, 0.2E+21, 0.8E+03, 0.2E+21, 0.12E+04, 0.2E+21, 0.16E+04, 0.2E+21, 0.2E+04, 0.2E+21, 0.24E+04, 0.2E+21, 0.28E+04, 0.2E+21, 0.32E+04, 0.2E+21, 0.36E+04, 0.2E+21, 0.4E+04, 0.5E+04, 0.1E+01, 0.2E+21, 0.98814229249012E+00, 0.2E+21, 0.97513408093613E+00]

        file1 = EclFile("tests/data/SPE9.INIT")
        
        tab = file1.get("TAB")

        self.assertTrue(isinstance(tab, np.ndarray))
        self.assertEqual(tab.dtype, "float64")
            
        for i in range(0, len(refTabData)):
            self.assertLess(abs(1.0 - refTabData[i]/tab[i]), 1e-12 )


    def test_get_function_integer(self):

        refTabdims = [ 885, 1, 1, 1, 1, 1, 1, 67, 11, 2, 1, 78, 1, 78, 78, 0, 0, 0, 83, 1, 686, 40, 1, 86, 40, 1, 
                      286, 1, 80, 1 ]
        
        file1 = EclFile("tests/data/SPE9.INIT")
        tabdims = file1.get("TABDIMS")

        self.assertTrue(isinstance(tabdims, np.ndarray))
        self.assertEqual(tabdims.dtype, "int32")

        for i in range(0, len(refTabdims)):
            self.assertEqual(refTabdims[i], tabdims[i])


    def test_get_function_logi(self):

        file1 = EclFile("tests/data/9_EDITNNC.INIT")

        self.assertTrue(file1.hasArray("LOGIHEAD"))
        logih = file1.get("LOGIHEAD")

        self.assertEqual(len(logih), 121)
        self.assertEqual(logih[0], True)
        self.assertEqual(logih[2], False)
        self.assertEqual(logih[8], True)
        
    def test_get_function_char(self):
        
        file1 = EclFile("tests/data/9_EDITNNC.SMSPEC")

        self.assertTrue(file1.hasArray("KEYWORDS"))
        keyw = file1.get("KEYWORDS")
    
        self.assertEqual(len(keyw), 312)
        self.assertEqual(keyw[0], "TIME")
        self.assertEqual(keyw[16], "FWCT")
                 
if __name__ == "__main__":

    unittest.main()
    
