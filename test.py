import unittest
import valley_centerline
import geopandas as gpd

class TestCenterlineExtraction(unittest.TestCase):
    def test_centerline_geometry(self):
        whitewater = gpd.read_file("Input\SampleData\Whitewater\WhitewaterWalls.gpkg").explode()
        self.assertTrue(valley_centerline.extract_centerline(whitewater).is_valid[0])

if __name__ == 'main':
    unittest.main()