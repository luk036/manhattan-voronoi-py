import unittest
import random
import sys
sys.path.insert(0, '.')
from voronoi import (
    _segments_overlap, _point_on_line, _point_on_segment,
    _line_intersection, bisectorIntersection, findL1Bisector,
    samePoint, generateVoronoiPoints, cleanData, generateL1Voronoi,
    distance,
)


class TestSegmentsOverlap(unittest.TestCase):

    def test_collinear_overlap(self):
        self.assertTrue(_segments_overlap([-8, 3], [-4, -1], [-3, -2], [-6, 1]))

    def test_collinear_no_overlap(self):
        self.assertFalse(_segments_overlap([0, 0], [4, 4], [5, 5], [8, 8]))

    def test_not_parallel(self):
        self.assertFalse(_segments_overlap([0, 0], [4, 0], [0, 4], [4, 4]))

    def test_parallel_not_collinear(self):
        self.assertFalse(_segments_overlap([0, 0], [0, 4], [2, 0], [2, 4]))


class TestPointOnLine(unittest.TestCase):

    def test_diagonal(self):
        self.assertTrue(_point_on_line([2, 2], [0, 0], [4, 4]))

    def test_vertical(self):
        self.assertTrue(_point_on_line([2, 0], [2, -5], [2, 10]))

    def test_horizontal(self):
        self.assertTrue(_point_on_line([5, 3], [0, 3], [10, 3]))

    def test_off_diagonal(self):
        self.assertFalse(_point_on_line([3, 4], [0, 0], [4, 4]))


class TestPointOnSegment(unittest.TestCase):

    def test_on_diagonal(self):
        self.assertTrue(_point_on_segment([2, 2], [0, 0], [4, 4]))

    def test_beyond_diagonal(self):
        self.assertFalse(_point_on_segment([5, 5], [0, 0], [4, 4]))

    def test_on_vertical(self):
        self.assertTrue(_point_on_segment([2, 5], [2, 0], [2, 10]))

    def test_on_horizontal(self):
        self.assertTrue(_point_on_segment([2, 3], [0, 3], [10, 3]))


class TestLineIntersection(unittest.TestCase):

    def test_horizontal_cross_vertical(self):
        pt = _line_intersection([0, 0], [4, 0], [2, 2], [2, -2])
        self.assertTrue(samePoint(pt, [2, 0]))

    def test_two_diagonals(self):
        pt = _line_intersection([0, 0], [4, 4], [4, 0], [0, 4])
        self.assertTrue(samePoint(pt, [2, 2]))

    def test_vertical_cross_horizontal(self):
        pt = _line_intersection([2, 0], [2, 4], [0, 1], [4, 1])
        self.assertTrue(samePoint(pt, [2, 1]))

    def test_parallel_diagonals(self):
        self.assertIsNone(_line_intersection([0, 0], [2, 2], [4, 4], [6, 6]))

    def test_parallel_verticals(self):
        self.assertIsNone(_line_intersection([0, 0], [0, 4], [2, 0], [2, 4]))


class TestBisectorIntersection(unittest.TestCase):

    def _s(self, x, y):
        return {'site': [x, y], 'bisectors': []}

    def test_overlap_resolution_identity(self):
        A = self._s(-4, 4); B = self._s(-1, 1)
        C1 = self._s(-8, -2); C2 = self._s(-8, -2)
        bAC = findL1Bisector(A, C1, 200, 200)
        bBC = findL1Bisector(B, C2, 200, 200)
        P = bisectorIntersection(bAC, bBC)
        self.assertTrue(samePoint(P, [-4, -1]))

    def test_overlap_resolution_rot90(self):
        A = self._s(4, 4); B = self._s(1, 1)
        C1 = self._s(-2, 8); C2 = self._s(-2, 8)
        bAC = findL1Bisector(A, C1, 200, 200)
        bBC = findL1Bisector(B, C2, 200, 200)
        P = bisectorIntersection(bAC, bBC)
        self.assertTrue(samePoint(P, [-1, 4]))

    def test_overlap_resolution_rot180(self):
        A = self._s(4, -4); B = self._s(1, -1)
        C1 = self._s(8, 2); C2 = self._s(8, 2)
        bAC = findL1Bisector(A, C1, 200, 200)
        bBC = findL1Bisector(B, C2, 200, 200)
        P = bisectorIntersection(bAC, bBC)
        self.assertTrue(samePoint(P, [4, 1]))

    def test_overlap_resolution_rot270(self):
        A = self._s(-4, -4); B = self._s(-1, -1)
        C1 = self._s(2, -8); C2 = self._s(2, -8)
        bAC = findL1Bisector(A, C1, 200, 200)
        bBC = findL1Bisector(B, C2, 200, 200)
        P = bisectorIntersection(bAC, bBC)
        self.assertTrue(samePoint(P, [1, -4]))

    def test_normal_shared_site(self):
        A = self._s(4, 6); B = self._s(3, 10)
        C1 = self._s(10, 6); C2 = self._s(10, 6)
        bAC = findL1Bisector(A, C1, 30, 30)
        bBC = findL1Bisector(B, C2, 30, 30)
        P = bisectorIntersection(bAC, bBC)
        self.assertIsInstance(P, list)
        self.assertEqual(len(P), 2)


class TestGenerateVoronoiPoints(unittest.TestCase):

    def test_basic(self):
        pts = [[4, 6], [3, 10], [10, 6], [1, 2]]
        result = generateVoronoiPoints(pts, 30, 30, distance)
        self.assertEqual(len(result), 900)


class TestCleanData(unittest.TestCase):

    def test_nudge(self):
        data = [[4, 6], [6, 4]]
        result = cleanData(data)
        self.assertNotEqual(data[1], [6, 4])


class TestGenerateL1Voronoi(unittest.TestCase):

    def test_four_sites(self):
        sites = [[4, 6], [3, 10], [10, 6], [1, 2]]
        result = generateL1Voronoi(sites, 30, 30, nudgeData=True)
        self.assertEqual(len(result), 4)
        for site in result:
            self.assertIn('d', site)
            self.assertIn('neighbors', site)
            self.assertIn('polygonPoints', site)
            self.assertGreaterEqual(len(site['polygonPoints']), 3)

    def test_16_random_sites(self):
        random.seed(42)
        sites = [[int(random.random() * 400), int(random.random() * 400)] for _ in range(16)]
        result = generateL1Voronoi(sites, 400, 400, nudgeData=True)
        self.assertEqual(len(result), 16)
        for site in result:
            self.assertGreaterEqual(len(site['neighbors']), 1)


class TestRandomStress(unittest.TestCase):

    @staticmethod
    def _random_normal(sharpness):
        return sum(random.random() for _ in range(sharpness)) / sharpness

    def _gen_sites(self, seed, n):
        random.seed(seed)
        return [[int(self._random_normal(2) * 400),
                 int(self._random_normal(2) * 400)] for _ in range(n)]

    def test_89_sites(self):
        sites = self._gen_sites(1, 89)
        result = generateL1Voronoi(sites, 400, 400, nudgeData=False)
        self.assertEqual(len(result), 89)
        for site in result:
            self.assertGreaterEqual(len(site['neighbors']), 1)
            self.assertGreaterEqual(len(site['polygonPoints']), 3)
            self.assertIn('d', site)
            self.assertTrue(site['d'].startswith('M '))
            self.assertTrue(site['d'].endswith(' Z'))


if __name__ == '__main__':
    unittest.main()
