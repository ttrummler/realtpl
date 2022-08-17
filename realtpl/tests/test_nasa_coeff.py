import numpy as np
from unittest import TestCase

from realtpl import nasa


class TestCoeff(TestCase):

    def test_simple(self):
        data = [
            {
                'temp_start': 50,
                'temp_end': 100,
                'coeff': [1, 3, 5]
            }
        ]
        ns = nasa.NasaCoefficients('test', data)
        self.assertTupleEqual(ns.get_temp_range(), (50, 100))
        self.assertListEqual(
            list(ns.get_coeff(1, np.array([45, 60, 110]))),
            [3.0, 3.0, 3.0]
        )
        return

    def test_complex(self):
        data = [
            {
                'temp_start': 50,
                'temp_end': 100,
                'coeff': [1, 3, 5]
            },
            {
                'temp_start': 100,
                'temp_end': 150,
                'coeff': [8, 9, 11]
            }
        ]
        ns = nasa.NasaCoefficients('test', data)
        self.assertTupleEqual(ns.get_temp_range(), (50, 150))
        self.assertListEqual(
            list(ns.get_coeff(1, np.array([45, 50, 60, 100, 110, 160]))),
            [3.0, 3.0, 3.0, 3.0, 9.0, 9.0]
        )
        return

    def test_inconsistent_t(self):
        data = [
            {
                'temp_start': 50,
                'temp_end': 100,
                'coeff': [1, 3, 5]
            },
            {
                'temp_start': 50,
                'temp_end': 200,
                'coeff': [1, 3, 5]
            }
        ]
        with self.assertRaises(RuntimeError):
            nasa.NasaCoefficients('test', data)
        return

    def test_inconsistent_coeff(self):
        data = [
            {
                'temp_start': 50,
                'temp_end': 100,
                'coeff': [1, 3, 5]
            },
            {
                'temp_start': 100,
                'temp_end': 200,
                'coeff': [4, 6, 8, 9]
            }
        ]
        with self.assertRaises(RuntimeError):
            nasa.NasaCoefficients('test', data)
        return
