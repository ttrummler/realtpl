import numpy as np
import os
import yaml


class NasaCoefficients:

    def __init__(self, name: str, data: list):
        self.name = name
        self.temp_bin_edges, self.n_coeff, self.coeff = self._parse_data(data)

    def _parse_data(self, data):
        data = sorted(data, key=lambda x: x['temp_start'])
        temp_start = np.array([x['temp_start'] for x in data])
        temp_end = np.array([x['temp_end'] for x in data])
        if not np.all(temp_start[:-1] < temp_start[1:]) \
           or not np.all(temp_end[:-1] < temp_end[1:]) \
           or not np.all(np.equal(temp_start[1:], temp_end[:-1])):
            raise RuntimeError(
                f'Temperature ranges provided for {self.name} '
                'are not consistent, please fix this.'
            )
        temp_bin_edges = np.append(temp_start, temp_end[-1])

        n_coeff = set([len(x['coeff']) for x in data])
        if len(n_coeff) == 1:
            n_coeff = n_coeff.pop()
        else:
            raise RuntimeError(
                f'Inconsistent coefficient data provided for {self.name} '
                'please fix this.'
            )

        coeff = np.empty([n_coeff, len(data)])
        for i, x in enumerate(data):
            coeff[:, i] = x['coeff']

        return temp_bin_edges, n_coeff, coeff

    def get_temp_range(self):
        return (self.temp_bin_edges[0], self.temp_bin_edges[-1])

    def get_coeff(self, idx: int, temp: np.array):
        return self.coeff[idx, :][
                np.searchsorted(self.temp_bin_edges[1:-1], temp, side='left')
        ]

    @classmethod
    def from_name_and_coeff(cls, name: str, n_coeff: int):
        if n_coeff == 7:
            with open(os.path.join(os.path.dirname(__file__),
                                   'nasa_7.yaml'), 'r') as stream:
                data_nasa_all = yaml.safe_load(stream)
        elif n_coeff == 9:
            with open(os.path.join(os.path.dirname(__file__),
                                   'nasa_9.yaml'), 'r') as stream:
                data_nasa_all = yaml.safe_load(stream)
        else:
            raise ValueError(f'Unknown NASA coefficient number: {n_coeff}')

        data = data_nasa_all.get(name)
        if not data:
            raise ValueError(f'No nasa data found for {name}!\n For the'
                             f' computation you have to provide the values '
                             f'in nasa_{n_coeff}.yaml.\n Useful sources'
                             f' therefore are provided in the header of the '
                             f'nasa_X.yaml files.')
        return cls(name, data)
