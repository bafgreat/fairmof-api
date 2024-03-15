import os
import requests
import json
from ase import Atoms
import pint
ureg = pint.UnitRegistry()
import numpy as np


def convert_system_to_atoms(system):
    '''
    Convert the system to an ase Atoms object
    '''
    cell = np.array(system['lattice_vectors']) * ureg.meter
    if cell is not None:
        cell = cell.to(ureg.angstrom).magnitude
    if len(system['positions']) > 0:
        coords = np.array(system['positions']) * ureg.meter
        coords = coords.to(ureg.angstrom).magnitude
    return Atoms(positions=coords, numbers=system['species'], cell=cell, pbc=system['periodic'])


def download_archive(refcode, result_folder="FairMOFs", format='cif'):
    '''
    Download ase atoms provided a csd refcode is available

    '''
    url = 'https://nomad-lab.eu/prod/v1/api/v1/entries/archive/query'
    headers = {
        'accept': 'application/zip',
        'Content-Type': 'application/json'
    }
    data = {
        "owner": "public",
        "query": {
            "mainfile": f"MOFData/{refcode}_fair_op/{refcode}_fair_op.out"
        }
    }

    response = requests.post(url, headers=headers, data=json.dumps(data))

    if response.status_code == 200:
        response_dict = {
            "status_code": response.status_code,
            "content": json.loads(response.content.decode('utf-8'))
        }
        if not os.path.exists(result_folder):
            os.makedirs(result_folder)
        label = response_dict['content']['data'][0]['archive']['results']['material']['topology'][1]['label']
        print (response_dict['content']['data'][0]['archive']['results']['material']['topology'][3])
        system = response_dict['content']['data'][0]['archive']['run'][0]['system'][-1]['atoms']

        convert_system_to_atoms(system).write(f'{result_folder}/{refcode}_fair_op.{format}')





response_data = download_archive("EDUSIF")
