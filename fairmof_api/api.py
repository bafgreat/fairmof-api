import os
import json
import re
import glob
import requests
import numpy as np
from ase import Atoms
import pint
from  mofstructure import mofdeconstructor
from fairmof_api import filetyper
ureg = pint.UnitRegistry()


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

def find_secondary_building_units(system_topology, system_atoms):
    """
    Find the secondary building units in the system topology

    Parameters
    ----------
    system_topology : dict
        The system topology

    Returns
    -------
    secondary_building_units : list
        A list of secondary building units
    """
    sbu_data = {}
    for i in range(len(system_topology)):
        tmp_sbu = {}
        label = system_topology[i]['label']
        if check_sbu_or_ligand(label, 'metal_sbu') and system_topology[i]['structural_type'] == 'molecule':
            sbu_indices = system_topology[i]['indices'][0]
            ase_atom = mofdeconstructor.wrap_systems_in_unit_cell(system_atoms[sbu_indices])
            tmp_sbu['ase_atom'] = ase_atom
            tmp_sbu['sbu_type'] = system_topology[i]['sbu_type']
            tmp_sbu['sbu_coordination_number']=system_topology[i]['sbu_coordination_number']
            tmp_sbu['label'] = label
            sbu_data[label] = tmp_sbu
        if check_sbu_or_ligand(label, 'organic_sbu') and system_topology[i]['structural_type'] == 'molecule':
            sbu_indices = system_topology[i]['indices'][0]
            ase_atom = mofdeconstructor.wrap_systems_in_unit_cell(system_atoms[sbu_indices])
            tmp_sbu['ase_atom'] = ase_atom
            tmp_sbu['sbu_coordination_number']=system_topology[i]['sbu_coordination_number']
            tmp_sbu['label'] = label
            sbu_data[label] = tmp_sbu


        if check_sbu_or_ligand(label, 'ligand') and system_topology[i]['structural_type'] == 'molecule':
            sbu_indices = system_topology[i]['indices'][0]
            ase_atom = mofdeconstructor.wrap_systems_in_unit_cell(system_atoms[sbu_indices])
            tmp_sbu['ase_atom'] = ase_atom
            tmp_sbu['label'] = label
            sbu_data[label] = tmp_sbu
    return sbu_data


def check_sbu_or_ligand(string, sbu_or_ligand):
    """
    Check if the string is a sbu or ligand

    Parameters
    ----------
    """
    pattern = rf'{sbu_or_ligand}'
    match = re.search(pattern, string)
    if match:
        return True
    else:
        return False



def extract_mof_and_properties(system_atom, system_topology):
    """
    Extract the mof and properties from the system topology

    Parameters
    ----------
    systm_atom : ase Atoms object
        The ase Atoms object
    system_topology : dict
        The system topology

    Returns
    -------
    mof : ase Atoms object
        The mof
    properties : dict
        The properties
    """
    properties = {}
    mof_indices = system_topology['indices'][0]
    void_fraction = system_topology['void_fraction']
    pore_limiting_diameter = system_topology['pore_limiting_diameter']*ureg.meter
    largest_cavity_diameter = system_topology['largest_cavity_diameter']*ureg.meter
    largest_included_sphere_along_free_sphere_path = system_topology['largest_included_sphere_along_free_sphere_path']*ureg.meter
    accessible_surface_area = system_topology['accessible_surface_area']*ureg.meter**2
    accessible_volume = system_topology['accessible_volume']*ureg.meter**3
    n_channels = system_topology['n_channels']

    properties['void_fraction'] = void_fraction
    properties['n_channels'] = n_channels
    properties['pore_limiting_diameter_A'] = pore_limiting_diameter.to(ureg.angstrom).magnitude
    properties['largest_cavity_diameter_A'] = largest_cavity_diameter.to(ureg.angstrom).magnitude
    properties['largest_included_sphere_along_free_sphere_path_A'] = largest_included_sphere_along_free_sphere_path.to(ureg.angstrom).magnitude
    properties['accessible_surface_area_A^2'] = accessible_surface_area.to(ureg.angstrom**2).magnitude
    properties['accessible_volume_A^3'] = accessible_volume.to(ureg.angstrom**3).magnitude

    mof_atom = system_atom[mof_indices]
    return mof_atom, properties


def download_archive(refcode, result_folder="FAIR-MOFs", extension='cif'):
    '''
    Download ase atoms provided a csd refcode is available

    '''
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    mof_cif_path = os.path.join(result_folder, 'mofs_recognised_by_nomad')
    if not os.path.exists(mof_cif_path):
        os.makedirs(mof_cif_path)
    non_mof_cif_path = os.path.join(result_folder, 'mofs_not_recognised_by_nomad')
    if not os.path.exists(non_mof_cif_path):
        os.makedirs(non_mof_cif_path)

    mof_properties_filename = os.path.join(result_folder, 'mof_properties.json')
    mof_sbu_filename = os.path.join(result_folder, 'mof_sbu.json')
    # Create json files name
    if isinstance(mof_properties_filename, str) and os.path.exists(mof_properties_filename):
        tmp_dic = filetyper.load_data(mof_properties_filename)
    if isinstance(mof_sbu_filename, str) and os.path.exists(mof_sbu_filename):
        sbu_data = filetyper.load_data(mof_sbu_filename)


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

        print(refcode)

        try:
            system = response_dict['content']['data'][0]['archive']['run'][0]['system'][-1]['atoms']
            system_atom = convert_system_to_atoms(system)

            system_topology = response_dict['content']['data'][0]['archive']['results']['material']['topology']

            if system_topology[1]['label'] == 'MOF':
                tmp_dic = {}
                sbu_data = {}

                system_mof = system_topology[1]
                mof_atom, mof_properties = extract_mof_and_properties(system_atom, system_mof)
                tmp_dic[f'{refcode}_fair_op'] = mof_properties
                mof_atom.write(f'{mof_cif_path}/{refcode}_fair_op.{extension}')

                filetyper.append_json(tmp_dic, mof_properties_filename)

                sbu_data[f'{refcode}_fair_op'] = find_secondary_building_units(system_topology, system_atom)
                filetyper.append_json_atom(sbu_data, mof_sbu_filename)
            else:
                system_atom.write(f'{non_mof_cif_path}/{refcode}_fair_op.{extension}')
        except Exception:
            pass


def atoms_from_entry_id(entry_id_json_file, result_folder="FAIR-MOFs", extension='cif'):
    """
    Download ase atoms provided a csd refcode is available

    Parameters
    ----------
    entry_id_json_file : str
    result_folder : str, optional
    extension : str, optional

    """

    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    mof_cif_path = os.path.join(result_folder, 'mofs_recognised_by_nomad')
    if not os.path.exists(mof_cif_path):
        os.makedirs(mof_cif_path)
    non_mof_cif_path = os.path.join(result_folder, 'mofs_not_recognised_by_nomad')
    if not os.path.exists(non_mof_cif_path):
        os.makedirs(non_mof_cif_path)

    # mof_properties_filename = os.path.join(result_folder, 'mof_properties.json')
    # mof_sbu_filename = os.path.join(result_folder, 'mof_sbu.json')
    # # Create json files name
    # if isinstance(mof_properties_filename, str) and os.path.exists(mof_properties_filename):
    #     tmp_dic = filetyper.load_data(mof_properties_filename)
    # if isinstance(mof_sbu_filename, str) and os.path.exists(mof_sbu_filename):
    #     sbu_data = filetyper.load_data(mof_sbu_filename)


    entry_id = filetyper.load_data(entry_id_json_file)

    mainfile = entry_id['archive']['metadata']['mainfile']
    refcode = mainfile.split('/')[-1].split('.')[0]

    print(refcode)
    # print (entry_id['archive'].keys())
    # print (entry_id['archive'].keys())

    archive = entry_id['archive']

    system = archive['run'][0]['system'][-1]['atoms']
    # print (system)

    system_atom = convert_system_to_atoms(system)
    # print (system_atom)

    system_topology = archive['results']['material']['topology']

    if len(system_topology)>1:
        if system_topology[1]['label'] == 'MOF':
            tmp_dic = {}
            sbu_data = {}

            system_mof = system_topology[1]
            mof_atom, mof_properties = extract_mof_and_properties(system_atom, system_mof)
            tmp_dic[refcode] = mof_properties
            mof_atom.write(f'{mof_cif_path}/{refcode}.{extension}')

            # filetyper.append_json(tmp_dic, mof_properties_filename)

            # sbu_data[refcode] = find_secondary_building_units(system_topology, system_atom)
            # filetyper.append_json_atom(sbu_data, mof_sbu_filename)
        else:
            system_atom.write(f'{non_mof_cif_path}/{refcode}.{extension}')
    else:
        system_atom.write(f'{non_mof_cif_path}/{refcode}.{extension}')
    # except Exception:
    #     pass
