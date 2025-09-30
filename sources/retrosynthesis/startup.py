import os
from flask import request
from sources.retrosynthesis.classes import AiZynthArgs
from aizynthfinder import aizynthfinder
from aizynthfinder.interfaces import aizynthcli


def make_config():
    BASEDIR = os.path.dirname(__file__)

    stock_request = request.args.get('stock')

    stock_file_dict = {
        'zinc': 'zinc_stock_17_04_20.hdf5',
        'overlay': 'stock_additions_2025_08.hdf5',
        'emols': 'zinc_and_emol_inchi_key.bloom',
        'alcohols': 'dummy_alcohols.hdf5',
        'naturals': 'np_08_25.hdf5',
        'non_iso_naturals': 'ninp_08_25.hdf5'}


    AIZYNTH = {
        'expansion': {
            'full': {
                'type': 'template-based',
                'model': os.path.join(BASEDIR, 'config_files', 'uspto_model.onnx'),
                'template': os.path.join(BASEDIR, 'config_files', 'uspto_templates.csv.gz')
            }
        },

        'stock': {
            'zinc': os.path.join(BASEDIR, 'config_files', 'zinc_stock_17_04_20.hdf5'),
            'overlay': os.path.join(BASEDIR, 'config_files', 'stock_additions_2025_08.hdf5'),
            'alcohols': os.path.join(BASEDIR, 'config_files', 'dummy_alcohols.hdf5'),
            'naturals': os.path.join(BASEDIR, 'config_files', 'np_08_25.hdf5'),
            'non_iso_naturals': os.path.join(BASEDIR, 'config_files', 'ninp_08_25.hdf5'),
            'paroutes_n1': os.path.join(BASEDIR, 'config_files', 'n1_stock.hdf5'),
            'paroutes_n5': os.path.join(BASEDIR, 'config_files', 'n5_stock.hdf5'),
            'askos_stock': os.path.join(BASEDIR, 'config_files', 'askcos_stock.hdf5')
            },

        #'stock': {'bloom': os.path.join(BASEDIR, 'config_files', 'zinc_and_emol_inchi_key.bloom')},
        'config_file': os.path.join(BASEDIR, 'config_files', 'aizynthfinder_config.yml'),
        'properties': {
            'max_transforms': 10,
            'time_limit': 3600,
            'iteration_limit': 500,
        }
    }

    print("1")
    print(AIZYNTH)
    aizynth_config_dic = AIZYNTH
    print(2)
    args = AiZynthArgs("placeholder",  aizynth_config_dic['expansion'],
                       aizynth_config_dic['stock'])
    args.stock = stock_request or 'zinc,overlay'
    print(3)
    finder = aizynthfinder.AiZynthFinder(configdict=aizynth_config_dic)
    print(4)
    aizynthcli._select_stocks(finder, args)
    print(5)
    finder.expansion_policy.select(args.policy or finder.expansion_policy.items[0])
    print(6)
    finder.filter_policy.select(args.filter)
    print(7)
    return finder
