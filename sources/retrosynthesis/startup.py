import os
from flask import request
from sources.retrosynthesis.classes import AiZynthArgs
from aizynthfinder import aizynthfinder
from aizynthfinder.interfaces import aizynthcli


def make_config():
    BASEDIR = os.path.dirname(__file__)

    stock_request = request.args.get('stock')

    stock_file_dict = {'emolecules': 'emolecules.hdf5',
                       'zinc': 'zinc_stock_17_04_20.hdf5'}

    chosen_stock_file = stock_file_dict.get(stock_request)

    AIZYNTH = {
        'expansion': {
            'full': {
                'type': 'template-based',
                'model': os.path.join(BASEDIR, 'config_files', 'uspto_model.onnx'),
                'template': os.path.join(BASEDIR, 'config_files', 'uspto_templates.csv.gz')
            }
        },

        'stock': {'bloom': os.path.join(BASEDIR, 'config_files', 'zinc_and_emol_inchi_key.bloom')}
        ,
        'config_file': os.path.join(BASEDIR, 'config_files', 'aizynthfinder_config.yml'),
        'search': {
            'max_transforms': 10,
            'time_limit': 3600,
            'iteration_limit': 500,
            # 'excluded_stock': {'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)Cl)ccc3OCC)nc12',
            #                    'CCCc1nn(C)c(C(N)=O)c1NC(=O)c1cc(S(=O)(=O)N2CCN(C)CC2)ccc1OCC',
            #                    'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3O)nc12',
            #                    'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)O)ccc3OCC)nc12',
            #                    'CCCOc1ccc(S(=O)(=O)N2CCN(C)CC2)cc1-c1nc2c(CCC)nn(C)c2c(=O)[nH]1',
            #                    'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)Cl)ccc3OCC)nc12',
            #                    'CCOc1ccc(S(=O)(=O)N2CCN(C)CC2)cc1C=O',
            #                    'CCCc1nn(C)c(C(N)=O)c1NC(=O)c1cc(S(=O)(=O)N2CCN(C)CC2)ccc1OCC',
            #                    'CCOc1ccc(S(=O)(=O)N2CCN(C)CC2)cc1C(=O)O'
            #                    'CCOc1ccc(S(=O)(=O)N2CCN(C)CC2)cc1C#N'
            #                    'CCOc1ccc(S(=O)(=O)N2CCN(C)CC2)cc1C(=O)O'
            #                    # 'CCOc1ccc(S(=O)(=O)N2CCN(C)CC2)cc1C#N'
            #                    },

        }
    }


    print("1")
    print(AIZYNTH)
    aizynth_config_dic = AIZYNTH
    print(2)
    # initiate object containing all required arguments
    args = AiZynthArgs("placeholder",  aizynth_config_dic['expansion'],
                      aizynth_config_dic['stock'])
    print(3)
    # AiZynthFinder object contains results data
    finder = aizynthfinder.AiZynthFinder(configdict=aizynth_config_dic)
    print(4)
    # set up stocks, policies, then start single smiles process
    aizynthcli._select_stocks(finder, args)
    print(5)
    finder.expansion_policy.select(args.policy or finder.expansion_policy.items[0])
    print(6)
    finder.filter_policy.select(args.filter)
    print(7)
    return finder
