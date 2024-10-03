
import os
import pickle
import json
from datetime import datetime
from shapely import from_wkt

import ARIAtools

class RunLog:
    """
    """
    def __init__(self, workdir, verbose=None):
        """
        """
        # Record parameters
        self.workdir = os.path.abspath(workdir)

        self.verbose = verbose

        # Automatically establish log directory and files
        self.log_dir = os.path.join(self.workdir, 'config')
        os.makedirs(self.log_dir, exist_ok=True)

        self.log_name = os.path.join(self.log_dir, 'PICKLE.pkl')
        if not os.path.exists(self.log_name):
            with open(self.log_name, 'wb') as log_file:
                pickle.dump({}, log_file)

    def load(self):
        """
        """
        if self.verbose:
            print('Reading log file')

        with open(self.log_name, 'rb') as log_file:
            log_data = pickle.load(log_file)

        return log_data

    def write(self, atr_name, atr_value):
        """
        """
        # Recall existing log data
        log_data = self.load()

        # Write new log data
        log_data[atr_name] = atr_value
        with open(self.log_name, 'wb') as log_file:
            pickle.dump(log_data, log_file)

    def write_gunws(self, files):
        """
        """
        if self.verbose:
            print('Logging list of GUNWs')

        self.write('files', files)


    # def write_gunws(self):
    #     """
    #     """
    #     return None


    # def __dump_pickle__(self):
    #     """
    #     Dump information to machine-readable-only pickle file.
    #     """
    #     pickle_name = os.path.join(self.log_dir, 'PICKLE.pkl')
    #     if self.verbose:
    #         print(f"PICKLE file: {pickle_name}")


    # def __dump_jsons__(self):
    #     """
    #     Dump relevant information to human-readable json file.
    #     """
    #     # Run ID based on date/time
    #     run = datetime.now().strftime('%Y%m%d-%H%M%S')

    #     # Save info to files with common run ID
    #     self.dump_gunw_list(run)

    #     self.dump_config_params(run)

    #     self.dump_extracted_product_info(run)
    #     exit()

    # def dump_gunw_list(self, run):
    #     """
    #     """
    #     # GUNW list file name
    #     outname = f"gunw_list.json"
    #     outname = os.path.join(self.log_dir, outname)

    # def dump_config_params(self, run):
    #     """
    #     Record:
    #     ARIAtools version
    #     Date/time of run
    #     Nb products to extract
    #     Program (ariaExtract.py or ariaTSsetup.py)
    #     Specific input parameters:
    #         bbox
    #         minimum overlap threshold
    #         gunw version
    #         layers
    #         workdir
    #         lat/lon spacing
    #     """
    #     # Config file name
    #     outname = f"aria_config.json"
    #     outname = os.path.join(self.log_dir, outname)

    #     # Items
    #     config_params = {'aria_version': ARIAtools.__version__,
    #                      'run_time': run}

    #     # Write to file
    #     if self.verbose:
    #         print(f"Writing config params to {outname}")

    #     with open(outname, 'w') as json_file:
    #         json.dump(config_params, json_file)

    # def dump_extracted_product_info(self, run):
    #     """
    #     List of extracted products.
    #     """
    #     # Extracted prod list name
    #     outname = f"extracted_prod_list.json"
    #     outname = os.path.join(self.log_dir, outname)
