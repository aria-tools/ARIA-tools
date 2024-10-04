
import os
import pickle
import json
from datetime import datetime
import shapely

import ARIAtools
import ARIAtools.util.shp

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

        self.config_name = os.path.join(self.log_dir, 'config.json')
        if not os.path.exists(self.config_name):
            with open(self.config_name, 'w') as config_file:
                json.dump({}, config_file)

        self.file_list_name = os.path.join(self.log_dir, 'files.json')
        if not os.path.exists(self.file_list_name):
            with open(self.file_list_name, 'w') as list_file:
                json.dump({}, list_file)

        self.extracted_files_name = os.path.join(self.log_dir,
                                                 'extracted_files.json')
        if not os.path.exists(self.extracted_files_name):
            with open(self.extracted_files_name, 'w') as extr_file:
                json.dump({}, extr_file)


    def load(self):
        """
        """
        with open(self.log_name, 'rb') as log_file:
            log_data = pickle.load(log_file)

        return log_data

    def update(self, atr_name, atr_value):
        """
        """
        # Recall existing log data
        log_data = self.load()

        if self.verbose:
            print(f'Writing {atr_name} to log file')

        # Write new log data
        log_data[atr_name] = atr_value
        with open(self.log_name, 'wb') as log_file:
            pickle.dump(log_data, log_file)

        # Check if attributes should be written to JSON file
        config_params = ['aria_version', 'run_time', 'aria_routine',
                         'input_params', 'workdir', 'bbox']
        if atr_name in config_params:
            self.update_configs(atr_name, atr_value)

        if atr_name == 'files':
            self.write_file_list(atr_value)

        if atr_name == 'extracted_files':
            self.write_extracted_files(atr_value)

    def update_configs(self, atr_name, atr_value):
        """
        """
        # Convert shapely polygon to WKT
        if type(atr_value) == shapely.Polygon:
            atr_value = atr_value.wkt

        # Recall existing config data
        with open(self.config_name, 'r') as config_file:
            config_data = json.load(config_file)

        # Write new data
        config_data[atr_name] = atr_value
        with open(self.config_name, 'w') as config_file:
            json.dump(config_data, config_file)

    def write_file_list(self, files):
        """
        """
        # Write new data
        with open(self.file_list_name, 'w') as list_file:
            json.dump({'files': files}, list_file)

    def write_extracted_files(self, extracted_files):
        """
        """
        extr_dict = {'extracted_files': extracted_files,
                     'nb_extracted': len(extracted_files)}
        with open(self.extracted_files_name, 'w') as extr_file:
            json.dump(extr_dict, extr_file)


    def determine_update_mode(self, outname):
        """
        """
        # Recall log data
        log_data = self.load()

        # Compare shape areas
        if 'exist_prods_TOTbbox' in log_data.keys():
            exist_shp = log_data['exist_prods_TOTbbox']
            exist_area = ARIAtools.util.shp.shp_area(exist_shp, projection=int(log_data['projection']))
            current_shp = ARIAtools.util.shp.open_shp(log_data['prods_TOTbbox'])
            overlap_shp = exist_shp.intersection(current_shp)
            ARIAtools.util.shp.save_shp('xx.json', overlap_shp, projection=int(log_data['projection']))
            overlap_area = ARIAtools.util.shp.shp_area(overlap_shp, projection=int(log_data['projection']))
            area_ratio = overlap_area / exist_area
        else:
            area_ratio = 0

        if area_ratio == 1.0 and os.path.exists(outname+'.vrt'):
            update_mode = 'skip'
        elif area_ratio > 0.99 and os.path.exists(outname+'.vrt'):
            update_mode = 'crop_only'
        else:
            update_mode = 'full_extract'

        return update_mode
