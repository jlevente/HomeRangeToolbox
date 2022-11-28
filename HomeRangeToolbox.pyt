# -*- coding: utf-8 -*-

"""ArcGIS Pro toolbox to estimate home ranges.

@Author: Levente Juhász

Copyright (c) 2022 Levente Juhász.
For more, see LICENSE
"""

import arcpy
import os
import re

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Home Range Analysis Toolbox"
        self.alias = "HomeRangeAnalysis"

        # List of tool classes associated with this toolbox
        self.tools = [HomeRangeKDE, HomeRangeKDE_Batch, HomeRangeMCP]

class HomeRangeKDE(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Home Range Estimation using KDE"
        self.description = '''Home Range Estimation using Kernel Density Estimators. This tool calculates home
        ranges for observations stored in a Feature Class.
        '''
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        observations = arcpy.Parameter(
            displayName="Input observations",
            name="observations",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input"
        )

        observations.filter.list = ['Point', 'MultiPoint']

        barrier_features = arcpy.Parameter(
            displayName="Barrier Features",
            name="barrier_features",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input"
        )

        search_radius = arcpy.Parameter(
            displayName="Search radius",
            name="search_radius",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )

        out_cell_size = arcpy.Parameter(
            displayName="Output cell size",
            name="out_cell_size",
            datatype="analysis_cell_size",
            parameterType="Required",
            direction="Output"
        )

        out_folder = arcpy.Parameter(
            displayName="Output folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Output"
        )

        params = [observations, barrier_features, search_radius, out_cell_size, out_folder]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Load core functionality
        processor = HomeRangeCalc()

        arcpy.CheckOutExtension('spatial')
        arcpy.env.overwriteOutput = True

        param_values = {}
        params = {}
        for param in parameters:
            param_values[param.name] = param.valueAsText
            params[param.name] = param

        obs_path = arcpy.Describe(params['observations']).catalogPath

        param_values['out_folder'] = param_values['out_folder'].replace('\\', '/')
        os.makedirs(param_values['out_folder'], exist_ok=True)
        arcpy.env.workspace = param_values['out_folder']
        arcpy.AddMessage("Outputs will be saved in %s" % (param_values['out_folder']))

        if params['barrier_features'].value:
            barrier_path = arcpy.Describe(params['barrier_features'].valueAsText).catalogPath
            processor.compute_utilization(arcpy, point_fc=obs_path, suffix=None, barrier_features=barrier_path,\
                                            cell_size=params['out_cell_size'].valueAsText, search_radius=params['search_radius'].valueAsText)
        else:
            print('e')
            processor.compute_utilization(arcpy, point_fc=obs_path, suffix=None, barrier_features=None, \
                                            cell_size=params['out_cell_size'].valueAsText, search_radius=params['search_radius'].valueAsText)

        arcpy.AddMessage("All done.")
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class HomeRangeKDE_Batch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Batch Home Range Estimation using KDE"
        self.description = '''Batch Home Range Estimation using Kernel Density Estimators. This tool calculates home
        ranges for multiple individuals if observations are stored in one Feature Class.
        '''
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        observations = arcpy.Parameter(
            displayName="Input observations",
            name="observations",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input"
        )

        observations.filter.list = ['Point', 'MultiPoint']

        animal_id = arcpy.Parameter(
            displayName="Animal ID field",
            name="animal_id",
            datatype="Field",
            parameterType="Required",
            direction="Input"
        )

        animal_id.parameterDependencies = [observations.name]

        barrier_features = arcpy.Parameter(
            displayName="Barrier Features",
            name="barrier_features",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input"
        )

        home_cutoff = arcpy.Parameter(
            displayName="Home range percentage",
            name="home_cutoff",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )

        search_radius = arcpy.Parameter(
            displayName="Search radius",
            name="search_radius",
            datatype="Double",
            parameterType="Optional",
            direction="Input"
        )

        out_cell_size = arcpy.Parameter(
            displayName="Output cell size",
            name="out_cell_size",
            datatype="analysis_cell_size",
            parameterType="Required",
            direction="Output"
        )

        out_folder = arcpy.Parameter(
            displayName="Output folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Output"
        )

        params = [observations, animal_id, barrier_features, search_radius, out_cell_size, out_folder]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""

        # Load core functionality
        processor = HomeRangeCalc()

        arcpy.CheckOutExtension('spatial')
        arcpy.env.overwriteOutput = True

        param_values = {}
        params = {}
        for param in parameters:
            param_values[param.name] = param.valueAsText
            params[param.name] = param

        obs_path = arcpy.Describe(params['observations']).catalogPath

        param_values['out_folder'] = param_values['out_folder'].replace('\\', '/')
        os.makedirs(param_values['out_folder'], exist_ok=True)
        arcpy.env.workspace = param_values['out_folder']
        arcpy.AddMessage("Outputs will be saved in %s" % (param_values['out_folder']))

        individuals = []

        # Get absolute path of observations
        with arcpy.da.SearchCursor(obs_path, param_values['animal_id']) as cursor:
            for row in cursor:
                animal = row[0]
                if animal not in individuals:
                    individuals.append(animal)
        
        # Loop over each individual and calculate KDE
        arcpy.AddMessage("Found %s individuals." % len(individuals))
        i = 0
        for animal in individuals:
            i += 1
            if i % 4 == 0:
                arcpy.AddMessage("%s/%s individuals done." % (i, len(individuals)))
                break

            tmp_points = r'animal.shp'
                      
            where_txt = '"%s" = \'%s\'' % (param_values['animal_id'], animal)
           
            arcpy.Select_analysis(obs_path, tmp_points, where_clause=where_txt)
            tmp_points_path = arcpy.Describe(tmp_points).catalogPath


            # Check and set up parameters
            if params['barrier_features'].value:
                barrier_path = arcpy.Describe(params['barrier_features'].valueAsText).catalogPath
            else:
                barrier_path = None

            out_cell_size = params['out_cell_size'].valueAsText

            if params['search_radius'].value:
                bandwidth = params['search_radius'].valueAsText
            else:
                coord_cursor = arcpy.da.SearchCursor(tmp_points_path, ["SHAPE@XY"])
                bandwidth = processor.calculate_optimized_bandwidth(coord_cursor)
                arcpy.AddMessage('Optimal bandwidth: %s' % bandwidth)

            processor.compute_utilization(arcpy, point_fc=tmp_points_path, suffix=animal, barrier_features=barrier_path, \
                                            cell_size=out_cell_size, search_radius=bandwidth)


        arcpy.AddMessage("All done.")
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class HomeRangeMCP(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Home Range Estimation using MCP"
        self.description = "Home Range Estimation using Minimum Convex Polygons"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        params = None
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

# Class that contains the functionality
class HomeRangeCalc(object):
    def compute_utilization(self, arcpy, point_fc, suffix, barrier_features, cell_size, search_radius):
        # Build KDE raster
        if barrier_features:
            raster = arcpy.sa.KernelDensity(point_fc, None, cell_size, search_radius, in_barriers=barrier_path)
        else:
            raster = arcpy.sa.KernelDensity(point_fc, None, cell_size, search_radius)

        if suffix:
            out_raster_name = "kde_%s.tif" % suffix
        else:
            out_raster_name = "kde.tif"

        raster.save(out_raster_name)
        arcpy.management.CalculateStatistics(out_raster_name)
        out_raster_path = arcpy.Describe(out_raster_name).catalogPath

        # Extract home ranges and core areas
        tmp_raster_points =  r'raster_values.shp'
        value_list = []

        arcpy.sa.ExtractValuesToPoints(point_fc, out_raster_name, tmp_raster_points)
        with arcpy.da.SearchCursor(tmp_raster_points, 'RASTERVALU') as cursor:
            for val in cursor:
                value_list.append(val[0])

        value_list.sort(reverse=True)
        num_records = len(value_list)

        home_cutoff = int(num_records * 0.5)
        core_cutoff = int(num_records * 0.95)

        home_cutoff_value = value_list[home_cutoff - 1]
        core_cutoff_value = value_list[core_cutoff - 1]

        max_kernel_value = arcpy.GetRasterProperties_management(out_raster_name, 'MAXIMUM')

        # Reclass expression format: min_value max_value RECLASS_VALUE (separated by ;)
        reclass_expression_home = '0 %s NODATA; %s %s %s' % (home_cutoff_value, home_cutoff_value, max_kernel_value, 1)
        reclass_expression_core = '0 %s NODATA; %s %s %s' % (core_cutoff_value, core_cutoff_value, max_kernel_value, 1)


        home_raster = arcpy.sa.Reclassify(out_raster_path, 'Value', reclass_expression_home, 'DATA')
        core_raster = arcpy.sa.Reclassify(out_raster_path, 'Value', reclass_expression_core, 'DATA')

        if suffix:
            out_home_name = r"home_range_%s" % suffix
            out_core_name = r"core_%s" % suffix
        else:
            out_home_name = r'home_range'
            out_core_name = r'core'
        home_raster.save(out_home_name + '.tif')
        core_raster.save(out_core_name + '.tif')

        arcpy.AddMessage(out_home_name + '.shp')
        arcpy.RasterToPolygon_conversion(home_raster, out_home_name + '.shp', 'SIMPLIFY', 'Value')
        arcpy.RasterToPolygon_conversion(core_raster, out_core_name + '.shp', 'SIMPLIFY', 'Value')

        arcpy.Delete_management(out_raster_path)
        arcpy.Delete_management(tmp_raster_points)
        # Do not delete point feature class if not called from within a loop (i.e. suffix is None)
        #if suffix:
        #    arcpy.Delete_management(point_fc)

    # This function extracts lists of x and y coordinates from a search cursor,
    # then returns and optimized bandwidth value for those observations
    # based on XXX
    def calculate_reference_bandwidth(self, search_cursor):
        import math
        import statistics

        x_list = []
        y_list = []

        for feature in search_cursor:
            x, y = feature[0]
            x_list.append(x)
            y_list.append(y)

        mean_x = sum(x_list) / len(x_list)
        mean_y = sum(y_list) / len(y_list)
        n = len(x_list)

        distances = []
        for i in range(len(x_list)):
            sq_dist = ((x_list[i] - mean_x) ** 2)  + ((y_list[i] - mean_y) ** 2)
            distances.append(sq_dist)

        st_dist = math.sqrt(sum(distances) / n)
        d_m = math.sqrt(1/math.log(2)) * math.sqrt(statistics.median(distances))


        if st_dist < d_m:
            A = st_dist
        else:
            A = d_m

        h = 0.9 * A * len(x_list) ** (-1/5)

        return h

    def calculate_lscv_bandwidth(self, search_cursor):
        x_list = []
        y_list = []

        for feature in search_cursor:
            x, y = feature[0]
            x_list.append(x)
            y_list.append(y)

        h_ref = self.calculate_reference_bandwidth(search_cursor)


# def main():
#     tbx = Toolbox()

#     test_tool = HomeRangeKDE_Batch()
#     test_tool.execute(test_tool.getParameterInfo(), None)

# if __name__ == '__main__':
#     main()