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
        self.tools = [HomeRangeKDE_Batch, HomeRangeMCP]

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
        arcpy.CheckOutExtension('spatial')
        arcpy.env.overwriteOutput = True

        param_values = {}
        params = {}
        for param in parameters:
            param_values[param.name] = param.valueAsText
            params[param.name] = param

        obs_path = arcpy.Describe(params['observations']).catalogPath

        arcpy.AddMessage(arcpy.Describe(parameters[0].valueAsText).name)
        arcpy.AddMessage(parameters)

        param_values['out_folder'] = param_values['out_folder'].replace('\\', '\\\\')
        os.makedirs(param_values['out_folder'], exist_ok=True)
        arcpy.env.workspace = param_values['out_folder']
        arcpy.AddMessage("Outputs will be saved in %s" % (param_values['out_folder']))

        individuals = []

        # Get absolute path of observations
        obs_path = arcpy.Describe(params['observations'].valueAsText).catalogPath
        arcpy.AddMessage(arcpy.Describe(parameters[0]))
        
        with arcpy.da.SearchCursor(obs_path, param_values['animal_id']) as cursor:
            for row in cursor:
                animal = row[0]
                if animal not in individuals:
                    individuals.append(animal)
        
        arcpy.AddMessage("Found %s individuals." % len(individuals))
        i = 0
        for animal in individuals:
            i += 1
            if i % 10 == 0:
                arcpy.AddMessage("%s/%s individuals done." % (i, len(individuals)))
                break

            tmp_points = r'animal.shp'
            tmp_raster_points =  'raster_values.shp'
          
            where_txt = '"%s" = \'%s\'' % (param_values['animal_id'], animal)
           
            arcpy.Select_analysis(obs_path, tmp_points, where_clause=where_txt)
            tmp_points_path = arcpy.Describe(tmp_points).catalogPath

            # Build KDE raster
            if params['barrier_features'].value:
                barrier_path = arcpy.Describe(params['barrier_features'].valueAsText).catalogPath
                raster = arcpy.sa.KernelDensity(tmp_points, None, params['out_cell_size'].value, in_barriers=barrier_path)
            else:
                raster = arcpy.sa.KernelDensity(tmp_points, None, params['out_cell_size'].valueAsText, param_values['search_radius'])

            out_raster_name = "kde_%s.tif" % animal
            raster.save(out_raster_name)
            arcpy.management.CalculateStatistics(out_raster_name)
            out_raster_path = arcpy.Describe(out_raster_name).catalogPath

            # Extract home ranges and core areas
            value_list = []

            arcpy.AddMessage("%s \n %s \n %s" % (tmp_points_path, out_raster_name, tmp_raster_points))
            arcpy.sa.ExtractValuesToPoints(tmp_points_path, out_raster_name, 'C:/Users/Levente Juhasz/projects/tmp/tmpdata/raster_values.shp')
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

            arcpy.AddMessage(reclass_expression_core)

            home_raster = arcpy.sa.Reclassify(out_raster_path, 'Value', reclass_expression_home, 'DATA')
            core_raster = arcpy.sa.Reclassify(out_raster_path, 'Value', reclass_expression_core, 'DATA')


            home_raster.save('home_raster_' + animal + '.tif')
            core_raster.save('core_raster_' + animal + '.tif')

            arcpy.RasterToPolygon_conversion(home_raster, 'home_range_' + animal + '.shp', 'SIMPLIFY', 'Value')
            arcpy.RasterToPolygon_conversion(core_raster, 'core_' + animal + '.shp', 'SIMPLIFY', 'Value')

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

# def main():
#     tbx = Toolbox()

#     test_tool = HomeRangeKDE_Batch()
#     test_tool.execute(test_tool.getParameterInfo(), None)

# if __name__ == '__main__':
#     main()