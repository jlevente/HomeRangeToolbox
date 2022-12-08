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
        """ArcGIS Python Toolbox class."""
        self.label = "Home Range Analysis Toolbox"
        self.alias = "HomeRangeAnalysis"

        # List of tool classes associated with this toolbox
        self.tools = [HomeRangeKDE, HomeRangeKDE_Batch, HomeRangeMCP, HomeRangeMCP_Batch]

class HomeRangeKDE(object):
    """Class to represent Home Range Calculations using KDE in a Python toolbox."""
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
        radius_method = arcpy.Parameter(
                displayName='Bandwidth calculation method', 
                name='radius_method', 
                datatype='GPString', 
                parameterType='Optional', 
                direction="Input", 
                multiValue=False)

        radius_method.filter.type = "ValueList"
        radius_method.filter.list = ['Reference bandwidth (Silverman 1986)', 'Least Squares Cross Validation (Seaman & Powell 1996)', 'Manual']

        search_radius = arcpy.Parameter(
            displayName="Manual Search radius",
            name="search_radius",
            datatype="Double",
            parameterType="Optional",
            direction="Input",
            enabled=False
        )

        home_cutoff = arcpy.Parameter(
            displayName="Home cutoff percentage",
            name="home_cutoff",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )
        home_cutoff.filter.type = "Range"
        home_cutoff.filter.list = [0, 100]
        home_cutoff.value = 95.0

        core_cutoff = arcpy.Parameter(
            displayName="Core area cutoff percentage",
            name="core_cutoff",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )

        core_cutoff.filter.type = "Range"
        core_cutoff.filter.list = [0, 100]
        core_cutoff.value = 50

        out_cell_size = arcpy.Parameter(
            displayName="Output cell size",
            name="out_cell_size",
            datatype="analysis_cell_size",
            parameterType="Required",
            direction="Output"
        )

        area_unit = arcpy.Parameter(
            displayName="Area unit",
            name="area_unit",
            datatype="GPString",
            direction="Input",
            parameterType="Required"
        )
        area_unit.filter.type = 'ValueList'
        area_unit.filter.list = ['Square meters', 'Square kilometers', 'Hectares', 'Square feet (int)', 'Square miles (int)', 'Acres']
        area_unit.value = 'Square meters'

        out_folder = arcpy.Parameter(
            displayName="Output folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Output"
        )

        params = [observations, barrier_features, radius_method, search_radius, home_cutoff, core_cutoff, out_cell_size, area_unit, out_folder]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        if parameters[2].value == 'Manual':
            parameters[3].enabled = True
        else:
            parameters[3].enabled = False
        return

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

        # Check and set up parameters
        if params['barrier_features'].value:
            barrier_path = arcpy.Describe(params['barrier_features'].valueAsText).catalogPath
        else:
            barrier_path = None

        out_cell_size = params['out_cell_size'].valueAsText
        home_cutoff = params['home_cutoff'].valueAsText
        core_cutoff = params['core_cutoff'].valueAsText
        areal_unit = params['area_unit'].valueAsText

        if params['radius_method'].valueAsText == 'Reference bandwidth (Silverman 1986)':
            coord_cursor = arcpy.da.SearchCursor(obs_path, ["SHAPE@XY"])
            bandwidth = processor.calculate_reference_bandwidth(coord_cursor)
            arcpy.AddMessage('Reference bandwidth based on Silverman 1986: %s' % bandwidth)
        elif params['radius_method'].valueAsText == 'Least Squares Cross Validation (Seaman & Powell 1996)':
            coord_cursor = arcpy.da.SearchCursor(obs_path, ["SHAPE@XY"])
            bandwidth = processor.calculate_lscv_bandwidth(coord_cursor)
            arcpy.AddMessage('LSCV bandwidth based on on Seaman & Powell 1996: %s' % bandwidth)
        elif params['radius_method'].valueAsText == 'Manual':
            bandwidth = params['search_radius'].valueAsText
            arcpy.AddMessage('Bandwidth set manually: %s' % bandwidth)

        processor.compute_kde(arcpy, point_fc=obs_path, suffix=None, barrier_features=barrier_path, \
                                            cell_size=out_cell_size, search_radius=bandwidth, home_cutoff=home_cutoff, \
                                            core_cutoff=core_cutoff, areal_unit=areal_unit)
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

        radius_method = arcpy.Parameter(
                displayName='Bandwidth calculation method', 
                name='radius_method', 
                datatype='GPString', 
                parameterType='Optional', 
                direction="Input", 
                multiValue=False)

        radius_method.filter.type = "ValueList"
        radius_method.filter.list = ['Reference bandwidth (Silverman 1986)', 'Least Squares Cross Validation (Seaman & Powell 1996)', 'Manual']

        search_radius = arcpy.Parameter(
            displayName="Manual Search radius",
            name="search_radius",
            datatype="Double",
            parameterType="Optional",
            direction="Input",
            enabled=False
        )

        home_cutoff = arcpy.Parameter(
            displayName="Home cutoff percentage",
            name="home_cutoff",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )
        home_cutoff.filter.type = "Range"
        home_cutoff.filter.list = [0, 100]
        home_cutoff.value = 95.0

        core_cutoff = arcpy.Parameter(
            displayName="Core area cutoff percentage",
            name="core_cutoff",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )

        core_cutoff.filter.type = "Range"
        core_cutoff.filter.list = [0, 100]
        core_cutoff.value = 50

        out_cell_size = arcpy.Parameter(
            displayName="Output cell size",
            name="out_cell_size",
            datatype="analysis_cell_size",
            parameterType="Required",
            direction="Output"
        )

        area_unit = arcpy.Parameter(
            displayName="Area unit",
            name="area_unit",
            datatype="GPString",
            direction="Input",
            parameterType="Required"
        )
        area_unit.filter.type = 'ValueList'
        area_unit.filter.list = ['Square meters', 'Square kilometers', 'Hectares', 'Square feet (int)', 'Square miles (int)', 'Acres']
        area_unit.value = 'Square meters'

        out_folder = arcpy.Parameter(
            displayName="Output folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Output"
        )

        params = [observations, animal_id, barrier_features, radius_method, search_radius, home_cutoff, core_cutoff, out_cell_size, area_unit, out_folder]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""

        if parameters[3].value == 'Manual':
            parameters[4].enabled = True
        else:
            parameters[4].enabled = False
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
            if i % 5 == 0:
                arcpy.AddMessage("%s/%s individuals done." % (i, len(individuals)))
                #break

            tmp_points = r'tmp.shp'
                      
            where_txt = '"%s" = \'%s\'' % (param_values['animal_id'], animal)
           
            arcpy.analysis.Select(obs_path, tmp_points, where_clause=where_txt)
            tmp_points_path = arcpy.Describe(tmp_points).catalogPath


            # Check and set up parameters
            if params['barrier_features'].value:
                barrier_path = arcpy.Describe(params['barrier_features'].valueAsText).catalogPath
            else:
                barrier_path = None

            out_cell_size = params['out_cell_size'].valueAsText
            home_cutoff = params['home_cutoff'].valueAsText
            core_cutoff = params['core_cutoff'].valueAsText
            areal_unit = processor.area_unit_lookup(params['area_unit'].valueAsText)

            if params['radius_method'].valueAsText == 'Reference bandwidth (Silverman 1986)':
                coord_cursor = arcpy.da.SearchCursor(tmp_points_path, ["SHAPE@XY"])
                bandwidth = processor.calculate_reference_bandwidth(coord_cursor)
                arcpy.AddMessage('Reference bandwidth based on Silverman 1986: %s' % bandwidth)
            elif params['radius_method'].valueAsText == 'Least Squares Cross Validation (Seaman & Powell 1996)':
                coord_cursor = arcpy.da.SearchCursor(tmp_points_path, ["SHAPE@XY"])
                bandwidth = processor.calculate_lscv_bandwidth(coord_cursor)
                arcpy.AddMessage('LSCV bandwidth based on on Seaman & Powell 1996: %s' % bandwidth)
            elif params['radius_method'].valueAsText == 'Manual':
                bandwidth = params['search_radius'].valueAsText
                arcpy.AddMessage('Bandwidth set manually: %s' % bandwidth)

            processor.compute_kde(arcpy, point_fc=tmp_points_path, suffix=animal, barrier_features=barrier_path, \
                                            cell_size=out_cell_size, search_radius=bandwidth, home_cutoff=home_cutoff, \
                                            core_cutoff=core_cutoff, areal_unit=areal_unit, iteration=i)

            arcpy.management.Delete(tmp_points_path)

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
        observations = arcpy.Parameter(
            displayName="Input observations",
            name="observations",
            datatype="GPFeatureLayer",
            parameterType="Required",
            direction="Input"
        )

        observations.filter.list = ['Point', 'MultiPoint']

        home_cutoff = arcpy.Parameter(
            displayName="Home cutoff percentage",
            name="home_cutoff",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )
        home_cutoff.filter.type = "Range"
        home_cutoff.filter.list = [0, 100]
        home_cutoff.value = 95.0

        core_cutoff = arcpy.Parameter(
            displayName="Core area cutoff percentage",
            name="core_cutoff",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )

        core_cutoff.filter.type = "Range"
        core_cutoff.filter.list = [0, 100]
        core_cutoff.value = 50

        area_unit = arcpy.Parameter(
            displayName="Area unit",
            name="area_unit",
            datatype="GPString",
            direction="Input",
            parameterType="Required"
        )
        area_unit.filter.type = 'ValueList'
        area_unit.filter.list = ['Square meters', 'Square kilometers', 'Hectares', 'Square feet (int)', 'Square miles (int)', 'Acres']
        area_unit.value = 'Square meters'

        out_folder = arcpy.Parameter(
            displayName="Output folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Output"
        )

        params = [observations, home_cutoff, core_cutoff, area_unit, out_folder]
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

        areal_unit = processor.area_unit_lookup(params['area_unit'].valueAsText)

        processor.compute_mcp(arcpy, obs_path, suffix=None, home_cutoff=params['home_cutoff'].valueAsText, core_cutoff=params['core_cutoff'].valueAsText, areal_unit=areal_unit)

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class HomeRangeMCP_Batch(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Batch Home Range Estimation using MCP"
        self.description = "Batch Home Range Estimation using Minimum Convex Polygons"
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

        home_cutoff = arcpy.Parameter(
            displayName="Home cutoff percentage",
            name="home_cutoff",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )
        home_cutoff.filter.type = "Range"
        home_cutoff.filter.list = [0, 100]
        home_cutoff.value = 95.0

        core_cutoff = arcpy.Parameter(
            displayName="Core area cutoff percentage",
            name="core_cutoff",
            datatype="Double",
            parameterType="Required",
            direction="Input"
        )

        core_cutoff.filter.type = "Range"
        core_cutoff.filter.list = [0, 100]
        core_cutoff.value = 50

        out_folder = arcpy.Parameter(
            displayName="Output folder",
            name="out_folder",
            datatype="DEFolder",
            parameterType="Required",
            direction="Output"
        )

        area_unit = arcpy.Parameter(
            displayName="Area unit",
            name="area_unit",
            datatype="GPString",
            direction="Input",
            parameterType="Required"
        )
        area_unit.filter.type = 'ValueList'
        area_unit.filter.list = ['Square meters', 'Square kilometers', 'Hectares', 'Square feet (int)', 'Square miles (int)', 'Acres']
        area_unit.value = 'Square meters'

        params = [observations, animal_id, home_cutoff, core_cutoff, area_unit, out_folder]
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
            if i % 5 == 0:
                arcpy.AddMessage("%s/%s individuals done." % (i, len(individuals)))
                #break

            tmp_points = r'tmp.shp'
            out_home = r'home_range_mcp.shp'
            out_core = r'core_range_mcp.shp'

            where_txt = '"%s" = \'%s\'' % (param_values['animal_id'], animal)
           
            arcpy.analysis.Select(obs_path, tmp_points, where_clause=where_txt)
            tmp_points_path = arcpy.Describe(tmp_points).catalogPath

            areal_unit = processor.area_unit_lookup(params['area_unit'].valueAsText)

            processor.compute_mcp(arcpy, tmp_points_path, suffix=animal, home_cutoff=params['home_cutoff'].valueAsText, core_cutoff=params['core_cutoff'].valueAsText, areal_unit=areal_unit, iteration=i)

            arcpy.management.Delete(tmp_points_path)

        return

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return

class HomeRangeCalc(object):
    """Class that contains methods needed to calculate home ranges."""
    def area_unit_lookup(self, unit):
        units = {
                "Square meters": "SQUARE_METERS",
                "Square kilometers": "SQUARE_KILOMETERS",
                "Hectares": "HECTARES",
                "Square feet (int)": "SQUARE_FEET_INT",
                "Square miles (int)": "SQUARE_MILES_INT",
                "Acres": "ACRES_US"
        }
        return units[unit]

    def compute_kde(self, arcpy, point_fc, suffix, barrier_features, cell_size, search_radius, home_cutoff, core_cutoff, areal_unit, iteration):
        """Implements Kernel Density Estimation

        If the argument `suffix` is None, the method is called from the HomeRangeKDE function. If `suffix` is set,
        the HomeRangeKDE_Batch Tool is being used.

        Parameters
        ----------
        arcpy : module
            arcpy module that contains environmental settings and ArcGIS functionality to be called.
        point_fc : str
            Path to a point Feature Class that contains point observations. Type must be Point or MultiPoint
        suffix : str, optional
            Unique animal ID. Used when iterating through individuals stored in the same input file. Called
            from a loop that iterates over unique animal IDs.
        barrier_features: str, optional
            Path to a Feature Class with barriers.
        cell_size
            Spatial resolution of output KDE raster.
        search_radius : double
            Bandwidth of Kernel used in density estimation.
        home_cutoff : double
            Cutoff percentage value to determine home range. Values must be between 0 and 100. Default is 95%.
        core_cutoff : double
            Cutoff percentage value to determine core area. Values must be between 0 and 100. Default is 50%.
        """

        # Build KDE raster
        if barrier_features:
            raster = arcpy.sa.KernelDensity(point_fc, None, cell_size, search_radius, in_barriers=barrier_features)
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

        home_cutoff = int(num_records * float(home_cutoff)/100.0)
        core_cutoff = int(num_records * float(core_cutoff)/100.0)

        home_cutoff_value = value_list[home_cutoff - 1]
        core_cutoff_value = value_list[core_cutoff - 1]

        max_kernel_value = arcpy.management.GetRasterProperties(out_raster_name, 'MAXIMUM')

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
        arcpy.conversion.RasterToPolygon(home_raster, out_home_name + 'poly.shp', 'SIMPLIFY', 'Value')
        arcpy.conversion.RasterToPolygon(core_raster, out_core_name + 'poly.shp', 'SIMPLIFY', 'Value')

        arcpy.management.Dissolve(out_home_name + 'poly.shp', out_home_name + '.shp' )
        arcpy.management.Dissolve(out_core_name + 'poly.shp', out_core_name + '.shp' )

        # Add field for animal ID
        if suffix:
            arcpy.management.AddField(out_home_name, 'subject', 'TEXT')
            arcpy.management.AddField(out_core_name, 'subject', 'TEXT')

        # Add field for area calculation
        arcpy.management.AddField(out_home_name, 'area', 'DOUBLE')
        arcpy.management.AddField(out_core_name, 'area', 'DOUBLE')

        # Calculate fields
        if suffix:
            arcpy.management.CalculateField(out_home_name, 'subject', "'" + suffix + "'", 'PYTHON3')
            arcpy.management.CalculateField(out_core_name, 'subject', "'" + suffix + "'", 'PYTHON3')

        arcpy.management.CalculateGeometryAttributes(out_home_name, [['area', 'AREA']], area_unit=areal_unit)
        arcpy.management.CalculateGeometryAttributes(out_core_name, [['area', 'AREA']], area_unit=areal_unit)


        # Create final output files
        if iteration==1 and suffix:
            arcpy.management.Copy(out_home_name, 'home_range.shp')
            arcpy.management.Copy(out_core_name, 'core_range.shp')
            arcpy.management.Delete(out_home_name)
            arcpy.management.Delete(out_core_name)
        elif iteration > 1 and suffix:
            arcpy.management.Append(out_home_name, 'home_range.shp')
            arcpy.management.Append(out_core_name, 'core_range.shp')
            arcpy.management.Delete(out_home_name)
            arcpy.management.Delete(out_core_name)

        arcpy.management.Delete(out_raster_path)
        arcpy.management.Delete(tmp_raster_points)
        arcpy.management.Delete(out_home_name + '.tif')
        arcpy.management.Delete(out_core_name + '.tif')
        arcpy.management.Delete(out_home_name + 'poly.shp')
        arcpy.management.Delete(out_core_name + 'poly.shp')
        # Do not delete point feature class if not called from within a loop (i.e. suffix is None)
        if suffix:
            arcpy.management.Delete(point_fc)

        # Barrier lines
        # Central Feature
        # Calculate Eucledian distance
        # Extract values from distance raster
        # Order distances, discard points for home and core
        # Calculate MCP with points

    def compute_mcp(self, arcpy, point_fc, suffix, home_cutoff, core_cutoff, areal_unit, iteration):
        """Summary
        
        Parameters
        ----------
        arcpy : TYPE
            Description
        point_fc : TYPE
            Description
        suffix : TYPE
            Description
        home_cutoff : TYPE
            Description
        core_cutoff : TYPE
            Description
        """
        central = r'central_feature.shp'
        dist_raster = r'distance.tif'
        tmp_points = r'selection.shp'

        if suffix:
            mcp_home = r'mcp_home_%s.shp' % suffix
            mcp_core = r'mcp_core_%s.shp' % suffix
        else:
            mcp_home = r'mcp_home.shp'
            mcp_core = r'mcp_core.shp'
            out_raster_name = "kde.tif"

        arcpy.stats.CentralFeature(point_fc, central, 'EUCLIDEAN_DISTANCE')

        arcpy.analysis.Near(point_fc, central)

        distances = []

        with arcpy.da.SearchCursor(point_fc, 'NEAR_DIST') as cursor:
            for val in cursor:
                distances.append(val[0])

        distances.sort(reverse=False)
        num_records = len(distances)

        home_cutoff = int(num_records * float(home_cutoff)/100.0)
        core_cutoff = int(num_records * float(core_cutoff)/100.0)

        home_cutoff_value = distances[home_cutoff - 1]
        core_cutoff_value = distances[core_cutoff - 1]

        arcpy.AddMessage(home_cutoff_value)
        arcpy.AddMessage(core_cutoff_value)

        # Home
        where_clause_home = '"NEAR_DIST" < %s' % home_cutoff_value 
        arcpy.analysis.Select(point_fc, tmp_points, where_clause=where_clause_home)
        arcpy.management.MinimumBoundingGeometry(tmp_points, mcp_home, geometry_type='CONVEX_HULL')

        # Core
        where_clause_core = '"NEAR_DIST" < %s' % core_cutoff_value 
        arcpy.analysis.Select(point_fc, tmp_points, where_clause=where_clause_core)
        arcpy.management.MinimumBoundingGeometry(tmp_points, mcp_core, geometry_type='CONVEX_HULL')

        # Add field for animal ID
        if suffix:
            arcpy.management.AddField(mcp_home, 'subject', 'TEXT')
            arcpy.management.AddField(mcp_core, 'subject', 'TEXT')

        # Add field for area calculation
        arcpy.management.AddField(mcp_home, 'area', 'DOUBLE')
        arcpy.management.AddField(mcp_core, 'area', 'DOUBLE')

        # Calculate fields
        if suffix:
            arcpy.management.CalculateField(mcp_home, 'subject', "'" + suffix + "'", 'PYTHON3')
            arcpy.management.CalculateField(mcp_core, 'subject', "'" + suffix + "'", 'PYTHON3')

        arcpy.management.CalculateGeometryAttributes(mcp_home, [['area', 'AREA']], area_unit=areal_unit)
        arcpy.management.CalculateGeometryAttributes(mcp_core, [['area', 'AREA']], area_unit=areal_unit)

        # Create final output files
        if iteration==1 and suffix:
            arcpy.management.Copy(mcp_home, 'mcp_home.shp')
            arcpy.management.Copy(mcp_core, 'mcp_core.shp')
            arcpy.management.Delete(mcp_home)
            arcpy.management.Delete(mcp_core)
        elif iteration > 1 and suffix:
            arcpy.management.Append(mcp_home, 'mcp_home.shp')
            arcpy.management.Append(mcp_core, 'mcp_core.shp')
            arcpy.management.Delete(mcp_home)
            arcpy.management.Delete(mcp_core)

        arcpy.management.Delete(dist_raster)
        arcpy.management.Delete(tmp_points)
        arcpy.management.Delete(central)
        # Do not delete point feature class if not called from within a loop (i.e. suffix is None)
        if suffix:
            arcpy.management.Delete(point_fc)     


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
        from scipy.spatial import distance_matrix
        import numpy as np
        import math

        ITER = 100

        x_list = []
        y_list = []
        coords = []

        for feature in search_cursor:
            x, y = feature[0]
            coords.append(feature[0])
            x_list.append(x)
            y_list.append(y)

        n = len(x_list)

        search_cursor.reset()
        h_ref = self.calculate_reference_bandwidth(search_cursor)
        min_h = h_ref * 0.1
        max_h = h_ref * 2

        step = (max_h - min_h) / ITER

        value_range = [min_h + step * x for x in list(range(0,ITER+1))]

        f = distance_matrix(coords, coords)
        f = np.tril(f, k=-1) # keep lower triangle, no diagonal
        f_flat = f.flatten(order='F') # flatten matrix Fortran style (column major)
        f_flat = np.array([x for x in f_flat if x > 0])

        res = []
        for h in value_range:
            out = sum(numpy.exp(-f_flat * f_flat / (4 * h * h)) - 4 * numpy.exp(-f_flat * f_flat / (2 * h * h)))
            x = 1.0 / (math.pi * h * h * n) + (2 * out - 3 * n)/(math.pi * 4.0 * h * h * n * n)
            res.append(x)

        res = numpy.array(res)
        h_value = value_range[res.argmin()]
        return h_value


# def main():
#     tbx = Toolbox()

#     test_tool = HomeRangeKDE_Batch()
#     test_tool.execute(test_tool.getParameterInfo(), None)

# if __name__ == '__main__':
#     main()