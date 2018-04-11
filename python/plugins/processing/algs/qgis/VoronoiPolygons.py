# -*- coding: utf-8 -*-

"""
***************************************************************************
    VoronoiPolygons.py
    ---------------------
    Date                 : August 2012
    Copyright            : (C) 2012 by Victor Olaya
    Email                : volayaf at gmail dot com
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""
# For comparing equality in lat/long - 10^6 meters are < 1 degree
# so this gives mm resolution at lat=0. For finer resolution tasks
# a localised projection will likely be used, so this should be OK.
FLT_EPSILON = 10**-9

__author__ = 'Victor Olaya'
__date__ = 'August 2012'
__copyright__ = '(C) 2012, Victor Olaya'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = 'f378a23ed8e56d225ab5734766fc26fa5d832369'

import os

from qgis.PyQt.QtGui import QIcon

from qgis.core import (QgsFeatureRequest,
                       QgsFeatureSink,
                       QgsFeature,
                       QgsGeometry,
                       QgsPointXY,
                       QgsWkbTypes,
                       QgsProcessing,
                       QgsProcessingException,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterNumber)

from processing.algs.qgis.QgisAlgorithm import QgisAlgorithm

from . import voronoi

pluginPath = os.path.split(os.path.split(os.path.dirname(__file__))[0])[0]


class VoronoiPolygons(QgisAlgorithm):

    INPUT = 'INPUT'
    BUFFER = 'BUFFER'
    OUTPUT = 'OUTPUT'

    def icon(self):
        return QIcon(os.path.join(pluginPath, 'images', 'ftools', 'voronoi.png'))

    def group(self):
        return self.tr('Vector geometry')

    def groupId(self):
        return 'vectorgeometry'

    def __init__(self):
        super().__init__()

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(self.INPUT, self.tr('Input layer'), [QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterNumber(self.BUFFER, self.tr('Buffer region'), type=QgsProcessingParameterNumber.Double,
                                                       minValue=0.0, maxValue=9999999999, defaultValue=0.0))

        self.addParameter(QgsProcessingParameterFeatureSink(self.OUTPUT, self.tr('Voronoi polygons'), type=QgsProcessing.TypeVectorPolygon))

    def name(self):
        return 'voronoipolygons'

    def displayName(self):
        return self.tr('Voronoi polygons')

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        buf = self.parameterAsDouble(parameters, self.BUFFER, context)
        (sink, dest_id) = self.parameterAsSink(parameters, self.OUTPUT, context,
                                               source.fields(), QgsWkbTypes.Polygon, source.sourceCrs())

        outFeat = QgsFeature()
        extent = source.sourceExtent()
        extraX = extent.height() * (buf / 100.0)
        extraY = extent.width() * (buf / 100.0)
        extent.setXMinimum(extent.xMinimum() - extraX)
        extent.setYMinimum(extent.yMinimum() - extraY)
        extent.setXMaximum(extent.xMaximum() + extraX)
        extent.setYMaximum(extent.yMaximum() + extraY)

        height = extent.height()
        width = extent.width()
        c = voronoi.Context()
        pts = []
        ptDict = {}
        ptNdx = -1

        features = source.getFeatures()
        total = 100.0 / source.featureCount() if source.featureCount() else 0
        for current, inFeat in enumerate(features):
            if feedback.isCanceled():
                break
            geom = inFeat.geometry()
            point = geom.asPoint()
            x = point.x() - extent.xMinimum()
            y = point.y() - extent.yMinimum()
            pts.append((x, y))
            ptNdx += 1
            ptDict[ptNdx] = inFeat.id()
            feedback.setProgress(int(current * total))

        if len(pts) < 3:
            raise QgsProcessingException(
                self.tr('Input file should contain at least 3 points. Choose '
                        'another file and try again.'))

        uniqueSet = set(item for item in pts)
        ids = [pts.index(item) for item in uniqueSet]
        sl = voronoi.SiteList([voronoi.Site(i[0], i[1], sitenum=j) for (j,
                                                                        i) in enumerate(uniqueSet)])
        voronoi.voronoi(sl, c)
        inFeat = QgsFeature()

        current = 0
        if len(c.polygons) == 0:
            raise QgsProcessingException(
                self.tr('There were no polygons created.'))

        total = 100.0 / len(c.polygons)

        for (site, edges) in list(c.polygons.items()):
            if feedback.isCanceled():
                break

            request = QgsFeatureRequest().setFilterFid(ptDict[ids[site]])
            inFeat = next(source.getFeatures(request))
            lines = self.clip_voronoi(edges, c, width, height, extent)

            geom = QgsGeometry.fromMultiPointXY(lines)
            geom = QgsGeometry(geom.convexHull())
            outFeat.setGeometry(geom)
            outFeat.setAttributes(inFeat.attributes())
            sink.addFeature(outFeat, QgsFeatureSink.FastInsert)

            current += 1
            feedback.setProgress(int(current * total))

        return {self.OUTPUT: dest_id}

    def clip_voronoi(self, edges, c, width, height, extent):
        """Clip voronoi function based on code written for Inkscape.
        Copyright (C) 2010 Alvin Penner, penner@vaxxine.com
        """
        def float_eq( a, b ):
            # the initial == test is for efficiency
            return a == b or abs( a - b ) <= FLT_EPSILON

        def clip_line(x1, y1, x2, y2, w, h,):
            # where both x's are left of the clipping region no line needed
            if x1 < 0 and x2 < 0:
                return [0, 0, 0, 0]
            # where both x's are right of the clipping region no line needed
            if x1 > w and x2 > w:
                return [0, 0, 0, 0]
            # where both y's are below the clipping region no line needed
            if y1 < 0 and y2 < 0:
                return [0, 0, 0, 0]
            # where both y's are above the clipping region no line needed
            if y1 > h and y2 > h:
                return [0, 0, 0, 0]
            # where end points are too close to tell apart no line needed
            if float_eq( x1, x2 ) and float_eq( y1, y2 ):
                return [0, 0, 0, 0]
            # y = mx + b: calculate gradient but make sure is not a vertical line
            if not float_eq( x1, x2 ):
                m = (y2 - y1) / (x2 - x1)
            else:
                m = float("+inf") # vertical lines have infinite gradient
            # y = mx + b: calculate intercept, if horizontal m==0 so use first Y value
            if ( not float_eq( y1, y2 ) and not float_eq( x1, x2 ) ) or float_eq( y1, y2):
                b = y1 - m*x1 # b = y1 when y1 == y2
            else: # y1 != y2 and x1 == x2
                b = float("nan") # vertical lines have no intercept
            # for simplicity later define if the line is vertical
            is_vertical = ( m == float("+inf") and b == float("nan") )
            # where either of the x1,y1 or x2,y2 points are outside of the clipping region
            # need to move them to be on the border of the region; this is done by recalculating
            # the x or y value based on the equation of the line holding the other axis as the
            # value of the clipping rectangle border that is violated.
            if x1 < 0: # first point is left of region
                # impossible to have vertical line here, could be horizontal
                # but that doesn't matter as m = 0, b = y1 in that case
                x1 = 0
                y1 = m * x1 + b
            if x2 < 0: # second point is left of region
                # impossible to have vertical line here, could be horizontal
                # but that doesn't matter as m = 0, b = y1 in that case
                x2 = 0
                y2 = m * x2 + b
            if x1 > w: # first point is right of region
                # impossible to have vertical line here, could be horizontal
                # but that doesn't matter as m = 0, b = y1 in that case
                x1 = w
                y1 = m * x1 + b
            if x2 > w: # second point is right of region
                # impossible to have vertical line here, could be horizontal
                # but that doesn't matter as m = 0, b = y1 in that case
                x2 = w
                y2 = m * x2 + b
            if y1 < 0: # first point is below the region
                # impossible to have horizontal line here, could be vertical though
                # if is vertical then x value does not need to change
                y1 = 0
                if not is_vertical: # if veritical then X doesn't need to change
                    x1 = (y1 - b) / m
            if y2 < 0: # second point is below the region
                # impossible to have horizontal line here, could be vertical though
                # if is vertical then x value does not need to change
                y2 = 0
                if not is_vertical: # if vertical then X doesn't need to change
                    x2 = (y2 - b) / m
            if y1 > h: # first point is above the region
                # impossible to have horizontal line here, could be vertical though
                # if is vertical then x value does not need to change
                y1 = h

                if not is_vertical: # if vertical then X doesn't need to change
                    x1 = (y1 - b) / m
            if y2 > h: # second point is above the region
                # impossible to have horizontal line here, could be vertical though
                # if is vertical then x value does not need to change
                y2 = h

                if not is_vertical: # if vertical then X doesn't need to change
                    x2 = (y2 - b) / m
            # recheck point approximate equivalance because both points may have been
            # outside the region and thus were both clipped to a point on the boundary
            if float_eq( x1, x2 ) and float_eq( y1, y2 ):
                return [ 0, 0, 0, 0 ]
            # return the clipped region
            return [x1, y1, x2, y2]
        lines = []
        hasXMin = False
        hasYMin = False
        hasXMax = False
        hasYMax = False
        for edge in edges:
            if edge[1] >= 0 and edge[2] >= 0:
                # Two vertices
                [x1, y1, x2, y2] = clip_line(
                    c.vertices[edge[1]][0],
                    c.vertices[edge[1]][1],
                    c.vertices[edge[2]][0],
                    c.vertices[edge[2]][1],
                    width,
                    height,
                )
            elif edge[1] >= 0:
                # Only one vertex
                if c.lines[edge[0]][1] == 0:
                    # Vertical line
                    xtemp = c.vertices[edge[1]][0]
                    if c.vertices[edge[1]][1] > height / 2:
                        ytemp = height
                    else:
                        ytemp = 0.0
                else:
                    xtemp = width
                    ytemp = (c.lines[edge[0]][2] - xtemp *
                             c.lines[edge[0]][0]) / c.lines[edge[0]][1]
                [x1, y1, x2, y2] = clip_line(
                    c.vertices[edge[1]][0],
                    c.vertices[edge[1]][1],
                    xtemp,
                    ytemp,
                    width,
                    height,
                )
            elif edge[2] >= 0:
                # Only one vertex
                if c.lines[edge[0]][1] == 0:
                    # Vertical line
                    xtemp = c.vertices[edge[2]][0]
                    if c.vertices[edge[2]][1] > height / 2:
                        ytemp = height
                    else:
                        ytemp = 0.0
                else:
                    xtemp = 0.0
                    ytemp = (c.lines[edge[0]][2] - xtemp *
                             c.lines[edge[0]][0]) / c.lines[edge[0]][1]
                [x1, y1, x2, y2] = clip_line(
                    xtemp,
                    ytemp,
                    c.vertices[edge[2]][0],
                    c.vertices[edge[2]][1],
                    width,
                    height,
                )
            if x1 or x2 or y1 or y2:
                lines.append(QgsPointXY(x1 + extent.xMinimum(),
                                        y1 + extent.yMinimum()))
                lines.append(QgsPointXY(x2 + extent.xMinimum(),
                                        y2 + extent.yMinimum()))
                if 0.0 in (x1, x2):
                    hasXMin = True
                if 0.0 in (y1, y2):
                    hasYMin = True
                if height in (y1, y2):
                    hasYMax = True
                if width in (x1, x2):
                    hasXMax = True
        if hasXMin:
            if hasYMax:
                lines.append(QgsPointXY(extent.xMinimum(),
                                        height + extent.yMinimum()))
            if hasYMin:
                lines.append(QgsPointXY(extent.xMinimum(),
                                        extent.yMinimum()))
        if hasXMax:
            if hasYMax:
                lines.append(QgsPointXY(width + extent.xMinimum(),
                                        height + extent.yMinimum()))
            if hasYMin:
                lines.append(QgsPointXY(width + extent.xMinimum(),
                                        extent.yMinimum()))
        return lines
