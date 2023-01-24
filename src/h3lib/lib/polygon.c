/*
 * Copyright 2018-2021 Uber Technologies, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
/** @file polygon.c
 * @brief Polygon (GeoLoop) algorithms
 */

#include "polygon.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>

#include "bbox.h"
#include "constants.h"
#include "h3api.h"
#include "latLng.h"
#include "linkedGeo.h"

// Define macros used in polygon algos for GeoLoop
#define TYPE GeoLoop
#define INIT_ITERATION INIT_ITERATION_GEOFENCE
#define ITERATE ITERATE_GEOFENCE
#define IS_EMPTY IS_EMPTY_GEOFENCE

#include "polygonAlgos.h"

#undef TYPE
#undef INIT_ITERATION
#undef ITERATE
#undef IS_EMPTY

/**
 * Create a bounding box from a GeoPolygon
 * @param polygon Input GeoPolygon
 * @param bboxes  Output bboxes, one for the outer loop and one for each hole
 */
void bboxesFromGeoPolygon(const GeoPolygon *polygon, BBox *bboxes) {
    bboxFromGeoLoop(&polygon->geoloop, &bboxes[0]);
    for (int i = 0; i < polygon->numHoles; i++) {
        bboxFromGeoLoop(&polygon->holes[i], &bboxes[i + 1]);
    }
}

/**
 * pointInsidePolygon takes a given GeoPolygon data structure and
 * checks if it contains a given geo coordinate.
 *
 * @param geoPolygon The geoloop and holes defining the relevant area
 * @param bboxes     The bboxes for the main geoloop and each of its holes
 * @param coord      The coordinate to check
 * @return           Whether the point is contained
 */
bool pointInsidePolygon(const GeoPolygon *geoPolygon, const BBox *bboxes,
                        const LatLng *coord) {
    // Start with contains state of primary geoloop
    bool contains =
        pointInsideGeoLoop(&(geoPolygon->geoloop), &bboxes[0], coord);

    // If the point is contained in the primary geoloop, but there are holes in
    // the geoloop iterate through all holes and return false if the point is
    // contained in any hole
    if (contains && geoPolygon->numHoles > 0) {
        for (int i = 0; i < geoPolygon->numHoles; i++) {
            if (pointInsideGeoLoop(&(geoPolygon->holes[i]), &bboxes[i + 1],
                                   coord)) {
                return false;
            }
        }
    }

    return contains;
}

bool geoLoopInsidePolygon(const GeoPolygon *geoPolygon, const BBox *bboxes,
                          const GeoLoop *loop) {
    // Check for loop points that are outside of the polygon
    for (int i = 0; i < loop->numVerts; i++) {
        if (!pointInsidePolygon(geoPolygon, bboxes, &loop->verts[i])) {
            return false;
        }
    }
    if (loop->numVerts < 2) {
        return loop->numVerts != 0;
    }

    // Check for hole points inside the loop
    if (loop->numVerts > 2 && geoPolygon->numHoles > 0) {
        BBox loopBbox;
        bboxFromGeoLoop(loop, &loopBbox);

        for (int i = 0; i < geoPolygon->numHoles; i++) {
            const GeoLoop *hole = &(geoPolygon->holes[i]);
            for (int j = 0; j < hole->numVerts; j++) {
                if (pointInsideGeoLoop(loop, &loopBbox, &(hole->verts[j]))) {
                    return false;
                }
            }
        }
    }

    // Check for loop segments intersecting the outer loop or the holes
    for (int i = 0; i < loop->numVerts; i++) {
        const LatLng *p1 = &loop->verts[i];
        const LatLng *p2 = &loop->verts[(i + 1) % loop->numVerts];

        if (segmentIntersectsGeoLoop(&geoPolygon->geoloop, &bboxes[0], p1,
                                     p2)) {
            return false;
        }

        for (int j = 0; j < geoPolygon->numHoles; j++) {
            if (segmentIntersectsGeoLoop(&(geoPolygon->holes[j]),
                                         &bboxes[j + 1], p1, p2)) {
                return false;
            }
        }
    }

    return true;
}

bool geoLoopIntersectsPolygon(const GeoPolygon *geoPolygon, const BBox *bboxes,
                              const GeoLoop *loop) {
    // For each point in the loop check if it is contained in the outer loop
    // and in one of the holes.
    // Return true if
    // - a point is inside the outer loop and not in a hole;
    // - if two points are in different holes (or one is in a hole, and the
    //   other is outside the outer loop)
    int holeIndex = -1;
    for (int i = 0; i < loop->numVerts; i++) {
        const LatLng *coord = &loop->verts[i];

        bool contains =
            pointInsideGeoLoop(&geoPolygon->geoloop, &bboxes[0], coord);

        if (contains && geoPolygon->numHoles > 0) {
            for (int j = 0; j < geoPolygon->numHoles; j++) {
                if (pointInsideGeoLoop(&(geoPolygon->holes[j]), &bboxes[j + 1],
                                       coord)) {
                    if (j > 0 && holeIndex != j) {
                        // Previous point is either outside the outer loop or in
                        // a different hole
                        return true;
                    }
                    holeIndex = j;
                    contains = false;
                    break;
                }
            }
        }
        if (contains) {
            return true;
        }
    }

    if (loop->numVerts > 1) {
        // All points are either outside the outer loop or inside a single hole.
        // Check if segments intersect corresponding loop.

        const GeoLoop *geoloop;
        const BBox *bbox;
        if (holeIndex == -1) {
            geoloop = &geoPolygon->geoloop;
            bbox = &(bboxes[0]);
        } else {
            geoloop = &(geoPolygon->holes[holeIndex]);
            bbox = &(bboxes[holeIndex]);
        }

        for (int i = 0; i < loop->numVerts; i++) {
            const LatLng *p0 = &loop->verts[i];
            const LatLng *p1 = &loop->verts[(i + 1) % loop->numVerts];
            if (segmentIntersectsGeoLoop(geoloop, bbox, p0, p1)) {
                return true;
            }
        }
    }

    return false;
}
