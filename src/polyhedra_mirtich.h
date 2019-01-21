//
//  polyhedra_mirtich.h
//  Batting3D
//
//  Created by Matthew Gardner on 5/5/18.
//

#ifndef polyhedra_mirtich_h
#define polyhedra_mirtich_h

#include "core.h"

/*
 ============================================================================
 constants
 ============================================================================
 */

#define MAX_VERTS 100     /* maximum number of polyhedral vertices */
#define MAX_FACES 100     /* maximum number of polyhedral faces */
#define MAX_POLYGON_SZ 10 /* maximum number of verts per polygonal face */

namespace mirtich
{
    /*
     ============================================================================
     data structures
     ============================================================================
     */

    typedef struct {
        int numVerts;
        double norm[3];
        double w;
        int verts[MAX_POLYGON_SZ];
        struct polyhedron *poly;
    } FACE;

    typedef struct polyhedron {
        int numVerts, numFaces;
        double verts[MAX_VERTS][3];
        FACE faces[MAX_FACES];
    } POLYHEDRON;

    /*
     ============================================================================
     data structures
     ============================================================================
     */

    void readPolyhedron(char *name, POLYHEDRON *p);
    void compProjectionIntegrals(FACE *f);
    void compFaceIntegrals(FACE *f);
    void compVolumeIntegrals(POLYHEDRON *p);
    void getMassProperties(double m, double vol, mg::Matrix3 &I, mg::Vec3 &c);
}

#endif /* polyhedra_mirtich_h */
