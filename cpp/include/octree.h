/*
 * GPisMap - Online Continuous Mapping using Gaussian Process Implicit Surfaces
 * https://github.com/leebhoram/GPisMap
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License v3 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of any FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU General Public License v3 for more details.
 *
 * You should have received a copy of the GNU General Public License v3
 * along with this program; if not, you can access it online at
 * http://www.gnu.org/licenses/gpl-3.0.html.
 *
 * Authors: Bhoram Lee <bhoram.lee@gmail.com>
 */

#ifndef __OCTREE_H_
#define __OCTREE_H_

#include <vector>
#include <unordered_set>
#include <memory>
#include <iostream>
#include <cstdint>
#include "strct.h" // FLOAT is defined in this file
#include "OnGPIS.h"

class AABB3
{
    Point3<FLOAT> center;
    FLOAT halfLength;
    FLOAT halfLengthSq;
    FLOAT xmin;
    FLOAT xmax;
    FLOAT ymin;
    FLOAT ymax;
    FLOAT zmin;
    FLOAT zmax;

    Point3<FLOAT> ptNWF;
    Point3<FLOAT> ptNEF;
    Point3<FLOAT> ptSWF;
    Point3<FLOAT> ptSEF;
    Point3<FLOAT> ptNWB;
    Point3<FLOAT> ptNEB;
    Point3<FLOAT> ptSWB;
    Point3<FLOAT> ptSEB;
public:
    AABB3(){
        halfLength = 0.0;
        halfLengthSq = 0.0;
        xmin = 0.0;
        xmax = 0.0;
        ymin = 0.0;
        ymax = 0.0;
        zmin = 0.0;
        zmax = 0.0;
    }
    AABB3(Point3<FLOAT> _center, FLOAT _halfLength) {
        center = _center;
        halfLength = _halfLength;
        halfLengthSq = halfLength*halfLength;
        xmin = center.x - halfLength;
        xmax = center.x + halfLength;
        ymin = center.y - halfLength;
        ymax = center.y + halfLength;
        zmin = center.z - halfLength;
        zmax = center.z + halfLength;
        ptNWF = Point3<FLOAT>(xmin,ymax,zmax);
        ptNEF = Point3<FLOAT>(xmax,ymax,zmax);
        ptSWF = Point3<FLOAT>(xmin,ymin,zmax);
        ptSEF = Point3<FLOAT>(xmax,ymin,zmax);
        ptNWB = Point3<FLOAT>(xmin,ymax,zmin);
        ptNEB = Point3<FLOAT>(xmax,ymax,zmin);
        ptSWB = Point3<FLOAT>(xmin,ymin,zmin);
        ptSEB = Point3<FLOAT>(xmax,ymin,zmin);
    }
    AABB3(FLOAT x, FLOAT y, FLOAT z, FLOAT _halfLength) {
        center = Point3<FLOAT>(x,y,z);
        halfLength = _halfLength;
        halfLengthSq = halfLength*halfLength;
        xmin = center.x - halfLength;
        xmax = center.x + halfLength;
        ymin = center.y - halfLength;
        ymax = center.y + halfLength;
        zmin = center.z - halfLength;
        zmax = center.z + halfLength;
        ptNWF = Point3<FLOAT>(xmin,ymax,zmax);
        ptNEF = Point3<FLOAT>(xmax,ymax,zmax);
        ptSWF = Point3<FLOAT>(xmin,ymin,zmax);
        ptSEF = Point3<FLOAT>(xmax,ymin,zmax);
        ptNWB = Point3<FLOAT>(xmin,ymax,zmin);
        ptNEB = Point3<FLOAT>(xmax,ymax,zmin);
        ptSWB = Point3<FLOAT>(xmin,ymin,zmin);
        ptSEB = Point3<FLOAT>(xmax,ymin,zmin);
    }

    const Point3<FLOAT> getCenter(){return center;}
    FLOAT getHalfLength(){return halfLength;}
    FLOAT getHalfLengthSq(){return halfLengthSq;}
    FLOAT getXMinbound(){return xmin;}
    FLOAT getXMaxbound(){return xmax;}
    FLOAT getYMinbound(){return ymin;}
    FLOAT getYMaxbound(){return ymax;}
    FLOAT getZMinbound(){return zmin;}
    FLOAT getZMaxbound(){return zmax;}
    const Point3<FLOAT>& getNWF(){return ptNWF;}
    const Point3<FLOAT>& getNEF(){return ptNEF;}
    const Point3<FLOAT>& getSWF(){return ptSWF;}
    const Point3<FLOAT>& getSEF(){return ptSEF;}
    const Point3<FLOAT>& getNWB(){return ptNWB;}
    const Point3<FLOAT>& getNEB(){return ptNEB;}
    const Point3<FLOAT>& getSWB(){return ptSWB;}
    const Point3<FLOAT>& getSEB(){return ptSEB;}

    bool containsPoint(Point3<FLOAT> pt) {
        return ((pt.x > xmin) &&
            (pt.x < xmax) &&
            (pt.y > ymin) &&
            (pt.y < ymax) &&
            (pt.z > zmin) &&
            (pt.z < zmax)) ;
    }

    bool intersectsAABB(AABB3 aabb) {
        return !((aabb.getXMaxbound() < xmin) ||
                (aabb.getXMinbound() > xmax) ||
                (aabb.getYMaxbound() < ymin) ||
                (aabb.getYMinbound() > ymax) ||
                (aabb.getZMaxbound() < zmin) ||
                (aabb.getZMinbound() > zmax) );
    }
};

class OcTree
{
    // Arbitrary constant to indicate how many elements can be stored in this quad tree node
    const int CHILD_TYPE_NWF = 1;
    const int CHILD_TYPE_NEF = 2;
    const int CHILD_TYPE_SWF = 3;
    const int CHILD_TYPE_SEF = 4;
    const int CHILD_TYPE_NWB = 5;
    const int CHILD_TYPE_NEB = 6;
    const int CHILD_TYPE_SWB = 7;
    const int CHILD_TYPE_SEB = 8;

    // Axis-aligned bounding box stored as a center with half-dimensions
    // to represent the boundaries of this quad tree
    AABB3 boundary;

    static tree_param param;    // see strct.h for definition

    // Points in this quad tree node
    std::shared_ptr<Node3> node;
    std::shared_ptr<Node3> closestChildNode; // representative point (clasest to the center)
    std::shared_ptr<OnGPIS> gp;

    bool leaf;
    bool maxDepthReached;
    bool rootLimitReached;

    int32_t numNodes;

    // Children
    OcTree* northWestFront;
    OcTree* northEastFront;
    OcTree* southWestFront;
    OcTree* southEastFront;
    OcTree* northWestBack;
    OcTree* northEastBack;
    OcTree* southWestBack;
    OcTree* southEastBack;

    OcTree* par;

    void Subdivide(); // create four children that fully divide this quad into four quads of equal area
    void SubdivideExcept(int childType);
    void deleteChildren();
    bool InsertToParent(std::shared_ptr<Node3> n);

    OcTree(AABB3 _boundary, OcTree* const p =0 );
    OcTree(AABB3 _boundary, OcTree* const ch, int child_type);

protected:
    void setParent(OcTree* const p){par = p;}
    OcTree* const getParent(){return par;}
    bool IsLeaf(){return leaf;} // leaf is true if chlidren are initialized
    bool IsEmpty(){return (node == nullptr);} // empty if the data node is null
    bool IsEmptyLeaf(){
        return (leaf & (node == nullptr) ); // true if no data node no child
    }
public:
    // Methods
    OcTree():northWestFront(0),
            northEastFront(0),
            southWestFront(0),
            southEastFront(0),
            northWestBack(0),
            northEastBack(0),
            southWestBack(0),
            southEastBack(0),
            par(0),
            maxDepthReached(false),
            rootLimitReached(false),
            leaf(true),
            numNodes(0),
            node(nullptr),
            closestChildNode(nullptr),
            gp(nullptr){}

    OcTree(Point3<FLOAT> c);

    ~OcTree(){
        deleteChildren();
    }

    bool IsRoot(){
        if (par)
            return false;
        else
            return true;
    }

    OcTree* const getRoot();

    bool Insert(std::shared_ptr<Node3> n);
    bool Insert(std::shared_ptr<Node3> n, std::unordered_set<OcTree*>& quads);
    bool IsNotNew(std::shared_ptr<Node3> n);
    bool Update(std::shared_ptr<Node3> n);
    bool Update(std::shared_ptr<Node3> n, std::unordered_set<OcTree*>& quads);
    bool Remove(std::shared_ptr<Node3> n, std::unordered_set<OcTree*>& quads);

    void Update(std::shared_ptr<OnGPIS> _gp);
    std::shared_ptr<OnGPIS> const getGP(){return gp;}

    // Let's implement
    bool Remove(std::shared_ptr<Node3> n);
    void QueryRange(AABB3 range, std::vector<std::shared_ptr<Node3> >& nodes);
    void QueryNonEmptyLevelC(AABB3 range, std::vector<OcTree*>& quads);
    void QueryNonEmptyLevelC(AABB3 range, std::vector<OcTree*>& quads, std::vector<FLOAT>& sqdst);
    void QueryNonEmptyLevelC(AABB3 range, std::vector<OcTree*>& quads, std::vector<std::vector<std::shared_ptr<Node3> > >& nodes);


    int32_t getNodeCount(){return numNodes;}
    Point3<FLOAT> getCenter(){return boundary.getCenter();}
    FLOAT getHalfLength(){return boundary.getHalfLength();}
    FLOAT getXMaxbound(){return boundary.getXMaxbound();}
    FLOAT getYMinbound(){return boundary.getYMinbound();}
    FLOAT getYMaxbound(){return boundary.getYMaxbound();}
    FLOAT getZMinbound(){return boundary.getZMinbound();}
    FLOAT getZMaxbound(){return boundary.getZMaxbound();}
    Point3<FLOAT> getNWF(){return boundary.getNWF();}
    Point3<FLOAT> getNEF(){return boundary.getNEF();}
    Point3<FLOAT> getSWF(){return boundary.getSWF();}
    Point3<FLOAT> getSEF(){return boundary.getSEF();}
    Point3<FLOAT> getNWB(){return boundary.getNWB();}
    Point3<FLOAT> getNEB(){return boundary.getNEB();}
    Point3<FLOAT> getSWB(){return boundary.getSWB();}
    Point3<FLOAT> getSEB(){return boundary.getSEB();}

    void getAllChildrenNonEmptyNodes(std::vector<std::shared_ptr<Node3> >& nodes);

    void updateCount();
    void printNodes();
    void printBoundary();

};



#endif

