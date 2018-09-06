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
#ifndef __QUADTREE_H_
#define __QUADTREE_H_

#include <vector>
#include <unordered_set>
#include <memory>
#include <iostream>
#include <cstdint>
#include "strct.h" // FLOAT is defined in this file
#include "OnGPIS.h"

class AABB
{
    Point<FLOAT> center;
    FLOAT halfLength;
    FLOAT halfLengthSq;
    FLOAT xmin;
    FLOAT xmax;
    FLOAT ymin;
    FLOAT ymax;

    Point<FLOAT> ptNW;
    Point<FLOAT> ptNE;
    Point<FLOAT> ptSW;
    Point<FLOAT> ptSE;
public:
    AABB(){
        halfLength = 0.0;
        halfLengthSq = 0.0;
        xmin = 0.0;
        xmax = 0.0;
        ymin = 0.0;
        ymax = 0.0;
    }
    AABB(Point<FLOAT> _center, FLOAT _halfLength) {
        center = _center;
        halfLength = _halfLength;
        halfLengthSq = halfLength*halfLength;
        xmin = center.x - halfLength;
        xmax = center.x + halfLength;
        ymin = center.y - halfLength;
        ymax = center.y + halfLength;
        ptNW = Point<FLOAT>(xmin,ymax);
        ptNE = Point<FLOAT>(xmax,ymax);
        ptSW = Point<FLOAT>(xmin,ymin);
        ptSE = Point<FLOAT>(xmax,ymin);
    }
    AABB(FLOAT x, FLOAT y, FLOAT _halfLength) {
        center = Point<FLOAT>(x,y);
        halfLength = _halfLength;
        halfLengthSq = halfLength*halfLength;
        xmin = center.x - halfLength;
        xmax = center.x + halfLength;
        ymin = center.y - halfLength;
        ymax = center.y + halfLength;
        ptNW = Point<FLOAT>(xmin,ymax);
        ptNE = Point<FLOAT>(xmax,ymax);
        ptSW = Point<FLOAT>(xmin,ymin);
        ptSE = Point<FLOAT>(xmax,ymin);
    }

    const Point<FLOAT> getCenter(){return center;}
    FLOAT getHalfLength(){return halfLength;}
    FLOAT getHalfLengthSq(){return halfLengthSq;}
    FLOAT getXMinbound(){return xmin;}
    FLOAT getXMaxbound(){return xmax;}
    FLOAT getYMinbound(){return ymin;}
    FLOAT getYMaxbound(){return ymax;}
    const Point<FLOAT>& getNW(){return ptNW;}
    const Point<FLOAT>& getNE(){return ptNE;}
    const Point<FLOAT>& getSW(){return ptSW;}
    const Point<FLOAT>& getSE(){return ptSE;}

    bool containsPoint(Point<FLOAT> pt) {
        return ((pt.x > xmin) &&
            (pt.x < xmax) &&
            (pt.y > ymin) &&
            (pt.y < ymax)) ;
    }

    bool intersectsAABB(AABB aabb) {
        return !((aabb.getXMaxbound() < xmin) ||
            	(aabb.getXMinbound() > xmax) ||
                (aabb.getYMaxbound() < ymin) ||
            	(aabb.getYMinbound() > ymax) );
    }
};

class QuadTree
{
    // Arbitrary constant to indicate how many elements can be stored in this quad tree node
    //const int QT_NODE_CAPACITY = 1;
    const int CHILD_TYPE_NW = 1;
    const int CHILD_TYPE_NE = 2;
    const int CHILD_TYPE_SW = 3;
    const int CHILD_TYPE_SE = 4;

    // Axis-aligned bounding box stored as a center with half-dimensions
    // to represent the boundaries of this quad tree
    AABB boundary;

    static tree_param param;    // see strct.h for definition

    // Points in this quad tree node
    std::shared_ptr<Node> node;
    //std::shared_ptr<Node> closestChildNode; // representative point (clasest to the center)
    std::shared_ptr<OnGPIS> gp;

    bool leaf;
    bool maxDepthReached;
    bool rootLimitReached;

    int32_t numNodes;

    // Children
    QuadTree* northWest;
    QuadTree* northEast;
    QuadTree* southWest;
    QuadTree* southEast;

    QuadTree* par;

    QuadTree(AABB _boundary, QuadTree* const p =0 );
    QuadTree(AABB _boundary, QuadTree* const ch, int child_type);

    void Subdivide(); // create four children that fully divide this quad into four quads of equal area
    void SubdivideExcept(int childType);
    void deleteChildren();
    bool InsertToParent(std::shared_ptr<Node> n);
    void updateCount();
    void setParent(QuadTree* const p){par = p;}
protected:
    QuadTree* const getParent(){return par;}
    bool IsLeaf(){return leaf;} // leaf is true if chlidren are initialized
    bool IsEmpty(){return (node == nullptr);} // empty if the data node is null
    bool IsEmptyLeaf(){
        return (leaf & (node == nullptr) ); // true if no data node no child
    }
public:
    // Methods
    QuadTree():northWest(0),
            northEast(0),
            southWest(0),
            southEast(0),
            par(0),
            maxDepthReached(false),
            rootLimitReached(false),
            leaf(true),
            numNodes(0),
            node(nullptr),
            gp(nullptr){}

    QuadTree(Point<FLOAT> center);

    ~QuadTree(){
        deleteChildren();
    }

    bool IsRoot(){
        if (par)
            return false;
        else
            return true;
    }

    QuadTree* const getRoot();

    // Note: Call this function ONLY BEFORE creating an instance of tree
    static void setTreeParam(tree_param par){
       param = par;
    };

    bool Insert(std::shared_ptr<Node> n);
    bool Insert(std::shared_ptr<Node> n, std::unordered_set<QuadTree*>& quads);
    bool IsNotNew(std::shared_ptr<Node> n);
    bool Update(std::shared_ptr<Node> n);
    bool Update(std::shared_ptr<Node> n, std::unordered_set<QuadTree*>& quads);
    bool Remove(std::shared_ptr<Node> n, std::unordered_set<QuadTree*>& quads);

    void Update(std::shared_ptr<OnGPIS> _gp);
    std::shared_ptr<OnGPIS> const getGP(){return gp;}

    // Let's implement
    bool Remove(std::shared_ptr<Node> n);
    void QueryRange(AABB range, std::vector<std::shared_ptr<Node> >& nodes);
    void QueryNonEmptyLevelC(AABB range, std::vector<QuadTree*>& quads);
    void QueryNonEmptyLevelC(AABB range, std::vector<QuadTree*>& quads, std::vector<FLOAT>& sqdst);
    void QueryNonEmptyLevelC(AABB range, std::vector<QuadTree*>& quads, std::vector<std::vector<std::shared_ptr<Node> > >& nodes);

    int32_t getNodeCount(){return numNodes;}
    Point<FLOAT> getCenter(){return boundary.getCenter();}
    FLOAT getHalfLength(){return boundary.getHalfLength();}
    FLOAT getXMaxbound(){return boundary.getXMaxbound();}
    FLOAT getYMinbound(){return boundary.getYMinbound();}
    FLOAT getYMaxbound(){return boundary.getYMaxbound();}
    Point<FLOAT> getNW(){return boundary.getNW();}
    Point<FLOAT> getNE(){return boundary.getNE();}
    Point<FLOAT> getSW(){return boundary.getSW();}
    Point<FLOAT> getSE(){return boundary.getSE();}

    void getAllChildrenNonEmptyNodes(std::vector<std::shared_ptr<Node> >& nodes);

    void printNodes();
    void printBoundary();

};


#endif

