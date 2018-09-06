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
#include "strct.h"
#include "OnGPIS.h"

class AABB
{
    Point<float> center;
    float halfLength;
    float halfLengthSq;
    float xmin;
    float xmax;
    float ymin;
    float ymax;

    Point<float> ptNW;
    Point<float> ptNE;
    Point<float> ptSW;
    Point<float> ptSE;
public:
    AABB(){
        halfLength = 0.0;
        halfLengthSq = 0.0;
        xmin = 0.0;
        xmax = 0.0;
        ymin = 0.0;
        ymax = 0.0;
    }
    AABB(Point<float> _center, float _halfLength) {
        center = _center;
        halfLength = _halfLength;
        halfLengthSq = halfLength*halfLength;
        xmin = center.x - halfLength;
        xmax = center.x + halfLength;
        ymin = center.y - halfLength;
        ymax = center.y + halfLength;
        ptNW = Point<float>(xmin,ymax);
        ptNE = Point<float>(xmax,ymax);
        ptSW = Point<float>(xmin,ymin);
        ptSE = Point<float>(xmax,ymin);
    }
    AABB(float x, float y, float _halfLength) {
        center = Point<float>(x,y);
        halfLength = _halfLength;
        halfLengthSq = halfLength*halfLength;
        xmin = center.x - halfLength;
        xmax = center.x + halfLength;
        ymin = center.y - halfLength;
        ymax = center.y + halfLength;
        ptNW = Point<float>(xmin,ymax);
        ptNE = Point<float>(xmax,ymax);
        ptSW = Point<float>(xmin,ymin);
        ptSE = Point<float>(xmax,ymin);
    }

    const Point<float> getCenter(){return center;}
    float getHalfLength(){return halfLength;}
    float getHalfLengthSq(){return halfLengthSq;}
    float getXMinbound(){return xmin;}
    float getXMaxbound(){return xmax;}
    float getYMinbound(){return ymin;}
    float getYMaxbound(){return ymax;}
    const Point<float>& getNW(){return ptNW;}
    const Point<float>& getNE(){return ptNE;}
    const Point<float>& getSW(){return ptSW;}
    const Point<float>& getSE(){return ptSE;}

    bool containsPoint(Point<float> pt) {
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

    QuadTree(Point<float> center);

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

    bool Remove(std::shared_ptr<Node> n);
    void QueryRange(AABB range, std::vector<std::shared_ptr<Node> >& nodes);
    void QueryNonEmptyLevelC(AABB range, std::vector<QuadTree*>& quads);
    void QueryNonEmptyLevelC(AABB range, std::vector<QuadTree*>& quads, std::vector<float>& sqdst);
    void QueryNonEmptyLevelC(AABB range, std::vector<QuadTree*>& quads, std::vector<std::vector<std::shared_ptr<Node> > >& nodes);

    int32_t getNodeCount(){return numNodes;}
    Point<float> getCenter(){return boundary.getCenter();}
    float getHalfLength(){return boundary.getHalfLength();}
    float getXMaxbound(){return boundary.getXMaxbound();}
    float getYMinbound(){return boundary.getYMinbound();}
    float getYMaxbound(){return boundary.getYMaxbound();}
    Point<float> getNW(){return boundary.getNW();}
    Point<float> getNE(){return boundary.getNE();}
    Point<float> getSW(){return boundary.getSW();}
    Point<float> getSE(){return boundary.getSE();}

    void getAllChildrenNonEmptyNodes(std::vector<std::shared_ptr<Node> >& nodes);

};

#endif
