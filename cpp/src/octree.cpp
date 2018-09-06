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

#include "octree.h"

#define EPS 1e-12

static float sqdist(const Point3<float>& pt1, const Point3<float>& pt2)
{
    float dx = (pt1.x - pt2.x);
    float dy = (pt1.y - pt2.y);
    float dz = (pt1.z - pt2.z);

    return dx*dx + dy*dy + dz*dz;
}

OcTree::OcTree(Point3<float> c)
        :northWestFront(0),
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
        gp(nullptr){
    boundary = AABB3(c,OcTree::param.initroot_halfleng);
}

OcTree::OcTree(AABB3 _boundary, OcTree* const p )
        :northWestFront(0),
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
        gp(nullptr){
    boundary = _boundary;
    if (boundary.getHalfLength() < OcTree::param.min_halfleng)
        maxDepthReached = true;
    if (boundary.getHalfLength() > OcTree::param.max_halfleng)
        rootLimitReached = true;
    if (p!=0)
        par = p;
}

OcTree::OcTree(AABB3 _boundary,  OcTree* const ch,  int child_type)
        :northWestFront(0),
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
        numNodes(0),
        node(nullptr),
        closestChildNode(nullptr),
        gp(nullptr){
    boundary = _boundary;
    if (boundary.getHalfLength() < OcTree::param.min_halfleng)
        maxDepthReached = true;
    if (boundary.getHalfLength() > OcTree::param.max_halfleng)
        rootLimitReached = true;
    if (child_type == 0)
    {
        leaf = true;
    }
    else
    {
        leaf = false;
        SubdivideExcept(child_type);
        if (child_type == CHILD_TYPE_NWF)
            northWestFront = ch;
        if (child_type == CHILD_TYPE_NEF)
            northEastFront = ch;
        if (child_type == CHILD_TYPE_SWF)
            southWestFront = ch;
        if (child_type == CHILD_TYPE_SEF)
            southEastFront = ch;
        if (child_type == CHILD_TYPE_NWB)
            northWestBack = ch;
        if (child_type == CHILD_TYPE_NEB)
            northEastBack = ch;
        if (child_type == CHILD_TYPE_SWB)
            southWestBack = ch;
        if (child_type == CHILD_TYPE_SEB)
            southEastBack = ch;
    }

    Point3<float> c = boundary.getCenter();
}

void OcTree::deleteChildren()
{
    if (northWestFront) {delete northWestFront; northWestFront = 0;}
    if (northEastFront) {delete northEastFront; northEastFront = 0;}
    if (southWestFront) {delete southWestFront; southWestFront = 0;}
    if (southEastFront) {delete southEastFront; southEastFront = 0;}
    if (northWestBack) {delete northWestBack; northWestBack = 0;}
    if (northEastBack) {delete northEastBack; northEastBack = 0;}
    if (southWestBack) {delete southWestBack; southWestBack = 0;}
    if (southEastBack) {delete southEastBack; southEastBack = 0;}
}

OcTree* const OcTree::getRoot(){
    OcTree* p = this;
    OcTree* p1 = p->getParent();
    while (p1!=0){
        p = p1;
        p1 = p->getParent();
    }
    return p;
}

bool OcTree::InsertToParent(std::shared_ptr<Node3> n){
    float l = getHalfLength();
    Point3<float> c = getCenter();

    // Find out what type the current node is
    const Point3<float> np = n->getPos();

    Point3<float> par_c;
    int childType = 0;
    if (np.x < c.x && np.y > c.y && np.z > c.z){
        childType = CHILD_TYPE_SEB;
        par_c.x = c.x - l;
        par_c.y = c.y + l;
        par_c.z = c.z + l;
    }
    if (np.x > c.x && np.y > c.y && np.z > c.z){
        childType = CHILD_TYPE_SWB;
        par_c.x = c.x + l;
        par_c.y = c.y + l;
        par_c.z = c.z + l;
    }
    if (np.x < c.x && np.y < c.y && np.z > c.z){
        childType = CHILD_TYPE_NEB;
        par_c.x = c.x - l;
        par_c.y = c.y - l;
        par_c.z = c.z + l;
    }
    if (np.x > c.x && np.y < c.y && np.z > c.z){
        childType = CHILD_TYPE_NWB;
        par_c.x = c.x + l;
        par_c.y = c.y - l;
        par_c.z = c.z + l;
    }
    if (np.x < c.x && np.y > c.y && np.z < c.z){
        childType = CHILD_TYPE_SEF;
        par_c.x = c.x - l;
        par_c.y = c.y + l;
        par_c.z = c.z - l;
    }
    if (np.x > c.x && np.y > c.y && np.z < c.z){
        childType = CHILD_TYPE_SWF;
        par_c.x = c.x + l;
        par_c.y = c.y + l;
        par_c.z = c.z - l;
    }
    if (np.x < c.x && np.y < c.y && np.z < c.z){
        childType = CHILD_TYPE_NEF;
        par_c.x = c.x - l;
        par_c.y = c.y - l;
        par_c.z = c.z - l;
    }
    if (np.x > c.x && np.y < c.y && np.z < c.z){
        childType = CHILD_TYPE_NWF;
        par_c.x = c.x + l;
        par_c.y = c.y - l;
        par_c.z = c.z - l;
    }

    AABB3 parbb(par_c,2.0*l);
    par = new OcTree(parbb,this,childType);
    return par->Insert(n);
}

bool OcTree::Insert(std::shared_ptr<Node3> n){

    // Ignore objects that do not belong in this quad tree
    if (!boundary.containsPoint(n->getPos())){
        if (getParent() == 0) {
            if (rootLimitReached)  return false;
            else                   return InsertToParent(n);
        }
        return false; // object cannot be added
    }

    if (maxDepthReached){
        if (node == nullptr) {// If this is the first point in this quad tree, add the object here
            node = n;
            numNodes = 1;
            return true;
        }
        else // no more points accepted at this resolution
            return false;
    }

    if (IsLeaf()){

        if (boundary.getHalfLength() > OcTree::param.cluster_halfleng){
            Subdivide();
        }
        else{
            if (node == nullptr)
            // If this is the first point in this quad tree, add the object here
            {
                node = n;
                numNodes = 1;
                return true;
            }

            // Otherwise, subdivide and then add the point to whichever node will accept it
            if (sqdist(node->getPos(), n->getPos()) < OcTree::param.min_halfleng_sqr){
                return false;
            }

            Subdivide();
            if (northWestFront->Insert(node)){
                ;
            }
            else if (northEastFront->Insert(node)){
                ;
            }
            else if (southWestFront->Insert(node)) {
                ;
            }
            else if (southEastFront->Insert(node)) {
                ;
            }
            else if (northWestBack->Insert(node)){
                ;
            }
            else if (northEastBack->Insert(node)){
                ;
            }
            else if (southWestBack->Insert(node)) {
                ;
            }
            else if (southEastBack->Insert(node)) {
                ;
            }
            node = nullptr;
        }
    }

    if (northWestFront->Insert(n)) {updateCount(); return true;}
    if (northEastFront->Insert(n)) {updateCount(); return true;}
    if (southWestFront->Insert(n)) {updateCount(); return true;}
    if (southEastFront->Insert(n)) {updateCount(); return true;}
    if (northWestBack->Insert(n)) {updateCount(); return true;}
    if (northEastBack->Insert(n)) {updateCount(); return true;}
    if (southWestBack->Insert(n)) {updateCount(); return true;}
    if (southEastBack->Insert(n)) {updateCount(); return true;}

    return false;
}

bool OcTree::Insert(std::shared_ptr<Node3> n, std::unordered_set<OcTree*>& quads){
    // Ignore objects that do not belong in this quad tree
    if (!boundary.containsPoint(n->getPos())){
        if (getParent() == 0) {
            if (rootLimitReached)  return false;
            else                   return InsertToParent(n);
        }
        return false; // object cannot be added
    }

    if (maxDepthReached){
        if (node == nullptr) {// If this is the first point in this quad tree, add the object here
            node = n;
            numNodes = 1;
            return true;
        }
        else // no more points accepted at this resolution
            return false;
    }

    if (IsLeaf()){

        if (boundary.getHalfLength() > OcTree::param.cluster_halfleng){
            Subdivide();
        }
        else{
            if (node == nullptr)
            {
                node = n;
                numNodes = 1;
                if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-6)
                    quads.insert(this);
                return true;
            }

             // Otherwise, subdivide and then add the point to whichever node will accept it
            if (sqdist(node->getPos(), n->getPos()) < OcTree::param.min_halfleng_sqr){
                return false;
            }

            Subdivide();
            if (northWestFront->Insert(node,quads)){
                ;
            }
            else if (northEastFront->Insert(node,quads)){
                ;
            }
            else if (southWestFront->Insert(node,quads)) {
                ;
            }
            else if (southEastFront->Insert(node,quads)) {
                ;
            }
            else if (northWestBack->Insert(node,quads)){
                ;
            }
            else if (northEastBack->Insert(node,quads)){
                ;
            }
            else if (southWestBack->Insert(node,quads)) {
                ;
            }
            else if (southEastBack->Insert(node,quads)) {
                ;
            }
            node = nullptr;
        }
    }

    if (northWestFront->Insert(n,quads)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-6)
            quads.insert(this);
        updateCount(); return true;
    }

    if (northEastFront->Insert(n,quads)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-6)
            quads.insert(this);
        updateCount(); return true;
    }

    if (southWestFront->Insert(n,quads)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-6)
            quads.insert(this);
        updateCount(); return true;
    }

    if (southEastFront->Insert(n,quads)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-6)
            quads.insert(this);
        updateCount(); return true;
    }

    if (northWestBack->Insert(n,quads)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-6)
            quads.insert(this);
        updateCount(); return true;
    }

    if (northEastBack->Insert(n,quads)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-6)
            quads.insert(this);
        updateCount(); return true;
    }

    if (southWestBack->Insert(n,quads)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-6)
            quads.insert(this);
        updateCount(); return true;
    }

    if (southEastBack->Insert(n,quads)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-6)
            quads.insert(this);
        updateCount(); return true;
    }

    return false;

}

void OcTree::updateCount()
{
    if (leaf==false){
        numNodes = 0;
        numNodes += northWestFront->getNodeCount();
        numNodes += northEastFront->getNodeCount();
        numNodes += southWestFront->getNodeCount();
        numNodes += southEastFront->getNodeCount();
        numNodes += northWestBack->getNodeCount();
        numNodes += northEastBack->getNodeCount();
        numNodes += southWestBack->getNodeCount();
        numNodes += southEastBack->getNodeCount();
    }
}

bool OcTree::IsNotNew(std::shared_ptr<Node3> n)
{
    if (!boundary.containsPoint(n->getPos())){
        return false; // object cannot be added
    }

    if (IsEmptyLeaf())
        return false ;

    if (!IsEmpty() && (sqdist(node->getPos(), n->getPos()) < OcTree::param.min_halfleng_sqr))
    {
        return true;
    }

    if (IsLeaf())
        return false;

    if (northWestFront->IsNotNew(n)) return true;
    if (northEastFront->IsNotNew(n)) return true;
    if (southWestFront->IsNotNew(n)) return true;
    if (southEastFront->IsNotNew(n)) return true;
    if (northWestBack->IsNotNew(n)) return true;
    if (northEastBack->IsNotNew(n)) return true;
    if (southWestBack->IsNotNew(n)) return true;
    if (southEastBack->IsNotNew(n)) return true;

    return false;
}

bool OcTree::Remove(std::shared_ptr<Node3> n){
    // Ignore objects that do not belong in this quad tree
    if (!boundary.containsPoint(n->getPos())){
        return false; // object cannot be added
    }

    if (IsEmptyLeaf())
        return false ;

    if (!IsEmpty() && (sqdist(node->getPos(), n->getPos()) < EPS))
    {
        node = nullptr;
        numNodes = 0;
        return true;
    }

    if (IsLeaf())
        return false;

    bool res = northWestFront->Remove(n);
    res |= northEastFront->Remove(n);
    res |= southWestFront->Remove(n);
    res |= southEastFront->Remove(n);
    res |= northWestBack->Remove(n);
    res |= northEastBack->Remove(n);
    res |= southWestBack->Remove(n);
    res |= southEastBack->Remove(n);

    if (res)
    {
        bool res2 = northWestFront->IsEmptyLeaf();
        res2 &= northEastFront->IsEmptyLeaf();
        res2 &= southWestFront->IsEmptyLeaf();
        res2 &= southEastFront->IsEmptyLeaf();
        res2 &= northWestBack->IsEmptyLeaf();
        res2 &= northEastBack->IsEmptyLeaf();
        res2 &= southWestBack->IsEmptyLeaf();
        res2 &= southEastBack->IsEmptyLeaf();
        if (res2)
        {
            deleteChildren();
            leaf = true;
            numNodes = 0;
        }
    }
    updateCount();

    return res;
}

bool OcTree::Remove(std::shared_ptr<Node3> n,std::unordered_set<OcTree*>& octs){
    // Ignore objects that do not belong in this quad tree
    if (!boundary.containsPoint(n->getPos())){
        return false; // object cannot be added
    }

    if (IsEmptyLeaf())
        return false ;

    if (!IsEmpty() && (sqdist(node->getPos(), n->getPos()) < EPS))
    {
        node = nullptr;
        numNodes = 0;
        return true;
    }

    if (IsLeaf())
        return false;

    bool res = northWestFront->Remove(n,octs);
    if (!res) res |= northEastFront->Remove(n,octs);
    if (!res) res |= southWestFront->Remove(n,octs);
    if (!res) res |= southEastFront->Remove(n,octs);
    if (!res) res |= northWestBack->Remove(n,octs);
    if (!res) res |= northEastBack->Remove(n,octs);
    if (!res) res |= southWestBack->Remove(n,octs);
    if (!res) res |= southEastBack->Remove(n,octs);

    if (res)
    {
        bool res1 = northWestFront->IsEmptyLeaf();
        bool res2 = northEastFront->IsEmptyLeaf();
        bool res3 = southWestFront->IsEmptyLeaf();
        bool res4 = southEastFront->IsEmptyLeaf();
        bool res5 = northWestBack->IsEmptyLeaf();
        bool res6 = northEastBack->IsEmptyLeaf();
        bool res7 = southWestBack->IsEmptyLeaf();
        bool res8 = southEastBack->IsEmptyLeaf();
        if (res1 & res2 & res3 & res4 & res5 & res6 & res7 & res8)
        {
            octs.erase(northWestFront);
            octs.erase(northEastFront);
            octs.erase(southWestFront);
            octs.erase(southEastFront);
            octs.erase(northWestBack);
            octs.erase(northEastBack);
            octs.erase(southWestBack);
            octs.erase(southEastBack);
            deleteChildren();
            leaf = true;
            numNodes = 0;
        }
    }
    updateCount();

    return res;
}

void OcTree::Update(std::shared_ptr<OnGPIS> _gp)
{
    gp = _gp;
}

bool OcTree::Update(std::shared_ptr<Node3> n){
    // Ignore objects that do not belong in this quad tree
    if (!boundary.containsPoint(n->getPos())){
        return false; // object cannot be added
    }

    if (IsEmptyLeaf())
        return false ;

    if (!IsEmpty() && (sqdist(node->getPos(), n->getPos()) < EPS))
    {
        node = n;
        return true;
    }

    if (IsLeaf())
        return false;

    if(northWestFront->Update(n)) return true;
    if(northEastFront->Update(n)) return true;
    if(southWestFront->Update(n)) return true;
    if(southEastFront->Update(n)) return true;

    if(northWestBack->Update(n)) return true;
    if(northEastBack->Update(n)) return true;
    if(southWestBack->Update(n)) return true;
    if(southEastBack->Update(n)) return true;

    return false;

}

bool OcTree::Update(std::shared_ptr<Node3> n, std::unordered_set<OcTree*>& octs){
    // Ignore objects that do not belong in this quad tree
    if (!boundary.containsPoint(n->getPos())){
        return false; // object cannot be added
    }

    if (IsEmptyLeaf())
        return false ;

    if (!IsEmpty() && (sqdist(node->getPos(), n->getPos()) < EPS))
    {
        node = n;
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-3)
            octs.insert(this);
        return true;
    }

    if (IsLeaf())
        return false;

    if(northWestFront->Update(n,octs)){
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-3)
            octs.insert(this);
        return true;
    }
    if(northEastFront->Update(n,octs)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-3)
            octs.insert(this);
        return true;
    }
    if(southWestFront->Update(n,octs)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-3)
            octs.insert(this);
        return true;
    }
    if(southEastFront->Update(n,octs)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-3)
            octs.insert(this);
        return true;
    }

    if(northWestBack->Update(n,octs)){
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-3)
            octs.insert(this);
        return true;
    }
    if(northEastBack->Update(n,octs)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-3)
            octs.insert(this);
        return true;
    }
    if(southWestBack->Update(n,octs)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-3)
            octs.insert(this);
        return true;
    }
    if(southEastBack->Update(n,octs)) {
        if (fabs(getHalfLength()-OcTree::param.cluster_halfleng) < 1e-3)
            octs.insert(this);
        return true;
    }

    return false;
}

void OcTree::Subdivide()
{
    float l = boundary.getHalfLength()*0.5;
    Point3<float> c = boundary.getCenter();

    Point3<float> nwf_c = Point3<float>(c.x-l,c.y+l,c.z+l);
    AABB3 nwf(nwf_c,l);
    northWestFront = new OcTree(nwf,this);

    Point3<float> nef_c = Point3<float>(c.x+l,c.y+l,c.z+l);
    AABB3 nef(nef_c,l);
    northEastFront =  new OcTree(nef,this);

    Point3<float> swf_c = Point3<float>(c.x-l,c.y-l,c.z+l);
    AABB3 swf(swf_c,l);
    southWestFront  =  new OcTree(swf,this);

    Point3<float> sef_c = Point3<float>(c.x+l,c.y-l,c.z+l);
    AABB3 sef(sef_c,l);
    southEastFront  =  new OcTree(sef,this);

    Point3<float> nwb_c = Point3<float>(c.x-l,c.y+l,c.z-l);
    AABB3 nwb(nwb_c,l);
    northWestBack = new OcTree(nwb,this);

    Point3<float> neb_c = Point3<float>(c.x+l,c.y+l,c.z-l);
    AABB3 neb(neb_c,l);
    northEastBack =  new OcTree(neb,this);

    Point3<float> swb_c = Point3<float>(c.x-l,c.y-l,c.z-l);
    AABB3 swb(swb_c,l);
    southWestBack  =  new OcTree(swb,this);

    Point3<float> seb_c = Point3<float>(c.x+l,c.y-l,c.z-l);
    AABB3 seb(seb_c,l);
    southEastBack  =  new OcTree(seb,this);

    leaf = false;

    return;
}

void OcTree::SubdivideExcept(int childType)
{
    float l = boundary.getHalfLength()*0.5;
    Point3<float> c = boundary.getCenter();

    if (childType != CHILD_TYPE_NWF)
    {
        Point3<float> nwf_c = Point3<float>(c.x-l,c.y+l,c.z+l);
        AABB3 nwf(nwf_c,l);
        northWestFront = new OcTree(nwf,this);
    }

    if (childType != CHILD_TYPE_NEF)
    {
        Point3<float> nef_c = Point3<float>(c.x+l,c.y+l,c.z+l);
        AABB3 nef(nef_c,l);
        northEastFront =  new OcTree(nef,this);
    }

    if (childType != CHILD_TYPE_SWF)
    {
         Point3<float> swf_c = Point3<float>(c.x-l,c.y-l,c.z+l);
        AABB3 swf(swf_c,l);
        southWestFront  =  new OcTree(swf,this);
    }

    if (childType != CHILD_TYPE_SEF)
    {
        Point3<float> sef_c = Point3<float>(c.x+l,c.y-l,c.z+l);
        AABB3 sef(sef_c,l);
        southEastFront  =  new OcTree(sef,this);
    }

    if (childType != CHILD_TYPE_NWB)
    {
        Point3<float> nwb_c = Point3<float>(c.x-l,c.y+l,c.z-l);
        AABB3 nwb(nwb_c,l);
        northWestBack = new OcTree(nwb,this);
    }

    if (childType != CHILD_TYPE_NEB)
    {
        Point3<float> neb_c = Point3<float>(c.x+l,c.y+l,c.z-l);
        AABB3 neb(neb_c,l);
        northEastBack =  new OcTree(neb,this);
    }

    if (childType != CHILD_TYPE_SWB)
    {
        Point3<float> swb_c = Point3<float>(c.x-l,c.y-l,c.z-l);
        AABB3 swb(swb_c,l);
        southWestBack  =  new OcTree(swb,this);
    }

    if (childType != CHILD_TYPE_SEB)
    {
        Point3<float> seb_c = Point3<float>(c.x+l,c.y-l,c.z-l);
        AABB3 seb(seb_c,l);
        southEastBack  =  new OcTree(seb,this);
    }

    leaf = false;
}

 // Find all points that appear within a range
void OcTree::QueryRange(AABB3 range, std::vector<std::shared_ptr<Node3> >& nodes)
{
    // Automatically abort if the range does not intersect this quad
    if (!boundary.intersectsAABB(range) || IsEmptyLeaf()){
        return; // empty list
    }

    // Check objects at this quad level
    if (IsLeaf()){
        if (sqdist(node->getPos(), range.getCenter()) <  range.getHalfLengthSq()){
            nodes.push_back(node);
        }
        return;
    }

    // Otherwise, add the points from the children
    northWestFront->QueryRange(range,nodes);
    northEastFront->QueryRange(range,nodes);
    southWestFront->QueryRange(range,nodes);
    southEastFront->QueryRange(range,nodes);

    northWestBack->QueryRange(range,nodes);
    northEastBack->QueryRange(range,nodes);
    southWestBack->QueryRange(range,nodes);
    southEastBack->QueryRange(range,nodes);

    return ;
 }

void OcTree::getAllChildrenNonEmptyNodes(std::vector<std::shared_ptr<Node3> >& nodes)
{
     if (IsEmptyLeaf())
         return;

     if (IsLeaf())
     {
         nodes.push_back(node);
         return;
     }

     northWestFront->getAllChildrenNonEmptyNodes(nodes);
     northEastFront->getAllChildrenNonEmptyNodes(nodes);
     southWestFront->getAllChildrenNonEmptyNodes(nodes);
     southEastFront->getAllChildrenNonEmptyNodes(nodes);
     northWestBack->getAllChildrenNonEmptyNodes(nodes);
     northEastBack->getAllChildrenNonEmptyNodes(nodes);
     southWestBack->getAllChildrenNonEmptyNodes(nodes);
     southEastBack->getAllChildrenNonEmptyNodes(nodes);

     return;
}

void OcTree::QueryNonEmptyLevelC(AABB3 range, std::vector<OcTree*>& octs)
{
    // Automatically abort if the range does not intersect this quad
    if (!boundary.intersectsAABB(range) || IsEmptyLeaf()){
        return; // empty list
    }

    if (IsLeaf()){ // no children
        if (boundary.getHalfLength() > (OcTree::param.cluster_halfleng+0.0001)){
            return;
        }
    }

    if (boundary.getHalfLength() > (OcTree::param.cluster_halfleng+0.001)){
        // Otherwise, add the points from the children
        northWestFront->QueryNonEmptyLevelC(range,octs);
        northEastFront->QueryNonEmptyLevelC(range,octs);
        southWestFront->QueryNonEmptyLevelC(range,octs);
        southEastFront->QueryNonEmptyLevelC(range,octs);
        northWestBack->QueryNonEmptyLevelC(range,octs);
        northEastBack->QueryNonEmptyLevelC(range,octs);
        southWestBack->QueryNonEmptyLevelC(range,octs);
        southEastBack->QueryNonEmptyLevelC(range,octs);
    }
    else
    {
        octs.push_back(this);
    }

    return ;
}

void OcTree::QueryNonEmptyLevelC(AABB3 range, std::vector<OcTree*>& octs, std::vector<float>& sqdst)
{

    // Automatically abort if the range does not intersect this quad
    if (!boundary.intersectsAABB(range) || IsEmptyLeaf()){
        return; // empty list
    }

    if (IsLeaf()){ // no children
        if (boundary.getHalfLength() > OcTree::param.cluster_halfleng+0.001){
            return;
        }
    }

    if (boundary.getHalfLength() > OcTree::param.cluster_halfleng+0.001){
        // Otherwise, add the points from the children
        northWestFront->QueryNonEmptyLevelC(range,octs,sqdst);
        northEastFront->QueryNonEmptyLevelC(range,octs,sqdst);
        southWestFront->QueryNonEmptyLevelC(range,octs,sqdst);
        southEastFront->QueryNonEmptyLevelC(range,octs,sqdst);
        northWestBack->QueryNonEmptyLevelC(range,octs,sqdst);
        northEastBack->QueryNonEmptyLevelC(range,octs,sqdst);
        southWestBack->QueryNonEmptyLevelC(range,octs,sqdst);
        southEastBack->QueryNonEmptyLevelC(range,octs,sqdst);
    }
    else
    {
        sqdst.push_back(sqdist(getCenter(),range.getCenter()));
        octs.push_back(this);
    }

    return ;
}
