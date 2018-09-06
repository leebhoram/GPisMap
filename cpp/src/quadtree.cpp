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

#include "quadtree.h"

#define EPS 1e-12

float sqdist(const Point<float>& pt1, const Point<float>& pt2)
{
    float dx = (pt1.x - pt2.x);
    float dy = (pt1.y - pt2.y);
    return dx*dx + dy*dy;
}

QuadTree::QuadTree(Point<float> c):
        northWest(0),
        northEast(0),
        southWest(0),
        southEast(0),
        par(0),
        maxDepthReached(false),
        rootLimitReached(false),
        leaf(true),
        numNodes(0),
        node(nullptr),
        gp(nullptr)
{
     boundary = AABB(c,QuadTree::param.initroot_halfleng);
}

QuadTree::QuadTree(AABB _boundary, QuadTree* const p )
        :northWest(0),
         northEast(0),
         southWest(0),
         southEast(0),
               par(0),
        maxDepthReached(false),
        rootLimitReached(false),
        leaf(true),
        numNodes(0),
        node(nullptr),
        gp(nullptr){
    boundary = _boundary;
    if (boundary.getHalfLength() < QuadTree::param.min_halfleng)
        maxDepthReached = true;
    if (boundary.getHalfLength() > QuadTree::param.max_halfleng)
        rootLimitReached = true;
    if (p!=0)
        par = p;
}

QuadTree::QuadTree(AABB _boundary,  QuadTree* const ch,  int child_type)
        :northWest(0),
         northEast(0),
         southWest(0),
         southEast(0),
               par(0),
        maxDepthReached(false),
        rootLimitReached(false),
        numNodes(0),
        node(nullptr),
        gp(nullptr){
    boundary = _boundary;
    if (boundary.getHalfLength() < QuadTree::param.min_halfleng)
        maxDepthReached = true;
    if (boundary.getHalfLength() > QuadTree::param.max_halfleng)
        rootLimitReached = true;
    if (child_type == 0)
    {
        leaf = true;
    }
    else
    {
        leaf = false;
        SubdivideExcept(child_type);
        if (child_type == CHILD_TYPE_NW)
            northWest = ch;
        if (child_type == CHILD_TYPE_NE)
            northEast = ch;
        if (child_type == CHILD_TYPE_SW)
            southWest = ch;
        if (child_type == CHILD_TYPE_SE)
            southEast = ch;
    }
}

void QuadTree::deleteChildren()
{
    if (northWest) {delete northWest; northWest = 0;}
    if (northEast) {delete northEast; northEast = 0;}
    if (southWest) {delete southWest; southWest = 0;}
    if (southEast) {delete southEast; southEast = 0;}
    leaf = true;
}

QuadTree* const QuadTree::getRoot(){
    QuadTree* p = this;
    QuadTree* p1 = p->getParent();
    while (p1!=0){
        p = p1;
        p1 = p->getParent();
    }
    return p;
}

bool QuadTree::InsertToParent(std::shared_ptr<Node> n){
    float l = getHalfLength();
    Point<float> c = getCenter();

    // Find out what type the current node is
    const Point<float> np = n->getPos();

    Point<float> par_c;
    int childType = 0;
    if (np.x < c.x && np.y > c.y){
        childType = CHILD_TYPE_SE;
        par_c.x = c.x - l;
        par_c.y = c.y + l;
    }
    if (np.x > c.x && np.y > c.y){
        childType = CHILD_TYPE_SW;
        par_c.x = c.x + l;
        par_c.y = c.y + l;
    }
    if (np.x < c.x && np.y < c.y){
        childType = CHILD_TYPE_NE;
        par_c.x = c.x - l;
        par_c.y = c.y - l;
    }
    if (np.x > c.x && np.y < c.y){
        childType = CHILD_TYPE_NW;
        par_c.x = c.x + l;
        par_c.y = c.y - l;
    }

    AABB parbb(par_c,2.0*l);
    par = new QuadTree(parbb,this,childType);
    return par->Insert(n);
}

bool QuadTree::Insert(std::shared_ptr<Node> n){

    // Ignore objects that do not belong in this quad tree
    if (!boundary.containsPoint(n->getPos())){
        if (getParent() == 0) {
            if (rootLimitReached)  return false;
            else                   return InsertToParent(n);
        }
        return false; // object cannot be added
    }

    if (maxDepthReached){
        if ( node == nullptr) {// If this is the first point in this quad tree, add the object here
            node = n;
            numNodes = 1;
            return true;
        }
        else // no more points accepted at this resolution
            return false;
    }

    if (IsLeaf()){

        if (boundary.getHalfLength() > QuadTree::param.cluster_halfleng){
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
            //numNodes = 0;
            if (sqdist(node->getPos(), n->getPos()) < QuadTree::param.min_halfleng_sqr){
                return false;
            }

            Subdivide();
            if (northWest->Insert(node)){
                ;
            }
            else if (northEast->Insert(node)){
                ;
            }
            else if (southWest->Insert(node)) {
                ;
            }
            else if (southEast->Insert(node)) {
                ;
            }
            node = nullptr;
        }
    }

    if (northWest->Insert(n)) {updateCount(); return true;}
    if (northEast->Insert(n)) {updateCount(); return true;}
    if (southWest->Insert(n)) {updateCount(); return true;}
    if (southEast->Insert(n)) {updateCount(); return true;}

    return false;
}

bool QuadTree::Insert(std::shared_ptr<Node> n,  std::unordered_set<QuadTree*>& quads)
{
     // Ignore objects that do not belong in this quad tree
    if (!boundary.containsPoint(n->getPos())){
        if (getParent() == 0) {
            if (rootLimitReached)  return false;
            else                   return InsertToParent(n);
        }
        return false; // object cannot be added
    }

    if (maxDepthReached){
        if ( node == nullptr) {// If this is the first point in this quad tree, add the object here
            node = n;
            numNodes = 1;
            if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
                quads.insert(this);
            return true;
        }
        else // no more points accepted at this resolution
            return false;
    }

    if (IsLeaf()){

        if (boundary.getHalfLength() > QuadTree::param.cluster_halfleng){
            Subdivide();
        }
        else{
            if (node == nullptr)
            // If this is the first point in this quad tree, add the object here
            {
                node = n;
                numNodes = 1;
                if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
                    quads.insert(this);
                return true;
            }

            // Otherwise, subdivide and then add the point to whichever node will accept it
            if (sqdist(node->getPos(), n->getPos()) < QuadTree::param.min_halfleng_sqr){
                return false;
            }

            Subdivide();
            if (northWest->Insert(node,quads)){
                ;
            }
            else if (northEast->Insert(node,quads)){
                ;
            }
            else if (southWest->Insert(node,quads)) {
                ;
            }
            else if (southEast->Insert(node,quads)) {
                ;
            }
            node = nullptr;
        }
    }

    if (northWest->Insert(n,quads)) {
        if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
            quads.insert(this);
        updateCount();
        return true;
    }

    if (northEast->Insert(n,quads)) {
        if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
            quads.insert(this);
        updateCount();
        return true;
    }

    if (southWest->Insert(n,quads)) {
        if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
            quads.insert(this);
        updateCount();
        return true;
    }

    if (southEast->Insert(n,quads)) {
        if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
            quads.insert(this);
        updateCount();
        return true;
    }

}

void QuadTree::updateCount()
{
    if (leaf==false){
        numNodes = 0;
        numNodes += northWest->getNodeCount();
        numNodes += northEast->getNodeCount();
        numNodes += southWest->getNodeCount();
        numNodes += southEast->getNodeCount();
    }
}

bool QuadTree::IsNotNew(std::shared_ptr<Node> n)
{
    if (!boundary.containsPoint(n->getPos())){
        return false; // object cannot be added
    }

    if (IsEmptyLeaf())
        return false ;

    if (!IsEmpty() && (sqdist(node->getPos(), n->getPos()) < QuadTree::param.min_halfleng_sqr))
    {
        return true;
    }

    if (IsLeaf())
        return false;

    if (northWest->IsNotNew(n)) return true;
    if (northEast->IsNotNew(n)) return true;
    if (southWest->IsNotNew(n)) return true;
    if (southEast->IsNotNew(n)) return true;

    return false;
}

bool QuadTree::Remove(std::shared_ptr<Node> n){
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

    bool res = northWest->Remove(n);
    res |= northEast->Remove(n);
    res |= southWest->Remove(n);
    res |= southEast->Remove(n);

    if (res)
    {
        bool res2 = northWest->IsEmptyLeaf();
        res2 &= northEast->IsEmptyLeaf();
        res2 &= southWest->IsEmptyLeaf();
        res2 &= southEast->IsEmptyLeaf();
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

bool QuadTree::Remove(std::shared_ptr<Node> n,std::unordered_set<QuadTree*>& quads){
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

    bool res = northWest->Remove(n,quads);
    if (!res) res |= northEast->Remove(n,quads);
    if (!res) res |= southWest->Remove(n,quads);
    if (!res) res |= southEast->Remove(n,quads);

    if (res)
    {
        bool res2 = northWest->IsEmptyLeaf();
        res2 &= northEast->IsEmptyLeaf();
        res2 &= southWest->IsEmptyLeaf();
        res2 &= southEast->IsEmptyLeaf();
        if (res2)
        {
            quads.erase(northWest);
            quads.erase(northEast);
            quads.erase(southWest);
            quads.erase(southEast);
            deleteChildren();
            leaf = true;
            numNodes = 0;
        }
    }
    updateCount();

    return res;
}

void QuadTree::Update(std::shared_ptr<OnGPIS> _gp)
{
    gp = _gp;
}

bool QuadTree::Update(std::shared_ptr<Node> n){
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

    if(northWest->Update(n)) return true;
    if(northEast->Update(n)) return true;
    if(southWest->Update(n)) return true;
    if(southEast->Update(n)) return true;

}

bool QuadTree::Update(std::shared_ptr<Node> n, std::unordered_set<QuadTree*>& quads){
    // Ignore objects that do not belong in this quad tree
    if (!boundary.containsPoint(n->getPos())){
        return false; // object cannot be added
    }

    if (IsEmptyLeaf())
        return false ;

    if (!IsEmpty() && (sqdist(node->getPos(), n->getPos()) < EPS))
    {
        node = n;
        if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
            quads.insert(this);
        return true;
    }

    if (IsLeaf())
        return false;

    if(northWest->Update(n,quads)){
        if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
            quads.insert(this);
        return true;
    }
    if(northEast->Update(n,quads)) {
        if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
            quads.insert(this);
        return true;
    }
    if(southWest->Update(n,quads)) {
        if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
            quads.insert(this);
        return true;
    }
    if(southEast->Update(n,quads)) {
        if (fabs(getHalfLength()-QuadTree::param.cluster_halfleng) < 1e-3)
            quads.insert(this);
        return true;
    }

}

void QuadTree::Subdivide()
{
    float l = boundary.getHalfLength()*0.5;
    Point<float> c = boundary.getCenter();
    Point<float> nw_c = Point<float>(c.x-l,c.y+l);
    AABB nw(nw_c,l);
    northWest = new QuadTree(nw,this);

    Point<float> ne_c = Point<float>(c.x+l,c.y+l);
    AABB ne(ne_c,l);
    northEast =  new QuadTree(ne,this);

    Point<float> sw_c = Point<float>(c.x-l,c.y-l);
    AABB sw(sw_c,l);
    southWest =  new QuadTree(sw,this);

    Point<float> se_c = Point<float>(c.x+l,c.y-l);
    AABB se(se_c,l);
    southEast =  new QuadTree(se,this);

    leaf = false;

    return;
}

void QuadTree::SubdivideExcept(int childType)
{
    float l = boundary.getHalfLength()*0.5;
    Point<float> c = boundary.getCenter();

    if (childType != CHILD_TYPE_NW)
    {
        Point<float> nw_c = Point<float>(c.x-l,c.y+l);
        AABB nw(nw_c,l);
        northWest = new QuadTree(nw,this);
    }

    if (childType != CHILD_TYPE_NE)
    {
        Point<float> ne_c = Point<float>(c.x+l,c.y+l);
        AABB ne(ne_c,l);
        northEast = new QuadTree(ne,this);
    }

    if (childType != CHILD_TYPE_SW)
    {
        Point<float> sw_c = Point<float>(c.x-l,c.y-l);
        AABB sw(sw_c,l);
        southWest = new QuadTree(sw,this);
    }

    if (childType != CHILD_TYPE_SE)
    {
        Point<float> se_c = Point<float>(c.x+l,c.y-l);
        AABB se(se_c,l);
        southEast = new QuadTree(se,this);
    }

    leaf = false;
}

 // Find all points that appear within a range
void QuadTree::QueryRange(AABB range, std::vector<std::shared_ptr<Node> >& nodes)
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
    northWest->QueryRange(range,nodes);
    northEast->QueryRange(range,nodes);
    southWest->QueryRange(range,nodes);
    southEast->QueryRange(range,nodes);

    return ;
 }

void QuadTree::getAllChildrenNonEmptyNodes(std::vector<std::shared_ptr<Node> >& nodes)
{
     if (IsEmptyLeaf())
         return;

     if (IsLeaf())
     {
         nodes.push_back(node);
         return;
     }

     northWest->getAllChildrenNonEmptyNodes(nodes);
     northEast->getAllChildrenNonEmptyNodes(nodes);
     southWest->getAllChildrenNonEmptyNodes(nodes);
     southEast->getAllChildrenNonEmptyNodes(nodes);
     return;
}

void QuadTree::QueryNonEmptyLevelC(AABB range, std::vector<QuadTree*>& quads)
{
    // Automatically abort if the range does not intersect this quad
    if (!boundary.intersectsAABB(range) || IsEmptyLeaf()){
        return; // empty list
    }

    if (IsLeaf()){ // no children
        if (boundary.getHalfLength() > (QuadTree::param.cluster_halfleng+0.0001)){
            return;
        }
    }

    if (boundary.getHalfLength() > (QuadTree::param.cluster_halfleng+0.001)){
        // Otherwise, add the points from the children
        northWest->QueryNonEmptyLevelC(range,quads);
        northEast->QueryNonEmptyLevelC(range,quads);
        southWest->QueryNonEmptyLevelC(range,quads);
        southEast->QueryNonEmptyLevelC(range,quads);
    }
    else
    {
        quads.push_back(this);
    }

    return ;
}

void QuadTree::QueryNonEmptyLevelC(AABB range, std::vector<QuadTree*>& quads, std::vector<float>& sqdst)
{

    // Automatically abort if the range does not intersect this quad
    if (!boundary.intersectsAABB(range) || IsEmptyLeaf()){
        return; // empty list
    }

    if (IsLeaf()){ // no children
        if (boundary.getHalfLength() > QuadTree::param.cluster_halfleng+0.001){
            return;
        }
    }

    if (boundary.getHalfLength() > QuadTree::param.cluster_halfleng+0.001){
        // Otherwise, add the points from the children
        northWest->QueryNonEmptyLevelC(range,quads,sqdst);
        northEast->QueryNonEmptyLevelC(range,quads,sqdst);
        southWest->QueryNonEmptyLevelC(range,quads,sqdst);
        southEast->QueryNonEmptyLevelC(range,quads,sqdst);
    }
    else
    {
        sqdst.push_back(sqdist(getCenter(),range.getCenter()));
        quads.push_back(this);
    }

    return ;
}
