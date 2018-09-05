// strct.h
#include <vector>
#include <memory>
#include "strct.h"

Node::Node(Point<FLOAT> _pos, FLOAT _val, FLOAT _pose_sig, Point<FLOAT> _grad, FLOAT _grad_sig, NODE_TYPE n)
{
    pos = _pos;
    val = _val;
    grad = _grad;
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
    nt = n;
}

Node::Node(Point<FLOAT> _pos, NODE_TYPE _nt): val(0.0),
                                              grad(Point<FLOAT>()),
                                              pose_sig(0.0),
                                              grad_sig(0.0)
{
    pos = _pos;
    nt = _nt;
}

Node::Node():val(0.0), pose_sig(0.0), grad_sig(0.0)
{
    nt = NODE_TYPE::NONE;
}

void Node::updateData(FLOAT _val, FLOAT _pose_sig, Point<FLOAT> _grad, FLOAT _grad_sig, NODE_TYPE n){
    val = _val;
    grad = _grad;
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
    nt = n;
}


void Node::updateNoise( FLOAT _pose_sig, FLOAT _grad_sig){
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
}

//// Node3
Node3::Node3(Point3<FLOAT> _pos, FLOAT _val, FLOAT _pose_sig, Point3<FLOAT> _grad, FLOAT _grad_sig, NODE_TYPE n)
{
    pos = _pos;
    val = _val;
    grad = _grad;
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
    nt = n;
}

Node3::Node3(Point3<FLOAT> _pos, NODE_TYPE _nt): val(0.0),
                                                                grad(Point3<FLOAT>()),
                                                                pose_sig(0.0),
                                                                grad_sig(0.0)
{
    pos = _pos;
    nt = _nt;
}

Node3::Node3():val(0.0), pose_sig(0.0), grad_sig(0.0)
{
    nt = NODE_TYPE::NONE;
}

void Node3::updateData(FLOAT _val, FLOAT _pose_sig, Point3<FLOAT> _grad, FLOAT _grad_sig, NODE_TYPE n){
    val = _val;
    grad = _grad;
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
    nt = n;
}

void Node3::updateNoise( FLOAT _pose_sig, FLOAT _grad_sig){
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
}

