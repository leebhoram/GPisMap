// strct.h
#include <vector>
#include <memory>
#include "strct.h"

Node::Node(Point<float> _pos, float _val, float _pose_sig, Point<float> _grad, float _grad_sig, NODE_TYPE n)
{
    pos = _pos;
    val = _val;
    grad = _grad;
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
    nt = n;
}

Node::Node(Point<float> _pos, NODE_TYPE _nt): val(0.0),
                                              grad(Point<float>()),
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

void Node::updateData(float _val, float _pose_sig, Point<float> _grad, float _grad_sig, NODE_TYPE n){
    val = _val;
    grad = _grad;
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
    nt = n;
}

void Node::updateNoise( float _pose_sig, float _grad_sig){
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
}

Node3::Node3(Point3<float> _pos, float _val, float _pose_sig, Point3<float> _grad, float _grad_sig, NODE_TYPE n)
{
    pos = _pos;
    val = _val;
    grad = _grad;
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
    nt = n;
}

Node3::Node3(Point3<float> _pos, NODE_TYPE _nt): val(0.0),
                                                                grad(Point3<float>()),
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

void Node3::updateData(float _val, float _pose_sig, Point3<float> _grad, float _grad_sig, NODE_TYPE n){
    val = _val;
    grad = _grad;
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
    nt = n;
}

void Node3::updateNoise( float _pose_sig, float _grad_sig){
    pose_sig = _pose_sig;
    grad_sig = _grad_sig;
}

