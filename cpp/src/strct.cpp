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

