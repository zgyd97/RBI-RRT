#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <list>
#include "kdtree.h"
using namespace std;

struct State {
    vector<double> values;
    State() {}
    State(vector<double> vals) : values(vals) {}
};

struct Obstacle {
    double x; // 左下角 x 坐标
    double y; // 左下角 y 坐标
    double w; // 宽度
    double h; // 高度
};

struct rrtNode {
    State state;
    rrtNode* parent;   // 父节点 index
    list<rrtNode*> children; // 子节点 index
    double cost; // 从根节点到该节点的代价

    rrtNode() {
        this->parent = nullptr;
        this->children = {};
        this->cost = 0.0;
    }
    rrtNode(State state){
        this->state = state;
        this->parent = nullptr;
        this->children = {};
    }
    rrtNode(State state, rrtNode* parent){
        this->state = state;
        this->parent = parent;
        this->children = {};
    }
    rrtNode(State state, rrtNode* parent, list<rrtNode*> children){
        this->state = state;
        this->parent = parent;
        this->children = children;
    }
    rrtNode(State state, rrtNode* parent, list<rrtNode*> children, double cost){
        this->state = state;
        this->parent = parent;
        this->children = children;
        this->cost = cost;
    }
};

struct rrtTree {
    rrtNode* root;
    int num_rrtNodes_num, valid_start_tree_node_nums;
    vector<rrtNode*> tree;
    
    rrtTree()
    {
        root=nullptr;
        num_rrtNodes_num = 1000;

        valid_start_tree_node_nums = 0;
        tree.resize(num_rrtNodes_num); //预分配10000个节点
        
        for(int i=0;i<num_rrtNodes_num;i++)
        {
            tree[i] = new rrtNode();
        }
    }

    rrtTree(int num_rrtNodes)
    {
        root=nullptr;
        num_rrtNodes_num = num_rrtNodes;

        valid_start_tree_node_nums = 0;
        tree.resize(num_rrtNodes_num); //预分配10000个节点
        
        for(int i=0;i<num_rrtNodes_num;i++)
        {
            tree[i] = new rrtNode();
        }
    }

    double distance(const vector<double>& stat_a, const vector<double>& stat_b)
    {
        double sum = 0.0;
        for (size_t i = 0; i < stat_a.size(); ++i)
        {
            double diff = stat_a[i] - stat_b[i];
            sum += diff * diff;
        }
        return sqrt(sum);
    }

    rrtNode* addTreeNode(State state, rrtNode* parent)
    {
        rrtNode* newNode = tree[valid_start_tree_node_nums++];
        if(parent) parent->children.push_back(newNode);
        newNode->cost = parent ? parent->cost + distance(parent->state.values, state.values) : 0.0;
        newNode->state = state;
        newNode->parent = parent;
        return newNode;
    }
};
